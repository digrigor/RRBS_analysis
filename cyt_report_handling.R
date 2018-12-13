run_file = '/home/dio01/RRBS/dio01/Results/QMGC_060708091011_Per_sample_stats.txt'

run_metrics = read.csv(run_file, header=TRUE, stringsAsFactors=FALSE)

colnames(run_metrics)=c('run','sample','qc1','qc2','qc3','usable_reads','mapping_efficiency','cpgs', 'path')
run_metrics = run_metrics[,c(1,2,6,7,8,9)]

#Filter coverage 10 000 000
use10 = run_metrics[run_metrics$usable_reads>=10000000,]

#Remove UCL runs
use10 = use10[-grep('^\\w\\d+-2.+', use10$run),]

use10 = use10[order(use10$sample),]


templist = c()
counter=1
for(name in unique(use10$sample)){
    temp = use10[use10$sample==name,]
    templist[counter] = list(subset(temp, usable_reads == max(usable_reads)))
    counter=counter+1
}

templist=unlist(templist)
use10_fin = matrix(templist, byrow=TRUE, ncol=6)
use10_fin = as.data.frame(use10_fin, stringsAsFactors=FALSE)
colnames(use10_fin)=c('run','sample','usable_reads','mapping_efficiency','cpgs', 'path')

file_names = use10_fin$sample
file_paths = use10_fin$path
file_names <- data.frame(lab_no = file_names, path = file_paths)
pheno_all <- read.table("/store/LevelFour/ARTISTIC/METHYLATION/RRBS/mmf27/RELEASE/Data/RRBS_Sample_phenotypes_all.txt", header=T, sep="\t")
pheno_all = pheno_all[, c("lab_no", "Phenotype.RRBS")]
colnames(pheno_all)[2] = "pheno"
pheno_data = merge(file_names, pheno_all)
sampleID = as.list(as.character(pheno_data$lab_no))
samplePATH = as.list(as.character(file_paths))
treatment = as.numeric(pheno_data$pheno =="Case")

myobj=methRead(samplePATH, sample.id=sampleID, assembly="hg19", pipeline='bismarkCytosineReport', treatment=treatment, context="CpG")
filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL)
meth.min_5=unite(filtered.myobj,min.per.group=5L, mc.cores=4)
meth.min_5_df <- getData(meth.min_5)
meth.min_5_per = percMethylation(meth.min_5)
myDiff_chisq=calculateDiffMeth(meth.min_5, mc.cores=4)
myDiff_fish=calculateDiffMeth(meth.min_5, overdispersion='MN', test="F", mc.cores=4)
myDiff_3=calculateDiffMeth(meth.min_5, overdispersion='MN', test="Chisq", mc.cores=4)
diff_chi_df = getData(myDiff_chisq)
diff_fish_df = getData(myDiff_fish)
diff_chi2_df = getData(myDiff_3)

conts = cases = means_cont = means_case = medians_cont = medians_case = c()
cases = c()
for(i in c(1:nrow(diff_chi_df))){
    current_cont = sum(!is.na(meth.min_5_per[i, c(which(treatment==0))]))
    current_case = sum(!is.na(meth.min_5_per[i, c(which(treatment==1))]))
    mean_cont = mean(meth.min_5_per[i,which(treatment==0)], na.rm=TRUE)
    median_cont = median(meth.min_5_per[i,which(treatment==0)], na.rm=TRUE)
    mean_case = mean(meth.min_5_per[i,which(treatment==1)], na.rm=TRUE)
    median_case = median(meth.min_5_per[i,which(treatment==1)], na.rm=TRUE)
    conts[i]=current_cont
    cases[i]=current_case
    means_cont[i]= mean_cont
    means_case[i] = mean_case
    medians_cont[i] = median_cont
    medians_case[i] = median_case
}

main_diff = cbind(diff_chi_df[,c(1,2,3,4,5,6)],diff_fish_df[,c(5,6)], diff_chi2_df[,c(5,6,7)],
            means_cont, means_case, means_case-means_cont, medians_cont, medians_case)

colnames(main_diff)[c(5,6,7,8,9,10,12,13,14,15,16,17,18)] = c('chi_pvalue','chi_qvalue', 'mnfisher_pvalue','mnfisher_qvalue',
                                                            'mnchi_pvalue','mnchi_qvalue','mean_control', 'mean_case', 'mean_difference',
                                                            'median_control', 'median_case', 'covered_controls', 'covered_cases')

main_diff = main_diff[order(-abs(main_diff$mean_difference)),]

#Try to cluster the samples
dists = c("correlation", "euclidean", "maximum", "manhattan", "canberra", "minkowski")
meths = c("ward", "single", "complete", "average", "mcquitty", "median","centroid")
for(i in dists){
for(j in meths){
jpeg(paste0(i,j,'.jpg'))
try(clusterSamples(meth.min_5[as.numeric(rownames(subset(main_diff, mnchi_qvalue < 0.05 & covered_controls>=13 & covered_cases >=13 & abs(mean_difference) > 18))),], dist=i, method=j, plot=TRUE))
dev.off()
}
}

#Get a top list
top_dm <- subset(main_diff, mnchi_pvalue < 0.05 & covered_controls>=8 & covered_cases >=8 & abs(mean_difference) > 15)
top_dm$cpg_id <- paste0(top_dm$chr, top_dm$start)
top_for_gr <- top_dm[,c(1,2,3,4,9,10,14,17,18,19)]
top_gr <- makeGRangesFromDataFrame(top_for_gr, keep.extra.columns = TRUE)

#Do the annotation
dm_annotated = annotate_regions(
    regions = top_gr,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
df_dm_annotated = data.frame(dm_annotated)

#Transpose it to good format
tr_df_df_annotated <- transpose_annotation(df_dm_annotated)



png(filename="../results/chi2_pvolcano_plot.png", units="in", width=11.5, height=8, res=600)
par(mar=c(5.1,6,4.1,2.1))
plot(main_diff$mean_difference, -log10(main_diff$mnchi_pvalue), xlab = "Mean Difference of Methylation Percentage", ylab = "-log10 pvalue", pch=19, cex = 0.5, cex.lab=2.5, cex.axis=2, xaxt = "n")
axis(1, at=seq((-ceiling(abs(min(main_diff$mean_difference)))), ceiling(max(main_diff$mean_difference)), by=0.5), cex.axis=2)
abline(v=-18, col="black", lty=2)
abline(v=-25, col="black", lty=2)
abline(v=18, col="black", lty=2)
abline(v=25, col="black", lty=2)
abline(a=-log10(0.05), b=0, col="blue")
#legend("bottomright", legend=c("logCPM2 > -5","logCPM2 < -5"), col = c("red", "black"), pch=19, cex=1.2)
dev.off()

boxplot_cpgs <- function(maindiff, cpg, gene){
    png(filename=paste0("../results/",gene,".png"), units="in", width=5, height=5, res=600)
    row = main_diff[main_diff$cpg_id==cpg,]
    name = rownames(row)
    pv = format(round(row$mnchi_pvalue, 3), nsmall = 3)
    x = meth.min_5_per[name, which(treatment==0)]
    conts = sum(!is.na(x))
    y = meth.min_5_per[name, which(treatment==1)]
    cases = sum(!is.na(y))
    boxplot(x, y,
            main=gene,
            ylab = "Methylation Percentage (%)",
            ylim = c(0,100),
            xlab = "Phenotype",
            col = c("blue","purple"),
            names = c(paste0("Controls (N=",conts,")"), paste0("Cases (N=",cases,")"))
    )
    text(x= 1.5, y= -1.5, labels=paste0("p=",pv), cex=0.8)
    dev.off()
}

cpgsbox = read.csv('cpgstobox.csv', header=FALSE, stringsAsFactors=FALSE)
for(i in c(1:nrow(cpgsbox))){
    boxplot_cpgs(main_diff, cpgsbox[i,1], cpgsbox[i,2])
}