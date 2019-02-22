#library(doMC)
#registerDoMC(30)
library(randomForest)
library(mlbench)
library(caret)
library(methylKit)
library(rlist)
library(GenomicRanges)

source("data_analysis_functions.R")

ucsc_to_hgnc = read.csv('~/Methyl_Annotations/ucsc_gene_id_to_symbol.txt', sep='\t', row.names=1)

master_csv = "../../Methylkit_DM_analysis/Input/Masterfile_16.01.2019.csv"
master = read.csv(master_csv)
master_exome_csv = "../../Exome_Methylkit_DM_analysis"
metrics = readRDS("../../Methylkit_DM_analysis/Input/all_run_metrics.rds")

master10=filter_out_samples(metrics,master,7,function(x)x<100000)

globalmaster = merge(x=master10, y=metrics, by.x="lab_no", by.y="sample")
globalmaster = globalmaster[order(globalmaster$RRBS.Set_ID),]
# globalmaster = globalmaster[!globalmaster$lab_no == 'H18127',]
conveff = read.csv('/store/LevelFour/ARTISTIC/METHYLATION/RRBS/dio01/Spike_in_Controls_analysis/Results/Bedgraph/CpG/conversion_efficiency.csv')
colnames(conveff)[1]='run'

globalmaster2 = merge(globalmaster, conveff, by.x='run', by.y='run')
globalmasterOLD = globalmaster
globalmaster2 = globalmaster2[!globalmaster2$RRBS.Set_ID == globalmaster2[which(globalmaster2$Unethylated_control_CR<98),]$RRBS.Set_ID,] #Exclude samples <98% conversion efficiency and their matches

sample_ages = read.csv('../../Methylkit_DM_analysis/Input/Age.csv', header=FALSE)
age_matrix = data.frame(age=sample_ages$V2)
rownames(age_matrix)= sample_ages$V1
colnames(sample_ages) = c('lab_no', 'age')

globalmaster3 = merge(globalmaster2, sample_ages, by.x='lab_no', by.y='lab_no')
globalmaster = globalmaster3

pheno_matrix = data.frame(treat=globalmaster$Phenotype.RRBS)
rownames(pheno_matrix) = globalmaster$lab_no
hpv_matrix = data.frame(treat=paste0('hpv', globalmaster$HPV.type))
rownames(hpv_matrix) = globalmaster$lab_no

#diff_regs = readRDS('../../Output/class_30_01_2019/diff_regs.rds')
#meth_win_per = readRDS('../../Output/class_30_01_2019/meth_win_per.rds')


#Reads VS number of CpGs with more than 9 reads plot
areads = globalmaster$usable_reads*globalmaster$mapping_efficiency*0.01
par(mar=c(5.1,6.1,4.1,2.1))
plot(areads, globalmaster$cpgs10, xlab='Aligned Reads', ylab='CpGs covered by at least 10 reads', cex.lab=1.5, pch=19)

#Preparing methylkit's input
pheno = globalmaster$Phenotype.RRBS
sampleID = as.list(as.character(globalmaster$lab_no))
samplePATH = as.list(as.character(globalmaster$path))
treatment = as.numeric(pheno =="Case")
age = factor(globalmaster$age)

#Preparing methylkit's input
myobj=methRead(samplePATH, sample.id=sampleID, assembly="hg19", pipeline='bismarkCytosineReport', treatment=treatment, context="CpG", mincov=2)
#myobj.db = makeMethylDB(myobj)
filtered.myobj_2=filterByCoverage(myobj,lo.count=2,lo.perc=NULL,hi.count=NULL,hi.perc=99.5)
filtered.myobj_10=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99.5)
#saveRDS(myobj, "../../Output/myobj.rds")
#saveRDS(filtered.myobj, "../../Output/filtered.myobj.rds")
meth.min_5=unite(filtered.myobj,min.per.group=5L)
myDiff_DSS = calculateDiffMethDSS(meth.min_5, mc.cores = 6)
#saveRDS(meth.min_5,  "../../Output/meth.min_5.rds")
myDiff4 = calculateDiffMeth(meth.min_5, overdispersion = "none", effect="wmean", test = "F", mc.cores = 6)
myDiff5 = calculateDiffMeth(meth.min_5, overdispersion = "none", effect="wmean", test = "Chisq", mc.cores = 6)
myDiff6 = calculateDiffMeth(meth.min_5, overdispersion = "MN", effect="wmean", test = "Chisq", mc.cores = 1)

regs = regions_from_cpgs(diff_df=topdiff, dmrs=dmr_table, size=100)
regs = data.frame(regs, stringsAsFactors=FALSE)
combs_regs = best_comb(regs, 85)

myDiff_DSS = calculateDiffMethDSS(meth.min_5, mc.cores = 6)
meth.min_5_df <- getData(meth.min_5)
meth.min_5_per = percMethylation(meth.min_5)
myDiff_DSS_df = getData(myDiff_DSS)

maindiff = get_final_diff_matrix(diffs)

rank_wcox = get_the_rank <- function(tablez, 12)
rank_dss = get_the_rank <- function(tablez, 14)

min20_bed <- maindiff_mincontscases_20[,c(1,2,3,14,6)]
colnames(min20_bed) = c('chrom', 'start', 'end', 'pvalue', 't')
write.table(min20_bed, '../../Output/min20.bed', quote=FALSE, row.names=FALSE, sep='\t')

top_diff = subset(maindiff, conts>=20 & cases >= 20)
top_diff$cpg_id <- paste0(top_diff$chr, top_diff$start)
top_gr <- makeGRangesFromDataFrame(top_diff, keep.extra.columns = TRUE)

p_dm_annotated = annotate_regions(
    regions = p_diff_gr,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
df_p_dm_annotated = data.frame(p_dm_annotated, stringsAsFactors=FALSE)


myDiff4 = calculateDiffMeth(meth.min_5, overdispersion = "none", effect="wmean", test = "F", mc.cores = 6)
myDiff5 = calculateDiffMeth(meth.min_5, overdispersion = "none", effect="wmean", test = "Chisq", mc.cores = 6)
myDiff6 = calculateDiffMeth(meth.min_5, overdispersion = "MN", effect="wmean", test = "Chisq", mc.cores = 6)
saveRDS(myDiff4, "../../../Output/mydiff4.rds")
saveRDS(myDiff5, "../../../Output/mydiff5.rds")
saveRDS(myDiff6, "../../../Output/mydiff6.rds")

get_plot <- function(sample){

}
#Get a list with maximum effect sizes per gene
fin_list=c()
df_dm_annotated$gene_symbol = as.character(ucsc_to_hgnc[df_dm_annotated$annot.tx_id,])
unique_genes = unique(df_dm_annotated$gene_symbol)
for(ge in unique_genes){
    rows <- df_dm_annotated[which(df_dm_annotated$gene_symbol == ge),]
    fin_list = c(fin_list, rows[which.max(abs(rows$mediansdiff)),c('meansdiff','mediansdiff','wcoxs')])
}

library(Gviz)
library(rtracklayer)
library(trackViewer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
#draw lollipops
#EPB41L3
gr <- GRanges("chr18", IRanges(5392379, 5635642), strand="*")
gr_c <- c("chr18", 5392379, 5635642, strand="-")
to_lol = draw_cpgs(maindiff, gr=gr_c)
to_lol[[2]][start(to_lol[[2]]) %in% c(5543548, 5543559, 5543561),]$score = 0.01
names(to_lol[[2]]) = paste0('Controls N=',maindiff[names(to_lol[[2]]), 'conts'],' | Cases N=',maindiff[names(to_lol[[2]]), 'cases'])
leg_cols = c('#e6ffe6', '#269900', '#ffe6e6', '#ff0000', 'grey')
legend <- list(labels=c('Hyper p>0.05', 'Hyper p<0.05', 'Hypo p>0.05', 'Hypo p<0.05', 'Unchanged'), col=leg_cols, fill=leg_cols)
lolliplot(to_lol[[2]][to_lol[[2]]$score>0,], to_lol[[1]], type="circle", ranges=gr, legend=legend, rescale=TRUE, cex=.5)

#CADM1
gr <- GRanges("chr11", IRanges(115044344, 115380241), strand="*")
gr_c <- c('chr11', 115044344, 115380241, strand="-")
to_lol = draw_cpgs(maindiff, gr=gr_c)
#to_lol[[2]][start(to_lol[[2]]) %in% c(5543548, 5543559, 5543561),]$score=0.01
leg_cols = c('#e6ffe6', '#269900', '#ffe6e6')
names(to_lol[[2]])=NULL
legend <- list(labels=c('Hyper p>0.05', 'Hyper p<0.05', 'Hypo p>0.05', 'Hypo p<0.05'), col=leg_cols, fill=leg_cols)
lolliplot(to_lol[[2]][to_lol[[2]]$score>0,], to_lol[[1]], type="circle", ranges=gr, legend=legend, rescale=TRUE, cex=.5)

#FAM19A4
gr <- GRanges("chr3", IRanges(68780914, 68986761), strand="*")
gr_c <- c('chr3', 68780914, 68986761, strand="-")
to_lol = draw_cpgs(maindiff, gr=gr_c)
names(to_lol[[2]])=NULL
leg_cols = c('#e6ffe6', '#269900', '#ffe6e6')
legend <- list(labels=c('Hyper p>0.05', 'Hyper p<0.05', 'Hypo p>0.05', 'Hypo p<0.05'), col=leg_cols, fill=leg_cols)
lolliplot(to_lol[[2]][to_lol[[2]]$score>0,], to_lol[[1]], type="circle", ranges=gr, legend=legend, rescale=TRUE, cex=.5)

#POU4F3
gr <- GRanges("chr5", IRanges(145713586, 145720083), strand="*")
gr_c <- c('chr5', 145713586, 145720083, strand="+")
to_lol = draw_cpgs(maindiff, gr=gr_c)
names(to_lol[[2]])=NULL
leg_cols = c('#e6ffe6', '#269900', '#ffe6e6')
legend <- list(labels=c('Hyper p>0.05', 'Hyper p<0.05', 'Hypo p>0.05', 'Hypo p<0.05'), col=leg_cols, fill=leg_cols)
lolliplot(to_lol[[2]][to_lol[[2]]$score>0,], to_lol[[1]], type="circle", ranges=gr, legend=legend, rescale=TRUE, cex=.5)

#MAL
gr <- GRanges("chr2", IRanges(95686399, 95719737), strand="*")
gr_c <- c('chr2', 95686399, 95719737, strand="+")
to_lol = draw_cpgs(maindiff, gr=gr_c)
names(to_lol[[2]])=NULL
leg_cols = c('#e6ffe6', '#269900', '#ffe6e6')
legend <- list(labels=c('Hyper p>0.05', 'Hyper p<0.05', 'Hypo p>0.05', 'Hypo p<0.05'), col=leg_cols, fill=leg_cols)
lolliplot(to_lol[[2]][to_lol[[2]]$score>0,], to_lol[[1]], type="circle", ranges=gr, legend=legend, rescale=TRUE, cex=.5)

#PAX1
gr <- GRanges("chr20", IRanges(21681296, 21699124), strand="*")
gr_c <- c('chr20', 21681296, 21699124, strand="+")
to_lol = draw_cpgs(maindiff, gr=gr_c)
names(to_lol[[2]])=NULL
leg_cols = c('#e6ffe6', '#269900', '#ffe6e6')
legend <- list(labels=c('Hyper p>0.05', 'Hyper p<0.05', 'Hypo p>0.05', 'Hypo p<0.05'), col=leg_cols, fill=leg_cols)
lolliplot(to_lol[[2]][to_lol[[2]]$score>0,], to_lol[[1]], type="circle", ranges=gr, legend=legend, rescale=TRUE, cex=.5)

#SOX1
gr <- GRanges("chr13", IRanges(112721912, 112726020), strand="*")
gr_c <- c('chr13', 112721912, 112726020, strand="+")
to_lol = draw_cpgs(maindiff, gr=gr_c)
names(to_lol[[2]])=NULL
leg_cols = c('#e6ffe6', '#269900', '#ffe6e6')
legend <- list(labels=c('Hyper p>0.05', 'Hyper p<0.05', 'Hypo p>0.05', 'Hypo p<0.05'), col=leg_cols, fill=leg_cols)
lolliplot(to_lol[[2]][to_lol[[2]]$score>0,], to_lol[[1]], type="circle", ranges=gr, legend=legend, rescale=TRUE, cex=.5)

#MIR124-2
gr <- GRanges("chr8", IRanges(65286705, 65291814), strand="*")
gr_c <- c('chr8', 65286705, 65291814, strand="+")
to_lol = draw_cpgs(maindiff, gr=gr_c)
names(to_lol[[2]])=NULL
leg_cols = c('#e6ffe6', '#269900', '#ffe6e6')
legend <- list(labels=c('Hyper p>0.05', 'Hyper p<0.05', 'Hypo p>0.05', 'Hypo p<0.05'), col=leg_cols, fill=leg_cols)
lolliplot(to_lol[[2]][to_lol[[2]]$score>0,], to_lol[[1]], type="circle", ranges=gr, legend=legend, rescale=TRUE, cex=.5)

#ASCL1
gr <- GRanges("chr12", IRanges(103346452, 103354294), strand="*")
gr_c <- c('chr12', 103346452, 103354294, strand="+")
to_lol = draw_cpgs(maindiff, gr=gr_c)
names(to_lol[[2]])=NULL
leg_cols = c('#e6ffe6', '#269900', '#ffe6e6')
legend <- list(labels=c('Hyper p>0.05', 'Hyper p<0.05', 'Hypo p>0.05', 'Hypo p<0.05'), col=leg_cols, fill=leg_cols)
lolliplot(to_lol[[2]][to_lol[[2]]$score>0,], to_lol[[1]], type="circle", ranges=gr, legend=legend, rescale=TRUE, cex=.5)

#LHX8
gr <- GRanges("chr1", IRanges(75589119, 75627218), strand="*")
gr_c <- c('chr1', 75589119, 75627218, strand="+")
to_lol = draw_cpgs(maindiff, gr=gr_c)
names(to_lol[[2]]) = paste0('Controls N=',maindiff[names(to_lol[[2]]), 'conts'],' | Cases N=',maindiff[names(to_lol[[2]]), 'cases'])
leg_cols = c('#e6ffe6', '#269900', '#ffe6e6')
legend <- list(labels=c('Hyper p>0.05', 'Hyper p<0.05', 'Hypo p>0.05', 'Hypo p<0.05'), col=leg_cols, fill=leg_cols)
lolliplot(to_lol[[2]][to_lol[[2]]$score>0,], to_lol[[1]], type="circle", ranges=gr, legend=legend, rescale=TRUE, cex=.5)

#ST6GALNAC5
gr <- GRanges("chr1", IRanges(77328165, 77533231), strand="*")
gr_c <- c('chr1', 77328165, 77533231, strand="+")
to_lol = draw_cpgs(maindiff, gr=gr_c)
names(to_lol[[2]])=NULL
leg_cols = c('#e6ffe6', '#269900', '#ffe6e6')
legend <- list(labels=c('Hyper p>0.05', 'Hyper p<0.05', 'Hypo p>0.05', 'Hypo p<0.05'), col=leg_cols, fill=leg_cols)
lolliplot(to_lol[[2]][to_lol[[2]]$score>0,], to_lol[[1]], type="circle", ranges=gr, legend=legend, rescale=TRUE, cex=.5)
