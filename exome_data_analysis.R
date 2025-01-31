library(methylKit)
#library(doMC)
#registerDoMC(30)
library(randomForest)
library(mlbench)
library(caret)
source("data_analysis_functions.R")

ucsc_to_hgnc = read.csv('~/Methyl_Annotations/ucsc_gene_id_to_symbol.txt', sep='\t', row.names=1)

master_exome_csv = "../../Exome_Methylkit_DM_analysis/Input/MASTERFILE_10.12.2018.csv"
master_exome = read.csv(master_exome_csv)
master_csv = "../../Methylkit_DM_analysis/Input/Masterfile_16.01.2019.csv"
master = read.csv(master_csv)
metrics = readRDS("../../Methylkit_DM_analysis/Input/all_run_metrics.rds")

metrics = metrics[metrics$sample %in% as.character(master_exome$Lab.no),]
metrics = metrics[metrics$cpgs10>2000,]
#metrics = metrics[!metrics$sample=="H02551",]

globalmaster = merge(x=master_exome, y=metrics, by.x="Lab.no", by.y="sample")

# globalmaster = globalmaster[order(globalmaster$RRBS.Set_ID),]
globalmaster = globalmaster[!globalmaster$Lab.no == 'H18127',]
globalmaster = globalmaster[!globalmaster$Lab.no == 'H02551',]
#conveff = read.csv('/store/LevelFour/ARTISTIC/METHYLATION/RRBS/dio01/Spike_in_Controls_analysis/Results/Bedgraph/CpG/conversion_efficiency.csv')
#colnames(conveff)[1]='run'

#globalmaster = merge(globalmaster, conveff, by.x='run', by.y='run')

sample_ages = read.csv('../../Methylkit_DM_analysis/Input/Age.csv', header=FALSE)
colnames(sample_ages)=c('lab_no','age')
age_matrix = data.frame(age=sample_ages$age)
rownames(age_matrix)=sample_ages$lab_no


globalmaster2 = merge(globalmaster, sample_ages, by.x='Lab.no', by.y='lab_no')
globalmaster = globalmaster2

pheno_table = read.csv('../../Methylkit_DM_analysis/Input/RRBS_Sample_phenotypes_all.txt', sep='\t')
pheno_matrix = data.frame(pheno = pheno_table$Phenotype.RRBS)
rownames(pheno_matrix) = pheno_table$lab_no

#diff_regs = readRDS('../../Output/class_30_01_2019/diff_regs.rds')
#meth_win_per = readRDS('../../Output/class_30_01_2019/meth_win_per.rds')


#Reads VS number of CpGs with more than 9 reads plot
areads = globalmaster$usable_reads*globalmaster$mapping_efficiency*0.01
par(mar=c(5.1,6.1,4.1,2.1))
plot(areads, globalmaster$cpgs10, xlab='Aligned Reads', ylab='CpGs covered by at least 10 reads', cex.lab=1.5, pch=19)

#Preparing methylkit's input
pheno = globalmaster$Phenotype.RRBS
sampleID = as.list(as.character(globalmaster$Lab.no))
samplePATH = as.list(as.character(globalmaster$path))
treatment = as.numeric(pheno =="Case")
age = factor(globalmaster$age)

#Preparing methylkit's input
myobj=methRead(samplePATH, sample.id=sampleID, assembly="hg19", pipeline='bismarkCytosineReport', treatment=treatment, context="CpG", mincov=0)
#myobj.db = makeMethylDB(myobj)
#filtered.myobj_5=filterByCoverage(myobj,lo.count=5,lo.perc=NULL,hi.count=NULL,hi.perc=99.5)
filtered.myobj_10=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99.5)
#saveRDS(myobj, "../../Output/myobj.rds")
#saveRDS(filtered.myobj, "../../Output/filtered.myobj.rds")
meth.min_5=unite(filtered.myobj,min.per.group=5L)
#myDiff_DSS = calculateDiffMethDSS(meth.min_5, mc.cores = 6)
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

top_diff = subset(maindiff, conts>=20 & cases >= 20 & abs(meansdiff) >= 20 & wcox_pvalue < 0.05)
top_diff$cpg_id <- paste0(top_diff$chr, top_diff$start)
top_gr <- makeGRangesFromDataFrame(top_diff, keep.extra.columns = TRUE)

dm_annotated = annotate_regions(
    regions = top_gr,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
df_dm_annotated = data.frame(dm_annotated, stringsAsFactors=FALSE)


myDiff4 = calculateDiffMeth(meth.min_5, overdispersion = "none", effect="wmean", test = "F", mc.cores = 6)
myDiff5 = calculateDiffMeth(meth.min_5, overdispersion = "none", effect="wmean", test = "Chisq", mc.cores = 6)
myDiff6 = calculateDiffMeth(meth.min_5, overdispersion = "MN", effect="wmean", test = "Chisq", mc.cores = 6)
saveRDS(myDiff4, "../../../Output/mydiff4.rds")
saveRDS(myDiff5, "../../../Output/mydiff5.rds")
saveRDS(myDiff6, "../../../Output/mydiff6.rds")

get_plot <- function(sample){

}



