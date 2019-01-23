library(methylKit)
source("data_analysis_functions.R")

master_csv = "../../Input/Masterfile_16.01.2019.csv"
master = read.csv(master_csv)
metrics = readRDS("../../Input/all_run_metrics.rds")

master10=filter_out_samples(metrics,master,7,function(x)x<100000)

globalmaster = merge(x=master10, y=metrics, by.x="lab_no", by.y="sample")
globalmaster = globalmaster[order(globalmaster$RRBS.Set_ID),]
globalmaster = globalmaster[!globalmaster$lab_no == 'H18127',]

#Reads VS number of CpGs with more than 9 reads plot
areads = globalmaster$usable_reads*globalmaster$mapping_efficiency*0.01
par(mar=c(5.1,6.1,4.1,2.1))
plot(areads, globalmaster$cpgs10, xlab='Aligned Reads', ylab='CpGs covered by at least 10 reads', cex.lab=1.5, pch=19)

#Preparing methylkit's input
pheno = globalmaster$Phenotype.RRBS
sampleID = as.list(as.character(globalmaster$lab_no))
samplePATH = as.list(as.character(globalmaster$path))
treatment = as.numeric(pheno =="Case")

#Preparing methylkit's input
myobj=methRead(samplePATH, sample.id=sampleID, assembly="hg19", pipeline='bismarkCytosineReport', treatment=treatment, context="CpG")
myobj.db = makeMethylDB(myobj)
filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99.9)
#saveRDS(myobj, "../../Output/myobj.rds")
#saveRDS(filtered.myobj, "../../Output/filtered.myobj.rds")
meth.min_5=unite(filtered.myobj,min.per.group=5L)
#saveRDS(meth.min_5,  "../../Output/meth.min_5.rds")
myDiff4 = calculateDiffMeth(meth.min_5, overdispersion = "none", effect="wmean", test = "F", mc.cores = 6)
myDiff5 = calculateDiffMeth(meth.min_5, overdispersion = "none", effect="wmean", test = "Chisq", mc.cores = 6)
myDiff6 = calculateDiffMeth(meth.min_5, overdispersion = "MN", effect="wmean", test = "Chisq", mc.cores = 1)
myDiff7 = calculateDiffMethDSS(meth.min_5)




meth.min_5_df <- getData(meth.min_5)
meth.min_5_per = percMethylation(meth.min_5)

myDiff4 = calculateDiffMeth(meth.min_5, overdispersion = "none", effect="wmean", test = "F", mc.cores = 6)
myDiff5 = calculateDiffMeth(meth.min_5, overdispersion = "none", effect="wmean", test = "Chisq", mc.cores = 6)
myDiff6 = calculateDiffMeth(meth.min_5, overdispersion = "MN", effect="wmean", test = "Chisq", mc.cores = 6)
saveRDS(myDiff4, "../../../Output/mydiff4.rds")
saveRDS(myDiff5, "../../../Output/mydiff5.rds")
saveRDS(myDiff6, "../../../Output/mydiff6.rds")

get_plot <- function(sample){

}



