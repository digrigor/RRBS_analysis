library(methylKit)
source("data_analysis_functions.R")

master_csv = "../../../Input/Masterfile_20.12.2018.csv"
master = read.csv(master_csv)
metrics = readRDS("../../../Input/all_run_metrics.rds")

master10=filter_out_samples(metrics,master,7,function(x)x<100000)

globalmaster = merge(x=master10, y=metrics, by.x="lab_no", by.y="sample")
globalmaster = globalmaster[order(globalmaster$RRBS.Set_ID),]

pheno = globalmaster$Phenotype.RRBS
sampleID = as.list(as.character(globalmaster$lab_no))
samplePATH = as.list(as.character(globalmaster$path))
treatment = as.numeric(pheno =="Case")

myobj=methRead(samplePATH, sample.id=sampleID, assembly="hg19", pipeline='bismarkCytosineReport', treatment=treatment, context="CpG")
saveRDS(myobj, "/myobj.rds")
filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL)
saveRDS(filtered.myobj, "/filtered.myobj.rds")
meth.min_5=unite(filtered.myobj,min.per.group=5L, mc.cores=12)
saveRDS(meth.min_5,  "../../../Output/meth.min_5.rds")

meth.min_5_df <- getData(meth.min_5)
meth.min_5_per = percMethylation(meth.min_5)

myDiff4 = calculateDiffMeth(meth.min_5, overdispersion = "none", effect="wmean", test = "F", mc.cores = 6)
myDiff5 = calculateDiffMeth(meth.min_5, overdispersion = "none", effect="wmean", test = "Chisq", mc.cores = 6)
myDiff6 = calculateDiffMeth(meth.min_5, overdispersion = "MN", effect="wmean", test = "Chisq", mc.cores = 6)
saveRDS(myDiff4, "../../../Output/mydiff4.rds")
saveRDS(myDiff5, "../../../Output/mydiff5.rds")
saveRDS(myDiff6, "../../../Output/mydiff6.rds")



