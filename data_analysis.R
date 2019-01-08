library(methylKit)

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


