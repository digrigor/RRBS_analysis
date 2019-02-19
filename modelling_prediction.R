library(doMC)
library(methylKit)
registerDoMC(50)
library(randomForest)
library(mlbench)
library(caret)

maindiff = readRDS('../../Methylkit_DM_analysis/Output/RDS_files/maindiff.rds')
filtered.myobj_2 = readRDS('../../Methylkit_DM_analysis/Output/RDS_files/filtered.myobj_.rds')
dmr_table = read.table('../../pcomb_DMR_analysis/Output/wcox.anno.hg19.bed')
#diff_regs = readRDS('../../Methylkit_DM_analysis/Output/RDS_files/diff_regs.rds')
#meth_win_per = readRDS('../../Methylkit_DM_analysis/Output/RDS_files/meth_win_per.rds')

diff_regs = regions_from_cpgs(maindiff, dmr_table, filtered.myobj_2, 100, 20, 20, 7, 0.05, 5)

st = rownames(subset(diff_regs, nsamples>200 & abs(win_mediansdiff)>=13))
samples_lims = c(160, 170)
cand_shared_range = c(10,20,30,0)
all_combs = expand.grid(st, samples_lims, cand_shared_range)

fin_list = list()
  for(i in c(1:nrow(all_combs))){
      print('########################################################')
      print(paste0(i, ' out of ', nrow(all_combs)))
      print(paste0('Starting point: ',all_combs[i,1]))
      print(paste0('Minimum shared: ', all_combs[i,2]))
      print(paste0('Close range: ', all_combs[i,3]))
      fin_list[[i]]=c(feature_selection(diff_regs, starting_point=as.character(all_combs[i,1]), shared_samples_limit=as.numeric(all_combs[i,2]), cand_shared_range=as.numeric(all_combs[i,3])))
  }