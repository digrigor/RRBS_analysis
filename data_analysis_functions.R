#' Filter out samples that meet a condition and their case/control matches and get a new masterfile without them.
#'
#' @param metrics Dataframe. It is the metrics dataframe containing run information for each sample.
#' @param master Dataframe. It is the masterfile dataframe containing phenotype and other information for each sample.
#' @param metrics_col Integer. On which column of the metrics df should the filter be applied to.
#' @condition A function. The condition that the filtered out samples have to meet.
# #' @examples
#' filter_out_samples(metrics,master,7,function(x)x<100000) to filter out samples with less than 100000 cpgs with at least 10 reads
filter_out_samples <- function(metrics, master, metrics_col, condition){
    a <- metrics[,metrics_col]
    out_ind = which(condition(a))
    filtered_out_samples = metrics[out_ind,]$sample
    out_master_matched_ids = unique(master[master$lab_no %in% filtered_out_samples,]$RRBS.Set_ID)
    return(master[!master$RRBS.Set_ID %in% out_master_matched_ids,])
}

#' Create a stranded GRange object of all the CpG sites of the hg19 human genome based on the BSgenome.Hsapiens.UCSC.hg19 package.
#'
#' @param output String. Where the produced table will be stored. It should end with '/'
#' @param save Boolean. Whether to store or not the output files. Output files: RDS object of the CpGs dataframe, RDS object of the CpGs GRanges object and a text file of all CpGs.
# #' @examples
#' get_all_hg19_cpgs('../../../Output/',save=TRUE)
get_all_hg19_cpgs <- function(output, save=TRUE){
    library(BSgenome.Hsapiens.UCSC.hg19)
    library(GenomicRanges)
    chrs <- names(Hsapiens)[1:24]
    cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
    cpgr <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 2))))
    hg19_cpg_df <- as.data.frame(cpgr)
    reps = rep(hg19_cpg_df$start, each=2)
    inds = c(1:length(reps))
    reps[inds[c(FALSE,TRUE)]] = reps[inds[c(FALSE,TRUE)]]+1
    new_starts = reps
    strands = rep(c("+","-"),length(reps)/2)
    old_chrs = as.character(hg19_cpg_df$seqnames)
    chrs = rep(old_chrs, each=2)
    hg19_cpgs = data.frame(seqnames=chrs, start=new_starts, strand=strands, id=paste0(chrs,new_starts))
    hg19_cpgs_gr = makeGRangesFromDataFrame(hg19_cpgs, keep.extra.columns=TRUE, seqnames.field=c('seqnames'), start.field=c('start'), end.field=c('start'), strand.field=c('strand'))
    if(save==TRUE){
        write.table(hg19_cpgs, paste0(output,'hg19_cpgs.txt'), quote=FALSE, row.names=FALSE)
        saveRDS(hg19_cpgs, paste0(output,'hg19_cpgs.rds'))
        saveRDS(hg19_cpgs_gr, paste0(output,'hg19_cpgs_grange.rds'))
    }
    return(hg19_cpgs_gr)
}


#' Summarise all performed methylkit test (methylDiff objects) into a final dataframe also contatining information about the number of cases/controls covering each CpG, wilcoxon test pvalue and mean/median differences.
#'
#' @param mydiffs List of objects. A list of all methylkit methylDiff objects which will be included.
# #' @examples
#'

get_final_diff_matrix <- function(mydiffs=list()){
    library(coin)
    ps_and_es = list()
    for(i in c(1:length(mydiffs))){
        ps_and_es[[i]] = cbind(mydiffs[[i]]$pvalue, mydiffs[[i]]$qvalue)
    }
    ps_and_es_matrix = do.call(cbind, ps_and_es)
    conts = cases = means_cont = means_case = medians_cont = medians_case = meansdiff = mediansdiff = age_meansdiff = wcoxs = age_mediansdiff = hpv_kstest = perms = c()
    z=1
    main_mydiff = getData(mydiffs[[1]])
    #length(conts) = length(cases) = length(means_cont) = length(means_case) = length(medians_cont) = length(medians_case) = length(meansdiff) = length(mediansdiff) = length(wcoxs) = length(perms) = nrow(main_mydiff)
    for(i in as.numeric(rownames(main_mydiff))){
        cont_vals = as.numeric(na.omit(meth.min_5_per[i, c(which(treatment==0))]))
        case_vals = as.numeric(na.omit(meth.min_5_per[i, c(which(treatment==1))]))
        current_cont = sum(!is.na(meth.min_5_per[i, c(which(treatment==0))]))
        current_case = sum(!is.na(meth.min_5_per[i, c(which(treatment==1))]))
        age_conts = age_matrix[names(which(!is.na(meth.min_5_per[i, c(which(treatment==0))]))),]
        age_cases = age_matrix[names(which(!is.na(meth.min_5_per[i, c(which(treatment==1))]))),]
        conts[z]=current_cont
        cases[z]=current_case
        means_cont[z]= mean(cont_vals)
        means_case[z] = mean(case_vals)
        medians_cont[z] = median(cont_vals)
        medians_case[z] = median(case_vals)
        meansdiff[z] = means_case[z] - means_cont[z]
        mediansdiff[z] = medians_case[z] - medians_cont[z]
        age_meansdiff[z] = mean(age_cases)-mean(age_conts)
        age_mediansdiff[z] = median(age_cases)-median(age_conts)
        hpv_kstest[z] = KLD(table(hpv_matrix[names(which(!is.na(meth.min_5_per[i, c(which(treatment==1))]))),]), table(hpv_matrix[names(which(!is.na(meth.min_5_per[i, c(which(treatment==0))]))),]))$sum.KLD.py.px
        wcoxs[z] = wilcox.test(cont_vals, case_vals)$p.value
        vals = c(cont_vals, case_vals)
        treats = as.factor(c(rep(0, length(cont_vals)), rep(1, length(case_vals))))
        #perm_wcox = wilcox_test(vals ~ treats)
        #perms[z] = pvalue(perm_wcox)
        z=z+1
        if(z%%100000==0){print(z)}
    }
        wcoxs = ifelse(is.nan(wcoxs), 1, wcoxs)
        main_diff = cbind(main_mydiff[,c(1,2,3,4,7)], meansdiff, mediansdiff, conts, cases, means_cont, means_case, medians_cont, medians_case, ps_and_es_matrix, wcoxs, age_meansdiff, age_mediansdiff, hpv_kstest)
        return(main_diff)
        }

get_the_rank <- function(tablez, col_number) {
#Function that returns the rank of a matrix based on one column values
tablez <- tablez[order(tablez[,col_number]),]
rank <- 1:nrow(tablez)
return(rank)
}
regions_from_cpgs = function(maindiff_df, dmr_table, filtered_myobj, size, sample_lim, pval_col, eff_col, pval, eff_size){
    #diff_df = topdiff
    dmrs = dmr_table
    dmrs = dmrs[,c(1,2,3)]
    dmrs$strand='*'
    dmrs$cpg_ig = dmrs$id = paste0(dmrs$V1, dmrs$V2, dmrs$V3)
    dmrs$type = 'dmr'
    colnames(dmrs) = c('chr', 'start', 'end', 'strand', 'id', 'cpg_id', 'type')
    diff_df = subset(maindiff_df, cases>=sample_lim & conts>=sample_lim & maindiff_df[,pval_col]<pval & abs(maindiff_df[,eff_col])>=eff_size)
    print(nrow(diff_df))
    diff_win_df = data.frame(chr=diff_df$chr, start=diff_df$start-(size/2), end=diff_df$start+(size/2), strand=diff_df$strand, cpg_id=paste0(diff_df$chr, diff_df$start, diff_df$end))
    diff_win_df$id=paste0(diff_win_df$chr, diff_win_df$start, diff_win_df$end)
    diff_win_df$type = 'cpg'
    diff_win_df = rbind(diff_win_df, dmrs)
    diff_win_gr <<- makeGRangesFromDataFrame(diff_win_df, keep.extra.columns = TRUE)
    red_diff_win_gr <- reduce(diff_win_gr, ignore.strand=TRUE)
    red_diff_win_gr_df = as.data.frame(red_diff_win_gr)
    red_diff_win_gr_df$id = paste0(red_diff_win_gr_df$seqnames, red_diff_win_gr_df$start, red_diff_win_gr_df$end)
    filtered.myobj_win <- regionCounts(filtered_myobj, red_diff_win_gr, strand.aware=FALSE)
    meth_win <<- unite(filtered.myobj_win,min.per.group=20L)
    rownames(meth_win)=paste0(meth_win$chr, meth_win$start, meth_win$end)
    print(length(paste0(meth_win$chr, meth_win$start, meth_win$end))==length(unique(paste0(meth_win$chr, meth_win$start, meth_win$end))))
    meth_win_per <<- percMethylation(meth_win)
    treatment = as.numeric(pheno_matrix[colnames(meth_win_per),]=='Case')
    win_conts = win_cases = win_means_cont = win_means_case = win_cont_ids = win_case_ids = win_medians_cont = win_medians_case = win_meansdiff = win_mediansdiff = win_wcoxs = win_perms = c()
    z=1
    for(i in red_diff_win_gr_df$id){
        win_cont_vals = as.numeric(na.omit(meth_win_per[i, c(which(treatment==0))]))
        win_case_vals = as.numeric(na.omit(meth_win_per[i, c(which(treatment==1))]))
        win_cont_ids[z] = paste(names(na.omit(meth_win_per[i, c(which(treatment==0))])),collapse=',')
        win_case_ids[z] = paste(names(na.omit(meth_win_per[i, c(which(treatment==1))])),collapse=',')
        win_current_cont = sum(!is.na(meth_win_per[i, c(which(treatment==0))]))
        win_current_case = sum(!is.na(meth_win_per[i, c(which(treatment==1))]))
        win_conts[z]=win_current_cont
        win_cases[z]=win_current_case
        win_means_cont[z]= mean(win_cont_vals)
        win_means_case[z] = mean(win_case_vals)
        win_medians_cont[z] = median(win_cont_vals)
        win_medians_case[z] = median(win_case_vals)
        win_meansdiff[z] = win_means_case[z] - win_means_cont[z]
        win_mediansdiff[z] = win_medians_case[z] - win_medians_cont[z]
        z = z+1
}
    fin_matrix = cbind(red_diff_win_gr_df[,c(1,2,3,4,6)], win_meansdiff, win_mediansdiff, win_conts, win_cases, win_means_cont, win_means_case, win_cont_ids, win_case_ids)
    #windows_diff_df = merge(diff_df, fin_matrix, by.x='cpg_id', by.y='cpg_id')
    return(fin_matrix)
}

check_pred_accuracy <- function(df, id, cluster_ids, phenotype, cv=10, dsplit=0.5){
   #print(id)
    cluster_meth = meth_win_per[df[cluster_ids,'id'],]
    cand_meth = meth_win_per[df[id,'id'],]
    #names(cand_meth) = df[id,'id']
    cand_dataset = rbind(cluster_meth, cand_meth)
    rownames(cand_dataset) = df[c(cluster_ids,id),'id']
    #cand_dataset = cluster_meth
    #print('OK')
    cand_dataset = cand_dataset[ , colSums(is.na(cand_dataset)) == 0]
    cand_dataset = data.frame(t(cand_dataset))
    cand_dataset$class = pheno_matrix[rownames(cand_dataset),]
    #print('OK')
    #control <- trainControl(method="repeatedcv", number=1000, repeats=12, search="grid")
    x = cand_dataset[,c(1:ncol(cand_dataset)-1)]
    y = cand_dataset[,ncol(cand_dataset)]
    #tunegrid <- expand.grid(.mtry=c(1:30))
    accuracy_vec = c()
    res_vec = list()
    for(i in c(1:cv)){
    intrain = createDataPartition(y, p=dsplit, list=FALSE)
    training = cand_dataset[intrain,]
    test = cand_dataset[-intrain,]
    model_1 = randomForest(x=training[,-ncol(training)],y=training[,ncol(training)], ntree=1000)
    model_2 = randomForest(x=test[,-ncol(test)],y=test[,ncol(test)], ntree=1000)
    results_1 = confusionMatrix(predict(model_1, test[,-ncol(test)]), test[,ncol(test)])
    results_2 = confusionMatrix(predict(model_2, training[,-ncol(training)]), training[,ncol(training)])
    accuracy_vec[i] = mean(results_1$overall[1], results_2$overall[1])
    }
    #res_vec[[i]]=results
    #}
    accuracy = mean(accuracy_vec)
    #print(accuracy)
    return(c(as.character(id),accuracy))
}

get_pred_metrics <- function(cm){
    res_vec = c(cm$overall[1], cm$byClass[c(1,2,3)])
    names(res_vec) = c('Accuracy', 'Sensitivity', 'Specificity', 'PPV')
    return(res_vec)
}

data_split_cv <- function(df, cluster_ids, pheno_matrix, cv, part, add_to_list){
    cand_meth = meth_win_per[df[cluster_ids,'id'],]
    cand_dataset = cand_meth
    cand_dataset = cand_dataset[ , colSums(is.na(cand_dataset)) == 0]
    cand_dataset = data.frame(t(cand_dataset))
    cand_dataset$class = pheno_matrix[rownames(cand_dataset),]
    #data_age = age_matrix[rownames(cand_dataset),]
    x = cand_dataset[,c(1:ncol(cand_dataset)-1)]
    y = cand_dataset[,ncol(cand_dataset)]
    res_noage = res_age = list()
    cols = ncol(cand_dataset)
    for(i in c(1:cv)){
        intrain = createDataPartition(y, p=part, list=FALSE)
        training = cand_dataset[intrain,]
        test = cand_dataset[-intrain,]
        rf_noage <- randomForest(x=training[,-ncol(training)],y=training[,ncol(training)], ntree=1000)
        #rf_age <- train(x=training[,-(cols-1)],y=training[,'class'],method='rf',tuneGrid=tunegrid, trControl=train_control)
        res_noage[[i]] <- get_pred_metrics(confusionMatrix(predict(rf_noage, test[,-cols]),test[,'class']))
        #res_age[[i]] <- get_pred_metrics(confusionMatrix(predict(rf_age, test[,-c(cols-1)]),test[,'class']))
    }
    noage_fres = c(unlist(ci.mean(unlist(lapply(res_noage, function(x) x[1])))[c(1,3,4)]), unlist(ci.mean(unlist(lapply(res_noage, function(x) x[2])))[c(1,3,4)]),
                  unlist(ci.mean(unlist(lapply(res_noage, function(x) x[3])))[c(1,3,4)]),
                  unlist(ci.mean(unlist(lapply(res_noage, function(x) x[4])))[c(1,3,4)]))
    return(c(add_to_list,noage_fres))
}

feature_selection = function(diff_regs, starting_point, shared_samples_limit, cand_shared_range){
    outlist = c()
    print(paste0('Starting with starting point: ', starting_point))
    print('Processing Starting point')
    #starting_point=c(5255:5259)
    #shared_samples_limit = 170
    iter=0
    accuracy=0
    seed=324432
    set.seed(seed)
    cluster_idx = starting_point
    print(paste0('Samples on the current cluster: ',cluster_idx))
    samples_no = length(cluster_idx)
    cluster_conts = Reduce(intersect, lapply(as.character(diff_regs[cluster_idx,'win_cont_ids']), function(x) unlist(strsplit(x, ','))))
    cluster_cases = Reduce(intersect, lapply(as.character(diff_regs[cluster_idx,'win_case_ids']), function(x) unlist(strsplit(x, ','))))
    shared_samples_no = length(cluster_conts)+length(cluster_cases)
    outlist <- c(iter, paste(cluster_idx, collapse=','), shared_samples_no, accuracy)
    rems_idx = which(diff_regs$win_conts+diff_regs$win_cases>shared_samples_limit & !rownames(diff_regs) %in% cluster_idx)
    #print(rems_idx)
    print('Starting point Processed')
    while(shared_samples_no > shared_samples_limit){
        print('')
        print('----------------------------------------------------------')
        iter=iter+1
        print(paste0('Processing candidates on iteration no ',iter))
        cand_shared_vec <- unlist(lapply(lapply(lapply(as.character(diff_regs[rems_idx,'win_cont_ids']),function(x) unlist(strsplit(x, ','))),
                                         function(x) list(cluster_conts, x)), function(x) length(Reduce(intersect, x))))+unlist(lapply(lapply(lapply(as.character(diff_regs[rems_idx,'win_case_ids']),
                                         function(x) unlist(strsplit(x, ','))), function(x) list(cluster_cases, x)), function(x) length(Reduce(intersect, x))))
        if(cand_shared_range!=0){
        cand_idx_vec <- rownames(diff_regs[rems_idx,][which(cand_shared_vec %in% c((max(cand_shared_vec)-cand_shared_range):max(cand_shared_vec))),])
        }
        if(cand_shared_range==0){
            cand_idx_vec <- rems_idx
        }
        if(max(cand_shared_vec)<shared_samples_limit){break}
        #cand_idx_vec <- rems_idx
        #print(cand_idx_vec)
        print(paste0('Got ', length(cand_idx_vec), ' candidates!'))
        print('Running the machine learning algorithm on them!')
        oper = foreach(i=cand_idx_vec) %dopar% check_pred_accuracy(diff_regs, i, cluster_idx, pheno_matrix)
        best_idx = oper[[which.max(unlist(lapply(oper, function(x) x[2])))]][1]
        accuracy = oper[[which.max(unlist(lapply(oper, function(x) x[2])))]][2]
        print('Optimal feature found!')
        print(paste0('Accuracy: ', accuracy))
        #print(best_idx)
        cluster_idx = c(cluster_idx, best_idx)
        #print(cluster_idx)
        samples_no = length(cluster_idx)
        cluster_conts = Reduce(intersect, lapply(as.character(diff_regs[cluster_idx,'win_cont_ids']), function(x) unlist(strsplit(x, ','))))
        cluster_cases = Reduce(intersect, lapply(as.character(diff_regs[cluster_idx,'win_case_ids']), function(x) unlist(strsplit(x, ','))))
        #cluster_meth = meth_win_per[diff_regs[cluster_idx,'id'],]
        shared_samples_no = length(cluster_conts)+length(cluster_cases)
        print(paste0('Number of shared samples: ', shared_samples_no))
        rems_idx = which(diff_regs$win_conts+diff_regs$win_cases>shared_samples_limit & !rownames(diff_regs) %in% cluster_idx)
        #print(rems_idx)
        curlist <- c(iter, paste(cluster_idx, collapse=','), shared_samples_no, accuracy)
        #print(curlist)
        outlist <- c(outlist, iter, paste(cluster_idx, collapse=','), shared_samples_no, accuracy)
        print(paste0('Samples in cluster: ', samples_no))
        print(paste0('End of iter ', iter))
        print('----------------------------------------------------------')
        print('')
    }
    return(outlist)
}
#fin_report = feature_selection(diff_regs, starting_point=st, shared_samples_limit=170)

fin_list = list()
  for(i in c(1:nrow(all_combs))){
      print(paste0('Starting point: ',all_combs[i,1]))
      print(paste0('Minimum shared: ', all_combs[i,2]))
      print(paste0('Close range: ', all_combs[i,3]))
      fin_list[[i]]=c(feature_selection(diff_regs, starting_point=as.character(all_combs[i,1]), shared_samples_limit=as.numeric(all_combs[i,2]), cand_shared_range=as.numeric(all_combs[i,3])))
  }

evaluate_clusters <- function(final_list, meth_win_per, pheno_matrix, age_matrix){
    fin = data.frame(matrix(unlist(lapply(final_list, function(x) x[-c(1,2)])), ncol=4, byrow=TRUE), stringsAsFactors=FALSE)
    colnames(fin)=c("iter", "features", "samples", "accu")
    fin$no_of_ids = unlist(lapply(fin$features, function(x) length(unlist(strsplit(x, ',')))))
    select_fin = subset(fin, as.numeric(no_of_ids)>=3 & as.numeric(accu)>0.7 & samples>=160)
    oper = foreach(i=1:nrow(select_fin)) %dopar% data_split_cv(diff_regs, cluster_ids = unlist(strsplit(select_fin[i,2],',')), pheno_matrix, cv=1000, part=0.75, add_to_list=c(select_fin[i,2], select_fin[i,3], select_fin[i,5]))
    fin_df = data.frame(matrix(unlist(oper), ncol=15, byrow=TRUE), stringsAsFactors=FALSE)
    colnames(fin_df) = c("ids", "sample_size", "no_of_features", "Accuracy", "Acc_low", "Acc_high", "Sensitivity", "Sens_low", "Sens_high", "Specificity", "Spec low", "Spec high", "PPV", "PPV_low", "PPV_high")
    return(fin_df)
}
#best_clusters = evaluate_clusters(fin_list, meth_win_per, pheno_matrix, age_matrix)
# }
# mat = c()
# for(i in length(fin_list)){
#     mat = c(mat, fin_list[[i]])
#
# }


best_comb = function(windows_diff, minimum_samples){
    #windows_diff=fin_df
    windows_diff=regs
    best=subset(windows_diff, win_conts>=85 & win_cases>=85)
    rows = c(1:nrow(best))
    all_combs = do.call("c", lapply(seq_along(rows), function(i) combn(rows, i, FUN = list)))
    ncombs = combs = conts = cases = c()
    for(z in c(1:length(all_combs))){
        combs[z] = paste(unlist(all_combs[z]), collapse=',')
        ncombs[z] = length(combs[z])
        conts[z] = length(Reduce(intersect, lapply(lapply(unlist(lapply(best[unlist(all_combs[z]),'win_cont_ids'], function (x) as.character(x))), function(x) strsplit(x, ',')),function(x) unlist(x))))
        cases[z] = length(Reduce(intersect, lapply(lapply(unlist(lapply(best[unlist(all_combs[z]),'win_case_ids'], function (x) as.character(x))), function(x) strsplit(x, ',')),function(x) unlist(x))))
    }
    return(cbind(combs, ncombs, conts, cases))
}

get_the_OR <- function(gores, sigs, all){
	#Function that takes the result of the go enrichment analysis as its argument
	#and returns the Fisher's Exact Odds Ratio for each reported GO term.
	#Arguments:
		#gores: Output of the 'GenTable' function of the 'topgo' package.
		#sigs: Vector of gene names which are DE and annotated to GO terms.
		#all: Vector of the all the gene names included in the analysis.

	s1 <- gores$Significant
	s2 <- sigs - s1
	b1 <- gores$Annotated
	b2 <- all - b1
	OR <- (s1*b2)/(s2*b1)

return(OR)
}

age_b <- function(z){
    ageconts = age_matrix[Reduce(intersect,lapply(as.character(diff_regs[unlist(strsplit(selection[z,1], ',')), 12]), function(x) unlist(strsplit(x, ',')))),]
    agecases = age_matrix[Reduce(intersect,lapply(as.character(diff_regs[unlist(strsplit(selection[z,1], ',')), 13]), function(x) unlist(strsplit(x, ',')))),]
    boxplot(ageconts, agecases, ylim=c(1,max(c(ageconts,agecases))+10), names=c(paste0('Controls N=',length(ageconts),' mean=',format(round(mean(ageconts),2),nsmall=2)),paste0('Cases N=',length(agecases),' mean=', format(round(mean(agecases),2),nsmall=2))), col=c("#3cba54","orange"))
}

go_analysis <- function(allgenes, path, name){
	#Function that takes the output of the 'decide_important' or 'merge_gene_sets' function
	#as its argument and returns a list of 9 objects, 3 for each specific ontology (BP=Biological Process,
	#CC=Cellular Component, MF=Molecular Function). Four different combinations of the topgo data object
	#were used for each ontology: 'elim' algorithm with KS test, 'elim' algorithm with exact fisher's test,
	#'weight' algorithm with exact fisher's test and 'weight01algorithm' with exact fisher's test. The result
	#matrices are being saved as CSV files, one for each ontology.
	#First three elements: Three topgodata objects built for the analysis (see topgo vignette). BP, CC, MF respectively.
	#Second three elements: The result matrices of the 'GenTable' function (see topgo vignette.) BP, CC, MF respectively.
						   #Each result matrix contains the results from all different combinations of algorithms and
						   #statistical tests. These matrices also include columns reporting the names of significant genes of
						   #each GO term (a column for the up-regulated and a column for the down-regulated genes) and a column
						   #reporting the Fisher's test Odds Ratio calculated for each term.
	#Third three elemets: Lists of DE genes assigned to each GO term. BP, CC, MF respectively. BP, CC, MF respectively.
	#Arguments:
		#allgenes: output of the 'decide_important' or 'merge_gene_sets' function.
		#path: path to store the output CSV files.

	#dir.create(path=path, showWarnings = TRUE, recursive = FALSE, mode = "0777")

	type <- deparse(substitute(type))
	nameit <- name

	#Annotate genes to Biological Process (BP) GO terms
	topgo_bp  <- new("topGOdata", ontology = "BP", allGenes = allgenes,
	geneSel = function(k) k > 0, annot = annFUN.org, mapping="org.Hs.eg.db", ID="symbol")

	#Genes annotated to Biological Process (BP) GO terms
	sigs_bp <- as.numeric(table(topgo_bp@feasible[which(topgo_bp@allScores == 1)])["TRUE"])
	all_bp <- length(topgo_bp@feasible[topgo_bp@feasible==TRUE])

	#Annotate genes to Cellular Component (CC) GO terms
	topgo_cc  <- new("topGOdata", ontology = "CC", allGenes = allgenes,
	geneSel = function(k) k > 0, annot = annFUN.org, mapping="org.Hs.eg.db", ID="symbol")

	#Genes annotated to Cellular Component (CC) GO terms
	sigs_cc <- as.numeric(table(topgo_cc@feasible[which(topgo_bp@allScores == 1)])["TRUE"])
	all_cc <- length(topgo_cc@feasible[topgo_bp@feasible==TRUE])

	#Annotate genes to Molecular Function (MF) GO terms
	topgo_mf  <- new("topGOdata", ontology = "MF", allGenes = allgenes,
	geneSel = function(k) k > 0, annot = annFUN.org, mapping="org.Hs.eg.db", ID="symbol")

	#Genes annotated to Molecular Function (MF) GO terms
	sigs_mf <- as.numeric(table(topgo_mf@feasible[which(topgo_bp@allScores == 1)])["TRUE"])
	all_mf <- length(topgo_mf@feasible[topgo_bp@feasible==TRUE])

	#Gemes reported up- or down- regulated
	#ups <- rownames(subset(allgenes[[3]], log2FoldChange > 0))
	#downs <- rownames(subset(allgenes[[3]], log2FoldChange < 0))

	#Run the tests for BP, CC and MF
	topgo_t1_bp <- runTest(topgo_bp, algorithm = "elim", statistic = "ks")
	topgo_t1_cc <- runTest(topgo_cc, algorithm = "elim", statistic = "ks")
	topgo_t1_mf <- runTest(topgo_mf, algorithm = "elim", statistic = "ks")
	topgo_t2_bp <- runTest(topgo_bp, algorithm = "elim", statistic = "fisher")
	topgo_t2_cc <- runTest(topgo_cc, algorithm = "elim", statistic = "fisher")
	topgo_t2_mf <- runTest(topgo_mf, algorithm = "elim", statistic = "fisher")
	topgo_t3_bp <- runTest(topgo_bp, algorithm = "weight", statistic = "fisher")
	topgo_t3_cc <- runTest(topgo_cc, algorithm = "weight", statistic = "fisher")
	topgo_t3_mf <- runTest(topgo_mf, algorithm = "weight", statistic = "fisher")
	topgo_t4_bp <- runTest(topgo_bp, algorithm = "weight01", statistic = "fisher")
	topgo_t4_cc <- runTest(topgo_cc, algorithm = "weight01", statistic = "fisher")
	topgo_t4_mf <- runTest(topgo_mf, algorithm = "weight01", statistic = "fisher")

	#Get the genes annotated per term
	all_GO_bp <- genesInTerm(topgo_bp)
	all_GO_cc <- genesInTerm(topgo_cc)
	all_GO_mf <- genesInTerm(topgo_mf)

	#Get the significant genes annotated per term
	SAM_bp <- lapply(all_GO_bp,function(x) x[x %in% names(allgenes[allgenes==1])] )
	SAM_cc <- lapply(all_GO_cc,function(x) x[x %in% names(allgenes[allgenes==1])] )
	SAM_mf <- lapply(all_GO_mf,function(x) x[x %in% names(allgenes[allgenes==1])] )

	#Summarise the BP results
	topgo_res_bp <- GenTable(topgo_bp, elimKS = topgo_t1_bp, elimFisher = topgo_t2_bp,
		weightFisher = topgo_t3_bp, weight01Fisher = topgo_t4_bp, orderBy = "weight01Fisher",
		ranksOf = "weight01Fisher", topNodes = length(topgo_bp@graph@nodes))
	topgo_res_bp$OR <- get_the_OR(topgo_res_bp, sigs_bp, all_bp) #add the OR as well

	#Find which genes are up- or down- regulated for each reported BP term
	# upgenes <- c()
	# downgenes <- c()
    #
	# for(i in topgo_res_bp[,1]){
    #
	# 	upgenes[i]	<- paste(c(SAM_bp[[i]][SAM_bp[[i]] %in% ups]), collapse=", ")
	# 	downgenes[i] <- paste(c(SAM_bp[[i]][SAM_bp[[i]] %in% downs]), collapse=", ")
	# }
    #
	# topgo_res_bp$upgenes <- upgenes
	# topgo_res_bp$downgenes <- downgenes

	#Summarise the CC results
	topgo_res_cc <- GenTable(topgo_cc, elimKS = topgo_t1_cc, elimFisher = topgo_t2_cc,
		weightFisher = topgo_t3_cc, weight01Fisher = topgo_t4_cc, orderBy = "weight01Fisher",
		ranksOf = "weight01Fisher", topNodes = length(topgo_cc@graph@nodes))
	topgo_res_cc$OR <- get_the_OR(topgo_res_cc, sigs_cc, all_cc) #add the OR as well

	# #Find which genes are up- or down- regulated for each reported CC term
	# upgenes <- c()
	# downgenes <- c()
    #
	# for(i in topgo_res_cc[,1]){
	# 	upgenes[i] <- paste(c(SAM_cc[[i]][SAM_cc[[i]] %in% ups]), collapse=", ")
	# 	downgenes[i] <- paste(c(SAM_cc[[i]][SAM_cc[[i]] %in% downs]), collapse=", ")
    #
	# }
    #
	# topgo_res_cc$upgenes <- upgenes
	# topgo_res_cc$downgenes <- downgenes

	#Summarise the MF results
	topgo_res_mf <- GenTable(topgo_mf, elimKS = topgo_t1_mf, elimFisher = topgo_t2_mf,
		weightFisher = topgo_t3_mf, weight01Fisher = topgo_t4_mf, orderBy = "weight01Fisher",
		ranksOf = "weight01Fisher", topNodes = length(topgo_mf@graph@nodes))
	topgo_res_mf$OR <- get_the_OR(topgo_res_mf, sigs_mf, all_mf) #add the OR as well

	# #Find which genes are up- or down- regulated for each reported MF term
	# upgenes <- c()
	# downgenes <- c()
	# z <- 1
	# for(i in topgo_res_mf[,1]){
	# 	upgenes[i] <- paste(c(SAM_mf[[i]][SAM_mf[[i]] %in% ups]), collapse=", ")
	# 	downgenes[i] <- paste(c(SAM_mf[[i]][SAM_mf[[i]] %in% downs]), collapse=", ")
	# 	z=z+1
	# }
    #
	# topgo_res_mf$upgenes <- upgenes
	# topgo_res_mf$downgenes <- downgenes

	#Write the output CSV files
	# write.csv(data.frame(lapply(topgo_res_bp, as.character), stringsAsFactors=FALSE), file=paste(path,nameit,"_bp.csv",sep=""))
	# write.csv(data.frame(lapply(topgo_res_mf, as.character), stringsAsFactors=FALSE), file=paste(path,nameit,"_mf.csv",sep=""))
	# write.csv(data.frame(lapply(topgo_res_cc, as.character), stringsAsFactors=FALSE), file=paste(path,nameit,"_cc.csv",sep=""))

	return(list(topgo_bp, topgo_cc, topgo_mf, topgo_res_bp, topgo_res_cc, topgo_res_mf, SAM_bp, SAM_cc, SAM_mf))
	}

go_barplot <- function(topgo, path, name, test, mars, wi, he, cen, cela, cele, legpos, save=TRUE){
	#Function that takes the output of the 'go_analysis' function as its argument and returns a plot
	#visualising the number of each of the significant genes of the top GO terms for all three ontologies
	#(BP, CC, MF). top GO terms are the GO terms fullifilling these criteria: test pvalue < 0.05,
	#OR>10, Significant genes >= 2).
	#Arguments:
		#go_res: Output of the 'go_analysis' function
		#save: Boolean. Whether to save or not the output matrix as a CSV file.
		#path: Path to save the CSV file.
		#test: For which test should the pvalue of a term be less than 0.05 in order
			  #for it to be included in the plot. Could be 'Weight01Fisher', 'Weight' or 'elim'.
		#mars: margins of the plot. (to be used in par(mar=c())).
		#wi: Weight parameter of the PNG
		#he: Height parameter of the PNG
		#cen: Font size of yaxis names and inside-the-bars text. To be passed to cex graphical parameter
		#cela: Font size of xaxis label. To be passed to cex graphical parameter
		#cele: Font size of the legend. To be passed to cex graphical parameter
		#legpos:Keyword to be used to position the legend. Accepted keywords:
		       #"bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right"
		       #and "center"

	#Name and prefix to be used for naming the output PNG

	#Top terms of the BP ontology analysis
	bp <- subset(topgo[[4]], as.numeric(Significant) >=2 & as.numeric(topgo[[4]][,test]) < 0.05 & as.numeric(Significant)>=4)
	bp <- bp[order(-as.numeric(bp[,test])),]

	#Top terms of the CC ontology analysis
	cc <- subset(topgo[[5]], as.numeric(Significant) >=2 & as.numeric(topgo[[5]][,test]) < 0.05 & as.numeric(Significant)>=3) #& as.numeric(OR) >= 10)
	cc <- cc[order(-as.numeric(cc[,test])),]

	#Top terms of the MF ontology analysis
	mf <- subset(topgo[[6]], as.numeric(Significant) >=2 & as.numeric(topgo[[6]][,test]) < 0.05 & as.numeric(Significant)>=3)
	mf <- mf[order(-as.numeric(mf[,test])),]

	bps <- nrow(bp)
	ccs <- nrow(cc)
	mfs <- nrow(mf)


	goterms <- Term(GOTERM)

	#Extract information from 'go_analysis' output object.
	#Create the empty lists for the loops
	gos_bp <- c()
	gos_val_bp <- c()
	gos_ann_bp <- c()
	gos_or_bp <- c()
	gos_p_bp <- c()

	gos_mf <- c()
	gos_val_mf <- c()
	gos_ann_mf <- c()
	gos_or_mf <- c()
	gos_p_mf <- c()

	gos_cc <- c()
	gos_val_cc <- c()
	gos_ann_cc <- c()
	gos_or_cc <- c()
	gos_p_cc <- c()

	#BP information
	for(i in 1:bps){
	gos_bp[i] <- goterms[[bp[i,1]]]
	gos_val_bp[i] <- as.numeric(bp[i,4])
	gos_ann_bp[i] <- as.numeric(bp[i,3])
	gos_or_bp[i] <- round(as.numeric(bp[i,10]), digits=0)
	gos_p_bp[i] <- format(as.numeric(bp[i,test]), scientific=TRUE)
	}

	#MF information
	for(i in 1:mfs){
	gos_mf[i] <- goterms[[mf[i,1]]]
	gos_val_mf[i] <- as.numeric(mf[i,4])
	gos_ann_mf[i] <- as.numeric(mf[i,3])
	gos_or_mf[i] <- round(as.numeric(mf[i,10]), digits=0)
	gos_p_mf[i] <- format(as.numeric(mf[i,test]), scientific=TRUE)
	}

	#CC information
  if(nrow(cc)>0){
  	for(i in 1:ccs){
     	gos_cc[i] <- goterms[[cc[i,1]]]
    	gos_val_cc[i] <- as.numeric(cc[i,4])
    	gos_ann_cc[i] <- as.numeric(cc[i,3])
    	gos_or_cc[i] <- round(as.numeric(cc[i,10]), digits=0)
    	gos_p_cc[i] <- format(as.numeric(cc[i,test]), scientific=TRUE)
    }
  } else {
    gos_cc <- c()
  	gos_val_cc <- c()
  	gos_ann_cc <- c()
  	gos_or_cc <- c()
  	gos_p_cc <- c()
	}


	#Combine all the information
	gos <- c(gos_cc, gos_mf, gos_bp)
	gos_val <- c(gos_val_cc, gos_val_mf, gos_val_bp)
	gos_ann <- c(gos_ann_cc, gos_ann_mf, gos_ann_bp)
	gos_or <- c(gos_or_cc, gos_or_mf, gos_or_bp)
	gos_p <- c(gos_p_cc, gos_p_mf, gos_p_bp)

  max_x <- max(gos_val)

  if(missing(wi)){
    if(max_x %in% c(1:5)){
      wi <- max_x+24
    } else if(max_x %in% c(6:9)){
      wi <- max_x+11
    } else if(max_x %in% c(10:14)){
      wi <- max_x+12
    } else {
      wi <- (max_x*25)/11
    }
	}

	out_of <- paste(gos_val, gos_ann, sep="/")


	#Create the plot for less than 25 BP topgo terms

	if(bps < 25){
		if(save==TRUE){
			png(filename=paste(path, name, ".png", sep=""), units="in", width=wi, height=he, res=600)
		}

		par(mar=mars)
		y <- barplot(gos_val, horiz=TRUE, names.arg=gos, col=c(rep("#4885ed", ccs),rep("#3cba54", mfs), rep("orange", bps)),
		cex.names=cen, xaxt="n", las=2, xlab="Number of differentially expressed genes", cex.lab=cela)
		x <- 0.5*gos_val
		stats <- paste(out_of, gos_or, gos_p, sep=", ")
		text(x, y, stats, cex=cen)
		axis(1, at=seq(0, (max_x+1), by=1), labels=c(0:(max_x+1)), xlim=c(0,(max_x+1)), cex.axis=1.7)
		#legend(legpos, pch=15, col=c("orange", "#3cba54", "#4885ed"), legend=c("Biological Process", "Molecular Function", "Cellular Component"), cex=cele, bty="n")
		abline(v=0)

		if(save==TRUE){
			dev.off()
		}
	}


	if(bps > 25){
    if(save==TRUE){
			png(filename=paste(path, name, ".png", sep=""), units="in", width=wi, height=he, res=600)
		}

    par(mar=mars)
		y <- barplot(log2(gos_val), horiz=TRUE, names.arg=gos, col=c(rep("#4885ed", ccs),rep("#3cba54", mfs), rep("orange", bps)),
		cex.names=cen, xaxt="n", las=2, xlab="Number of differentially expressed genes", cex.lab=cela)
		x <- 0.5*log2(gos_val)
		stats <- paste(out_of, gos_or, gos_p, sep=", ")
		text(x, y, stats, cex=cen)
		axis(1)
		#legend(legpos, pch=15, col=c("orange", "#3cba54", "#4885ed"), legend=c("Biological Process", "Molecular Function", "Cellular Component"), cex=cele, bty="n")
		abline(v=0)

		if(save==TRUE){
			dev.off()
		}
	}

	} #End of function

check_raw_data <- function(lab_nos=as.character(globalmaster$lab_no), paths=globalmaster$path, chr, start, end){
    raw_meths=c()
    for(i in c(1:length(paths))){
        p = paths[i]
        #print(p)
        raw_vals = fread(p, sep='\t', header=FALSE)
        raw_spec = subset(raw_vals, V1==chr & V2==start & V2==end)
        #rm(p)
        if(nrow(raw_spec)==0){
            raw_meth=NA
        } else {
            if(as.numeric(raw_spec$V4)+as.numeric(raw_spec$V5)<=9){
                raw_meth=NA
            } else {
            raw_meth = (raw_spec$V4/(raw_spec$V5+raw_spec$V4))*100
            }
        }
        raw_meths[i]=raw_meth
    }
    names(raw_meths) = lab_nos
    return(rbind(raw_meths, meth.min_5_per[paste0(chr,start),names(raw_meths)]))
}

create_bed_for_combp <- function(df, col_for_p, col_for_eff, eff_limit, path_name){
    bedp <- df[,c(1,2,3,col_for_p, col_for_eff)]
    colnames(bedp) <- c('chrom', 'start', 'end', 'pvalue', 't')
    bedp$end = bedp$start+1
    newp = ifelse(bedp$t<eff_limit, 1, bedp$pvalue)
    bedp$pvalue <- newp
    write.table(bedp, path_name, sep='\t', row.names=FALSE, quote=FALSE)
    return(bedp)
}

draw_cpgs = function(mdiff, as.score='median', gr, cols=c("chr","start","end","strand","mediansdiff","meansdiff","wcoxs")){
    epb = subset(mdiff, start %in% c(as.numeric(gr[2]):as.numeric(gr[3])) & chr==gr[1])
    epb = epb[,cols]
    epb_gr = makeGRangesFromDataFrame(data.frame(epb, stringsAsFactors=TRUE), keep.extra.columns=TRUE)
    epb_gr$mediansdiff = as.numeric(as.character(epb_gr$mediansdiff))
    epb_gr$meansdiff = as.numeric(as.character(epb_gr$meansdiff))
    epb_gr$color = 'gray'
    if(as.score=='median'){
    try(epb_gr[as.numeric(as.character(epb_gr$wcoxs))>0.05 & epb_gr$mediansdiff>0,]$color <- '#e6ffe6', silent=TRUE)
    try(epb_gr[as.numeric(as.character(epb_gr$wcoxs))<0.05 & epb_gr$mediansdiff>0,]$color <- '#269900', silent=TRUE)
    try(epb_gr[as.numeric(as.character(epb_gr$wcoxs))>0.05 & epb_gr$mediansdiff<0,]$color <- '#ffe6e6', silent=TRUE)
    try(epb_gr[as.numeric(as.character(epb_gr$wcoxs))<0.05 & epb_gr$mediansdiff<0,]$color <- '#ff0000', silent=TRUE)
    epb_gr$score = abs(epb_gr$mediansdiff)
    }
    if(as.score=='mean'){
    try(epb_gr[as.numeric(as.character(epb_gr$wcoxs))>0.05 & epb_gr$meansdiff>0,]$color <- '#e6ffe6', silent=TRUE)
    try(epb_gr[as.numeric(as.character(epb_gr$wcoxs))<0.05 & epb_gr$meansdiff>0,]$color <- '#269900', silent=TRUE)
    try(epb_gr[as.numeric(as.character(epb_gr$wcoxs))>0.05 & epb_gr$meansdiff<0,]$color <- '#ffe6e6', silent=TRUE)
    try(epb_gr[as.numeric(as.character(epb_gr$wcoxs))<0.05 & epb_gr$meansdiff<0,]$color <- '#ff0000', silent=TRUE)
    epb_gr$score = abs(epb_gr$meansdiff)
    }
    trs = geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg19.knownGene, org.Hs.eg.db, chrom=gr[1], start=as.numeric(gr[2]), end=as.numeric(gr[3]), strand=gr[4])
    summ_to_transcripts = data.frame(matrix(unlist(lapply(trs, function(y) lapply(unique(y$dat$feature),
                          function(x) c(levels(seqnames(y$dat[y$dat$feature==x,]))[1], start(y$dat[y$dat$feature==x,])[1],
                          tail(end(y$dat[y$dat$feature==x,]),n=1), y$dat[y$dat$feature==x,]$feature[1],
                          y$dat[y$dat$feature==x,]$transcript[1])))), ncol=5, byrow=TRUE), stringsAsFactors=FALSE)
    colnames(summ_to_transcripts) = c('seqnames','start','end','feature','transcript')
    summ_to_transcripts2 = data.frame(seqnames=summ_to_transcripts$seqnames, start=apply(summ_to_transcripts[, 2:3], 1, min),
                           end=apply(summ_to_transcripts[, 2:3], 1, max), feature=summ_to_transcripts$feature,
                           transcript=summ_to_transcripts$transcript)
    summ_to_transcripts_gr = makeGRangesFromDataFrame(summ_to_transcripts2, keep.extra.columns=TRUE)
    flevs = c('CDS','utr3','utr5','ncRNA')
    feats = factor(summ_to_transcripts_gr$feature, levels=flevs)
    summ_to_transcripts_gr$featureLayerID = summ_to_transcripts_gr$transcript
    summ_to_transcripts_gr$height = unit(12, "points")
    cols = c("#51C6E6", "#FF8833", "#DFA32D", "darkgray")
    summ_to_transcripts_gr$fill = cols[feats]
    names(summ_to_transcripts_gr) = summ_to_transcripts_gr$feature
    return(list(summ_to_transcripts_gr, epb_gr))
}

case_vs_control_box <- function(methper, rows, parfrow){
    par(mfrow=parfrow, cex.lab=5)
    for(row in rows){
        cases = as.numeric(na.omit(methper[row, c(which(treatment==1))]))
        conts = as.numeric(na.omit(methper[row, c(which(treatment==0))]))
        ncases = length(cases)
        nconts = length(conts)
        labcon = paste0('Controls', '(N=', nconts, ')')
        labcas = paste0('Cases', '(N=', ncases, ')')
        boxplot(conts, cases, col=c("#3cba54", "orange"), names=c(labcon, labcas), cex.names=5)
    }
}

