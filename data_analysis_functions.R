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
    conts = cases = means_cont = means_case = medians_cont = medians_case = meansdiff = mediansdiff = wcoxs = perms = c()
    z=1
    main_mydiff = getData(mydiffs[[1]])
    #length(conts) = length(cases) = length(means_cont) = length(means_case) = length(medians_cont) = length(medians_case) = length(meansdiff) = length(mediansdiff) = length(wcoxs) = length(perms) = nrow(main_mydiff)
    for(i in as.numeric(rownames(main_mydiff))){
        cont_vals = as.numeric(na.omit(meth.min_5_per[i, c(which(treatment==0))]))
        case_vals = as.numeric(na.omit(meth.min_5_per[i, c(which(treatment==1))]))
        current_cont = sum(!is.na(meth.min_5_per[i, c(which(treatment==0))]))
        current_case = sum(!is.na(meth.min_5_per[i, c(which(treatment==1))]))
        conts[z]=current_cont
        cases[z]=current_case
        means_cont[z]= mean(cont_vals)
        means_case[z] = mean(case_vals)
        medians_cont[z] = median(cont_vals)
        medians_case[z] = median(case_vals)
        meansdiff[z] = means_case[z] - means_cont[z]
        mediansdiff[z] = medians_case[z] - medians_cont[z]
        wcoxs[z] = wilcox.test(cont_vals, case_vals)$p.value
        vals = c(cont_vals, case_vals)
        treats = as.factor(c(rep(0, length(cont_vals)), rep(1, length(case_vals))))
        #perm_wcox = wilcox_test(vals ~ treats)
        #perms[z] = pvalue(perm_wcox)
        z=z+1
        if(z%%100000==0){print(z)}
    }
        wcoxs = ifelse(is.nan(wcoxs), 1, wcoxs)
        main_diff = cbind(main_mydiff[,c(1,2,3,4,7)], meansdiff, mediansdiff, conts, cases, means_cont, means_case, ps_and_es_matrix, wcoxs)
        return(main_diff)
        }

get_the_rank <- function(tablez, col_number) {
#Function that returns the rank of a matrix based on one column values
tablez <- tablez[order(tablez[,col_number]),]
rank <- 1:nrow(tablez)
return(rank)
}

regions_from_cpgs = function(maindiff_df, dmr_table, size, sample_lim, pval, eff_size){
    #diff_df = topdiff
    dmrs = dmr_table
    size=100
    dmrs = dmrs[,c(1,2,3)]
    dmrs$strand='*'
    dmrs$cpg_ig = dmrs$id = paste0(dmrs$V1, dmrs$V2, dmrs$V3)
    colnames(dmrs) = c('chr', 'start', 'end', 'strand', 'id', 'cpg_id')
    diff_df = subset(maindiff_df, cases>=sample_lim & conts>=sample_lim & wcox_pvalue<pval & abs(meansdiff)>=eff_size)
    diff_win_df = data.frame(chr=diff_df$chr, start=diff_df$start-(size/2), end=diff_df$start+(size/2), strand=diff_df$strand, cpg_id=paste0(diff_df$chr, diff_df$start))
    diff_win_df$id=paste0(diff_win_df$chr, diff_win_df$start, diff_win_df$end)
    diff_win_df = rbind(diff_win_df, dmrs)
    diff_win_gr <- makeGRangesFromDataFrame(diff_win_df, keep.extra.columns = TRUE)
    filtered.myobj_win <<- regionCounts(filtered.myobj_5, diff_win_gr, strand.aware=FALSE)
    meth_win <<- unite(filtered.myobj_win,min.per.group=20L)
    rownames(meth_win)=paste0(meth_win$chr, meth_win$start, meth_win$end)
    meth_win_per <<- percMethylation(meth_win)
    win_conts = win_cases = win_means_cont = win_means_case = win_cont_ids = win_case_ids = win_medians_cont = win_medians_case = win_meansdiff = win_mediansdiff = win_wcoxs = win_perms = c()
    z=1
    for(i in diff_win_df$id){
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
    fin_matrix = cbind(diff_win_df, win_meansdiff, win_mediansdiff, win_conts, win_cases, win_means_cont, win_means_case, win_cont_ids, win_case_ids)
    #windows_diff_df = merge(diff_df, fin_matrix, by.x='cpg_id', by.y='cpg_id')
    return(fin_matrix)
}

check_pred_accuracy <- function(df, id, cluster_ids, phenotype, cv=10){
   #print(id)
    cluster_meth = meth_win_per[df[cluster_ids,'id'],]
    cand_meth = meth_win_per[df[id,'id'],]
    #names(cand_meth) = df[id,'id']
    cand_dataset = rbind(cluster_meth, cand_meth)
    rownames(cand_dataset) = diff_regs[c(cluster_ids,id),'id']
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
    for(i in c(1:100)){
        intrain = createDataPartition(y, p=0.75, list=FALSE)
        training = cand_dataset[intrain,]
        test = cand_dataset[-intrain,]
        model = randomForest(x=training[,-ncol(training)],y=training[,ncol(training)], ntree=1000)
        results = confusionMatrix(predict(model, test[,-ncol(test)]), test[,ncol(test)])
        accuracy_vec[i] = results$overall[1]
        res_vec[[i]]=results
    }
    mean(accuracy_vec)
    #print(accuracy)
    return(accuracy)
}

accs=c()
sams=c()
for(cl in cluster_lists){
if(length(cl)>1){
acc = check_pred_accuracy(diff_regs, id=NULL, cl, pheno_matrix, cv=100)
accs = c(accs, acc)
sams = c(sams, paste(cl, collapse=','))
print(cl)
print(acc)
}
}


feature_selection = function(diff_regs, starting_point, shared_samples_limit){
    outlist = c()
    print('Processing Starting point')
    #starting_point=c(5255:5259)
    #shared_samples_limit = 170
    iter=0
    accuracy=0
    cluster_idx = starting_point
    print(cluster_idx)
    samples_no = length(cluster_idx)
    print(samples_no)
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
        cand_idx_vec <- rownames(diff_regs[rems_idx,][which(cand_shared_vec %in% c((max(cand_shared_vec)-30):max(cand_shared_vec))),])
        if(max(cand_shared_vec)<shared_samples_limit){break}
        cand_idx_vec <- rems_idx
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
    fin_list[[i]]=feature_selection(diff_regs, starting_point=unlist(all_combs[i,1]), shared_samples_limit=unlist(all_combs[i,2]))
}
mat = c()
for(i in length(fin_list)){
    mat = c(mat, fin_list[[i]])

}


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


# library(doParallel)
# cl <- makeCluster(12)
# registerDoParallel(cl)
# top_ones_names = as.numeric(rownames(subset(windows_diff, win_conts>=114 & win_cases>=114)))
# remainings = windows_diff[-top_ones_names,]
# pre_set_cases = as.character(subset(windows_diff, win_conts>=114 & win_cases>=114)$win_case_ids)
# pre_set_conts = as.character(subset(windows_diff, win_conts>=114 & win_cases>=114)$win_cont_ids)
# all_rows = as.numeric(rownames(remainings))
# shared_cases=shared_conts=combinat=nocomb=c()
# for(i in c(1:3000000)){
#     if(i<=330000){z=10}
#     if(i %in% c(330000:(660000-1))){z=20}
#     if(i %in% c(660000:(1000000-1))){z=30}
#     if(i %in% c(1000000:(1330000-1))){z=40}
#     if(i %in% c(1330000:(1660000-1))){z=50}
#     if(i %in% c(1660000:(2000000-1))){z=100}
#     if(i %in% c(2000000:(2250000-1))){z=150}
#     if(i %in% c(2250000:(2500000-1))){z=200}
#     if(i %in% c(2500000:(2750000-1))){z=250}
#     if(i %in% c(2750000:(3000000-1))){z=300}
#     #print('OK1')
#     idxs = sample(all_rows, z)
#     shared_conts[i] = length(Reduce(intersect,lapply(c(as.character(remainings[idxs,]$win_cont_ids), pre_set_conts), function(x) unlist(strsplit(x, ',')))))
#     shared_cases[i] = length(Reduce(intersect,lapply(c(as.character(remainings[idxs,]$win_case_ids), pre_set_cases), function(x) unlist(strsplit(x, ',')))))
#     combinat[i] = paste(idxs,collapse=',')
#     nocomb[i] = length(idxs)
#     if(i%%100000==0){print(i)}
#     #print('OK2')
# }
