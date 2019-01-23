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
    #out_master_matched_ids = unique(master[master$lab_no %in% filtered_out_samples,]$RRBS.Set_ID)
    return(master[!master$lab_no %in% filtered_out_samples,])
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
        #wcoxs[z] = wilcox.test(cont_vals, case_vals)$p.value
        #vals = c(cont_vals, case_vals)
        #treats = as.factor(c(rep(0, length(cont_vals)), rep(1, length(case_vals))))
        #perm_fisher = oneway_test(vals ~ treats)
        #perms[z] = pvalue(perm_fisher)
        z=z+1
        if(z%%100000==0){print(z)}
    }
        wcoxs = ifelse(is.nan(wcoxs), 1, wcoxs)
        main_diff = cbind(main_mydiff[,c(1,2,3,4,7)], meansdiff, mediansdiff, conts, cases, means_cont, means_case, ps_and_es_matrix)
        return(main_diff)
        }
