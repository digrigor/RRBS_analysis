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


