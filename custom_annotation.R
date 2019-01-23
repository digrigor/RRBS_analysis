#This script annotates CpG sites with annotatr and converts the annotation file to
#a more readable format

chf <- function(df0){
d=''
for(j in 1:nrow(df0)){
  b <- paste(lapply(df0[j, c(20,21,22,24)], as.character), collapse=':')
  a <- paste(lapply(df0[j, c(24,23)], as.character), collapse=':')
  c <- paste(a,b,sep=':')
  d <- paste(c,d, sep=';')
  }
d <- paste(unique(unlist(strsplit(unlist(d),';'))), collapse=';')
return(d)
}

chf2 <- function(df2){
d=''
for(j in 1:nrow(df2)){
  b <- paste(lapply(df2[j, c(20,21,22,24)], as.character), collapse=':')
  d <- paste(b,d, sep=';')
  }
d <- paste(unique(unlist(strsplit(unlist(d),';'))), collapse=';')
return(d)
}

chf3 <- function(df3){
d=''
for(j in 1:nrow(df3)){
  b <- paste(lapply(df3[j,c(20,21,22,24)], as.character), collapse=':')
  a <- paste(lapply(df3[j, c(23,24)], as.character), collapse=':')
  c <- paste(a,b,sep=':')
  if(d==''){
    d <- paste(c,d, sep='')
    } else {
    d <- paste(c,d, sep=';')
    }
  }
  d <- paste(unique(unlist(strsplit(unlist(d),';'))), collapse=';')
  return(d)
}

chf4 <- function(df4, i){
d=''
for(j in 1:nrow(df4)){
  a <- unlist(strsplit(as.character(df4[j,25]),'Hmec-'))[2]
  b <- paste(lapply(df4[j,c(20,21,22,24)], as.character), collapse=':')
  c <- paste(a,b,sep=':')
  if(d==''){
    d <- paste(c,d, sep='')
    } else {
    d <- paste(c,d, sep=';')
    }
  }
  d <- paste(unique(unlist(strsplit(unlist(d),';'))), collapse=';')
  return(d)
}


transpose_annotation <- function(df){

  count = 1
  finlist = list()
  lenids <- length(unique(df$cpg_id))
  lenannots <- length(unique(df$annot.type)) - sum(startsWith(unique(df$annot.type), 'hg19_chromatin')) + 1

  for(d in c(1:lenids)) {

    chrout=''
    templist = rep('', lenannots)
    rows = df[df$cpg_id == unique(df$cpg_id)[d],]
    cpglist <- unname(unlist(lapply(rows[1,][c(15,1:14)], as.character)))

    for(type in unique(rows$annot.type)){

      rowstype = rows[rows$annot.type==type,]
      if(startsWith(type, 'hg19_chromatin') == TRUE){

        func <- 'chf4'
        tempout <- do.call(func, list(rowstype))
        chrout <- paste(chrout, tempout, sep='')

      } else {

        func <- as.character(anndict$V2[anndict$V1==type])
        ind <- as.numeric(anndict$V3[anndict$V1==type])
        out <- do.call(func, list(rowstype))
        templist[ind]=out

      }

      templist[17] <- chrout

    }

    rowlist <- append(cpglist, templist)
    finlist <- unlist(append(finlist, rowlist))
    count <- count+1
    if(count%%1000==0){print(count)}
  }
return(finlist)
}
finmatrix = matrix(finlist, byrow=TRUE, ncol=28)
fdf = as.data.frame(finmatrix)
colnames(fdf) <- append(colnames(df)[c(11,1:10)], anndict$V4)
return(fdf)
return(fdf)
}

#Load required library
library(annotatr)

#Set the annotation path
annotation_dir <- '/home/dio01/Methyl_Annotations/'
anndict <- readRDS(paste(annotation_dir,'anndict.rds',sep=''))

#Prepare the annotation
TFBS_file <- paste(annotation_dir, 'wgEncodeRegTfbsClusteredV3.bed.gz', sep='')
Dnase_file <- paste(annotation_dir, 'wgEncodeRegDnaseClusteredV3.bed.gz', sep='')
read_annotations(con=TFBS_file, genome='hg19', name='TfbsClusteredV3',format='bed', extraCols=c(gene_id="character", symbol="numeric"))
read_annotations(con=Dnase_file, genome='hg19', name='DnaseClusteredV3',format='bed')
annots <- builtin_annotations()[c(4,14,24,34,44,54,64,74,84,94,104,112,120,128,136,251,255,261,269)]
annots <- append(annots, 'hg19_Hmec-chromatin')
annots <- append(annots, annotatr_cache$list_env())
annots <- annots[-17]
annots <- annots[-c(8,9)]

annotations = build_annotations(genome = 'hg19', annotations = annots)

