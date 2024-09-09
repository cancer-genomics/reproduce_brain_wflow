# helper functions derived here
.libPaths('/dcs04/scharpf/data/pipeline-hub/pipeline-lib/R-libs/3.17.bioc-release')

library(tidyverse)
library(openxlsx)
library(pbapply)
library(pROC)
library(effsize)

#---------------------------------------#
extract_comparison_data <- function(path, labels){
    tfs <- data.table::fread(path)

    if (colnames(tfs)[1] == 'id'){
       colnames(tfs)[1] = 'V1'
       }
       

    tf.mat <- tfs %>% data.frame() %>% dplyr::rename(id = V1) %>% column_to_rownames('id') %>%   t() %>% as.matrix()
    rownames(tf.mat) <- colnames(tfs)[-1]

    meta = data.frame(id = c(labels$healthy, labels$cancer),
    	              type = c(rep('healthy', length(labels$healthy)), rep('cancer', length(labels$cancer))))
		      
    cols <- data.frame(id  = seq(ncol(tf.mat)), tumor.pgdx.id = colnames(tf.mat))
    cols <- merge(cols, meta, by.x = 'tumor.pgdx.id', by.y = 'id', all.x = TRUE)
    cols <- cols %>% filter(! is.na(type)) %>% arrange(id)

    tf.mat <- tf.mat[,as.character(cols$tumor.pgdx.id)]
    tf.cat <- cols$type
    names(tf.cat) <- cols$tumor.pgdx.id

    out <- list(tf.mat = tf.mat, tf.cat = tf.cat)
    return(out)
}
#---------------------------------------#
extract_cohort_data = function(path, comparisons){
    out = list()
    for (id in names(comparisons)){
         out[[id]] = extract_comparison_data(path, comparisons[[id]])
    }
    return(out)
}
#---------------------------------------#
extract_template_data <- function(x){
    cwd <- getwd()
    setwd('/dcs04/scharpf/data/nniknafs/delfi-brain/analysis/tfbs-expr/manuscript')

#     tf.mat <- tf.mat[,as.character(cols$tumor.pgdx.id)]
#     tf.cat <- cols$type
#    names(tf.cat) <- cols$tumor.pgdx.id

    setwd(cwd)
#    out <- list(tf.mat = tf.mat, tf.cat = tf.cat)
#    return(out)
}

#---------------------------------------#
extract_rel_cov_stats <- function(cohort_data){

    tf.mat = cohort_data$tf.mat
    tf.cat = cohort_data$tf.cat

    train_data <- tf.mat
    train_labels <- tf.cat
    cancer.ind = which(train_labels == 'cancer')
    healthy.ind = which(train_labels == 'healthy')

    aucs <- lapply(rownames(train_data),
            function(x)  auc(pROC::roc(train_labels, train_data[x,], levels = c('healthy', 'cancer'), quiet = TRUE))) %>% 
            unlist()
    names(aucs) <- rownames(train_data)


    d <- lapply(rownames(train_data),
         function(x) cohen.d(train_data[x,cancer.ind], train_data[x,healthy.ind])$estimate) %>% unlist()
    names(d) <- rownames(train_data)

    base <- data.frame(expr  = rownames(train_data), auc = aucs, cohen.d = d)
    return(base)
}
#---------------------------------------#
get_cor = function(x, a, b) data.frame(r = cor.test(x[,a], x[,b], method = 'pearson')$estimate, rho = cor.test(x[,a], x[,b], method = 'spearman')$estimate)
#---------------------------------------#
get_bootstrap_sample <- function(cohort_data){
    tf.mat = cohort_data$tf.mat
    tf.cat = cohort_data$tf.cat
    
    cancer.ind = which(tf.cat == 'cancer')
    cancer.mat = tf.mat[,cancer.ind]
    cancer.cat = tf.cat[cancer.ind]
    
    healthy.ind = which(tf.cat == 'healthy')
    healthy.mat = tf.mat[,healthy.ind]
    healthy.cat = tf.cat[healthy.ind]


    ind = seq_along(cancer.cat)
    
    selected <- sample(ind, replace = TRUE) %>% sort()

    selected.tf.cat = cancer.cat[selected]
    selected.tf.mat = cancer.mat[,selected]

    ids <- paste('s', seq_along(selected.tf.cat), sep ='')
    names(selected.tf.cat) = ids
    colnames(selected.tf.mat) = ids

    out.mat = cbind(selected.tf.mat, healthy.mat)
    out.cat = c(selected.tf.cat, healthy.cat)
    
    out = list(tf.cat = out.cat, tf.mat = out.mat)
    return(out)
}
#---------------------------------------#
get_permute_sample <- function(cohort_data){
    tf.mat = cohort_data$tf.mat

    tf.cat = cohort_data$tf.cat
    ind = sample(seq_along(tf.cat), replace = FALSE)
    out.cat = tf.cat[ind]
    names(out.cat) = names(tf.cat)

    out = list(tf.cat = out.cat, tf.mat = tf.mat)
    return(out)
}
#---------------------------------------#
get_null_sample <- function(cohort_data){
    # start by only looking at healthy samples
    healthy.ids = names(cohort_data$tf.cat[cohort_data$tf.cat == 'healthy']) %>% as.character()
    # randomly assign a subset to be cancer for generating a null
    # keep the cancer/ healthy ratio the same as the original data set
    
    counts = table(cohort_data$tf.cat)
    cancer.frac = counts['cancer'] / sum(counts)
    cancer.count = (cancer.frac * length(healthy.ids)) %>% floor()

    tf.mat = cohort_data$tf.mat[,healthy.ids]
    cancer.ind = sort(sample(seq_along(healthy.ids), cancer.count, replace = FALSE))

    tf.cat = rep('healthy', ncol(tf.mat))
    tf.cat[cancer.ind] = 'cancer'
    names(tf.cat) = healthy.ids

    out = list(tf.cat = tf.cat, tf.mat = tf.mat)
    return(out)
}
#---------------------------------------#
get_subset_sample <- function(cohort_data, n){
    tf.mat = cohort_data$tf.mat
    tf.cat = cohort_data$tf.cat
    
    cancer.ind = which(tf.cat == 'cancer')
    cancer.mat = tf.mat[,cancer.ind]
    cancer.cat = tf.cat[cancer.ind]
    
    healthy.ind = which(tf.cat == 'healthy')
    healthy.mat = tf.mat[,healthy.ind]
    healthy.cat = tf.cat[healthy.ind]


    ind = seq_along(cancer.cat)
    
    selected <- sample(ind, size = n, replace = FALSE) %>% sort()

    selected.tf.cat = cancer.cat[selected]
    selected.tf.mat = cancer.mat[,selected]

    ids <- paste('s', seq_along(selected.tf.cat), sep ='')
    names(selected.tf.cat) = ids
    colnames(selected.tf.mat) = ids

    out.mat = cbind(selected.tf.mat, healthy.mat)
    out.cat = c(selected.tf.cat, healthy.cat)
    
    out = list(tf.cat = out.cat, tf.mat = out.mat)
    return(out)
}
#---------------------------------------#
get_test_data <- function(n){
   tf.mat = matrix(sort(rep(seq(n), 4)), ncol = n, byrow = FALSE)
   tf.cat = c(rep('healthy', n/2), rep('cancer', n/2))
   colnames(tf.mat) = letters[1:n]
   rownames(tf.mat) = c('W', 'X', 'Y','Z')
   names(tf.cat) = letters[1:n]

   cohort_data= list(tf.mat = tf.mat, tf.cat = tf.cat)
   extract_rel_cov_stats(cohort_data)
   return(cohort_data)
}
#---------------------------------------#
rep.row<-function(x,n){
   matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
   matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
#---------------------------------------#
extract_selected_ids <- function(x, ids){
   tf.mat = x$tf.mat
   tf.cat = x$tf.cat

   samples = colnames(tf.mat)
   samples = samples[samples %in% ids]

   tf.mat = tf.mat[,samples]
   tf.cat = tf.cat[samples]

   out = list(tf.mat = tf.mat, tf.cat = tf.cat)
   return(out)
}
