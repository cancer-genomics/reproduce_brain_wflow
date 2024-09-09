.libPaths('/dcs04/scharpf/data/pipeline-hub/pipeline-lib/R-libs/3.17-bioc-release/')

setwd('/dcs07/scharpf/data/nniknafs/delfi-brain/data/manuscript/analysis/tfbs-expr/submission-July-2024')

source('scripts/00_functions.r')
library(getopt)

opts <- getopt(matrix(c('seed', 's', 2, 'numeric'), ncol = 4, byrow = TRUE))

#---------------------------------------#

expr <- readRDS('data/meta/expression-stats-by-type.rds') 
mapping <- readRDS('data/meta/gene-id-mapping.rds')

types <-  sort(unique(expr$type)) %>% as.character()
expr <- expr %>% 
        select(gene, type, cohen.d) %>% 
	reshape2::dcast(., gene ~ type, value.var = 'cohen.d') %>% 
        merge(., mapping, by.x = 'gene', by.y = 'gene_id', all.x = TRUE) %>% 
        filter(gene_name != 'ZBED1')
expr <- expr %>% select(all_of(c('gene', 'gene_name',  types))) %>% unique() %>% dplyr::rename(tx = gene, gene = gene_name)
#---------------------------------------#

meta = readRDS('data/meta/TFBS_meta.rds')
meta$experiment = gsub(' ', '', meta$experiment)
colnames(meta)[4] = 'experiment.id'
meta$quantile = rank(meta$binding.sites) / length(meta$binding.sites)
#---------------------------------------#
print(getwd())
print(list.files('data/meta'))

comps = readRDS('results/main/comparison_sets.rds')
br.ids = comps[['brain']]
path = 'data/brain-relative-coverage/results/brain_novaseq.cohort_rel_cov.txt'
br = extract_cohort_data(path, comps['brain'])[[1]]

#---------------------------------------#
a = br$tf.cat
ids= names(a)
a = as.character(a)
names(a) = ids
br$tf.cat = a
#---------------------------------------#
set.seed(opts$seed)

btp <- get_null_sample(br)
print(dim(btp$tf.mat))
print(table(btp$tf.cat))

btp.stats <- extract_rel_cov_stats(btp)
rownames(btp.stats) <- NULL

stats <- merge(meta, btp.stats, by.x = 'experiment', by.y = 'expr', all.x = TRUE)
stats <- merge(stats, expr, by = 'gene', all.x = TRUE)

#---------------------------------------#
get_trace <- function(x){
    thresh <- seq(1,99) * 0.01
    
    o = lapply(types, function(z){

            f = lapply(thresh, 
                       function(y){get_cor(x %>% filter(quantile >= y), 
                                         'cohen.d', z)}) %>% 
                do.call(rbind, .) %>% 
                cbind(data.frame(type = z, quantile = thresh), .)
	    return(f)})
    o = do.call(rbind, o)            
            
    return(o)}


master <- get_trace(stats)

out.file <- file.path('results/null_brain/iters', paste0('run-', opts$seed, '.rds'))
saveRDS(master, file = out.file)

print(warnings())
