.libPaths('/dcl01/scharpf/data/pipeline-hub/pipeline-lib/R-libs/3.10-bioc-release/')
source('scripts/00_functions.r')

library(tidyverse)
library(data.table)


#---------------------------------------#
mapping <- readRDS('data/meta/gene-id-mapping.rds')

expr <- readRDS('data/meta/expression-stats-by-type.rds') 
types <-  sort(unique(expr$type)) %>% as.character()
# types <- setdiff(types, c('Adult tumor periphery:Astrocyte', "Adult temporal lobe:Neuron"))

expr <- expr %>% 
        select(gene, type, cohen.d) %>% 
	reshape2::dcast(., gene ~ type, value.var = 'cohen.d') %>% 
        merge(., mapping, by.x = 'gene', by.y = 'gene_id', all.x = TRUE) %>% 
        filter(gene_name != 'ZBED1')

#---------------------------------------#
meta = readRDS('data/meta/TFBS_meta.rds')
meta$experiment = gsub(' ', '', meta$experiment)
colnames(meta)[4] = 'experiment.id'
meta$quantile = rank(meta$binding.sites) / length(meta$binding.sites)
#---------------------------------------#
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
br.stats <- extract_rel_cov_stats(br)
rownames(br.stats) <- NULL
br.stats <- merge(meta, br.stats, by.x = 'experiment', by.y = 'expr', all.x = TRUE)
br.stats <- merge(br.stats, expr, by.x = 'gene', by.y = 'gene_name', all.x = TRUE)

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
#---------------------------------------#
br.trace = get_trace(br.stats)
rownames(br.trace) = NULL
br.trace <- br.trace %>% unite('key', type:quantile, sep = '_', remove = FALSE)
br.trace.split <- split(br.trace, f = br.trace$key)
#---------------------------------------#

files = paste('results/null_brain/iters/run-', seq(1,1000), '.rds', sep = '')
iter.count = length(files)
f = pblapply(files, readRDS)
null = rbindlist(f, idcol = 'iter') %>% 
       unite('key', type:quantile, sep = '_', remove = FALSE)

null.split = split(null, f = null$key)
#---------------------------------------#
r.pval = pblapply(names(null.split), 
                function(k) length(which(null.split[[k]]$r <= br.trace.split[[k]]$r)) / iter.count)

rho.pval = pblapply(names(null.split), 
                function(k) length(which(null.split[[k]]$rho <= br.trace.split[[k]]$rho)) / iter.count)

m = data.frame(key = names(null.split), r.pval = r.pval %>% unlist(), rho.pval = rho.pval %>% unlist())

br.trace <- merge(br.trace, m, by = 'key', all.x = TRUE)

saveRDS(br.trace, file = 'results/null_brain/summary/brain_cor_pvals.rds')
plot.data = list(null = null.split, obs = br.trace.split)
saveRDS(plot.data, file = 'results/null_brain/summary/brain_with_null_plotdata.rds')

null = plot.data$null

r     <- pblapply(null, function(x) x %>% pull(r) %>% quantile(., prob =c(0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90,  0.95, 0.99))) %>%
     pblapply(., function(x) data.table(t(x))) %>%
     do.call(rbind, .) %>%
     mutate(key = names(null)) %>%
     merge( rbindlist(null) %>% select(key, type, quantile) %>% unique(), . , by = 'key')

rho <- pblapply(null, function(x) x %>% pull(rho) %>% quantile(., prob =c(0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90,  0.95, 0.99))) %>%
     pblapply(., function(x) data.table(t(x))) %>%
     do.call(rbind, .) %>%
     mutate(key = names(null)) %>%
     merge( rbindlist(null) %>% select(key, type, quantile) %>% unique(), . , by = 'key')

out = list(obs = plot.data$obs, r = r, rho = rho)
saveRDS(out, file = 'results/null_brain/summary/brain_with_null_plotdata_stats.rds')


