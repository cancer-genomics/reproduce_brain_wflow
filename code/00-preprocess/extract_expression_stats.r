# extract whole blood samples from gtex (healthy individuals)

.libPaths('/dcs04/scharpf/data/pipeline-hub/pipeline-lib/R-libs/3.17-bioc-release')

library(tidyverse)
library(data.table)
library(pbapply)
library(effsize)

#----------------------------------------------#
# annotate gene ids
mapping = readRDS('data/04-decifer/gene-id-mapping.rds')

tf.ids = readRDS('data/04-decifer/all_tfs.rds')
mapping <- mapping %>% filter(gene_name %in% tf.ids)
map.array = matrix(mapping$gene_name, ncol = 1)
rownames(map.array) <- mapping$gene_id

#----------------------------------------------#
# read in tpm values
expr <- readRDS('data/expression/merged_tpm.rds')
meta <- expr$meta
mat <- expr$matrix

mat <- mat[rownames(map.array),]
#----------------------------------------------#
# annotate sample sources for aggregation
s = melt(mat)
colnames(s) = c('gene_id', 'sample_id', 'value')

s$source = meta[s$sample_id, 'source']
s$study = meta[s$sample_id, 'study']

s$gene_name = map.array[as.character(s$gene_id),1]
remove(list =c('expr'))
#----------------------------------------------#

gene.ids <- rownames(mat)
types <- meta$source %>% as.character() %>% unique() %>% sort()

master <- data.frame()
ref.ids <- meta %>% filter(source == 'gtex.whole_blood') %>% rownames()

for (type in setdiff(types, 'gtex.whole_blood')){
   type.ids <- meta %>% filter(source == type) %>% rownames() %>% as.character()
   print(type)
   for (gene in gene.ids){
       x = log2(mat[gene, ref.ids]+0.001)
       y = log2(mat[gene, type.ids]+0.001)
       master <- rbind(master, data.frame(type = type, gene = gene, 
                                          gtex.whole_blood.mean = mean(x),
                                          test.mean = mean(y),
                                          gtex.whole_blood.median = median(x),
                                          test.median = median(y),
                                          cohen.d = cohen.d(y, x)$estimate))}}




counts = meta %>% group_by(source) %>% summarize(n = n()) %>% ungroup() %>% arrange(n) %>% data.frame() %>% dplyr::rename(type = source, test.group.size = n)

master = merge(master, counts, by = 'type', all.x = TRUE)
saveRDS(master, file = 'data/04-decifer/expression-stats-by-type.rds')
