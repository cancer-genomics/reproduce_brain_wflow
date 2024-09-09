.libPaths('/dcs04/scharpf/data/pipeline-hub/pipeline-lib/R-libs/3.17.bioc-release')
library(tidyverse)
library(data.table)

# setwd('/dcs07/scharpf/data/nniknafs/delfi-brain/data/manuscript/analysis/tfbs-expr/submission-July-2024')
source('code/04-decifer/00_functions.r')

cov <- readRDS('output/04-decifer/main/cohort_stats_tfbs_rel_cov.rds')
expr <- readRDS('data/04-decifer/expression-stats-by-type.rds') 
mapping <- readRDS('data/04-decifer/gene-id-mapping.rds')

types = c('GBM',  'LGG', 'AML', 'BRCA', 'BLCA', 'COAD', 'KIRC', 'LIHC','LUSC',  'LUAD',  'PRAD',  'STAD', 'OV',
          'GSE60424.neutrophils', 'GSE60424.CD8', 'GSE60424.CD4', 'GSE60424.monocytes',
	  'Adult tumor core:Astrocyte', 'Adult temporal lobe:Astrocyte', 'Adult temporal lobe:Endothelial', 'Adult temporal lobe:Myeloid',
          'Adult temporal lobe:Oligodendrocyte', 'Adult temporal lobe:Whole cortex',
	  'lymphocytes')
	    
expr.used = expr %>% filter(type %in% types)
# saveRDS(expr.used, 'data/meta/expression-stats-by-type.rds')
expr = expr.used


types <-  sort(unique(expr$type)) %>% as.character()
# types <- setdiff(types, c('Adult tumor periphery:Astrocyte', "Adult temporal lobe:Neuron"))

expr <- expr %>% 
        select(gene, type, cohen.d) %>% 
	reshape2::dcast(., gene ~ type, value.var = 'cohen.d') %>% 
        merge(., mapping, by.x = 'gene', by.y = 'gene_id', all.x = TRUE) %>% 
        filter(gene_name != 'ZBED1')


cov <- merge(cov, expr, by.x = 'gene', by.y = 'gene_name', all.x = TRUE)
saveRDS(cov, file = 'output/04-decifer/main/cohort_rel_cov_tpm_merged_matrix_complete.rds')


cov <- split(cov, cov$cohort)
thresh <- seq(1,99) * 0.01
master <- data.frame()

for (cohort in names(cov)){
    for (tumor.type in types){
       print(paste0(cohort, ' - ', tumor.type))
       print(date())
       f = lapply(thresh, 
           function(x){get_cor(cov[[cohort]] %>% filter(quantile >= x), 
                                         'cohen.d', tumor.type)}) %>% 
           do.call(rbind, .) %>% 
           cbind(data.frame(cohort = cohort, type = tumor.type, quantile = thresh), .)
        master <- rbind(master, f)

}

}

saveRDS(master, file = 'output/04-decifer/main/cohort_rel_cov_tpm_merged_cor_complete.rds')
