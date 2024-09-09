
library(tidyverse)
library(pbapply)
library(openxlsx)
library(here)
#----------------------------------------------#
s = readRDS(here('data', '04-decifer', 'TFBS_meta.rds'))

out = s %>%
      select(experiment, expr, gene, specimen, binding.sites) %>%
      dplyr::rename(`Remap2020 Database ID` = experiment,
                    `Gene Symbol` = gene,
                    `Experiment ID` = expr,
                    `Specimen` = specimen,
                    `Number of Binding Sites` = binding.sites)

wb = createWorkbook()
addWorksheet(wb, 'Table S7')
writeData(wb, 'Table S7', out, startRow = 3, keepNA = TRUE)
#----------------------------------------------#
s = readRDS(here('output', '04-decifer', 'main', 'cohort_data_tfbs_rel_cov.rds'))

# elements to export are:
# s$brain, s$brain_validation_gbm, s$liver?

validation = s$brain_validation_gbm$tf.mat
discovery = s$brain$tf.mat

validation_gbm_samples = setdiff(colnames(validation), colnames(discovery))
out = cbind(discovery, validation[,validation_gbm_samples] )

addWorksheet(wb, 'Table S8')
writeData(wb, 'Table S8', out, rowNames = TRUE, startRow = 3, keepNA = TRUE)
#----------------------------------------------#
s = readRDS(here('data', '04-decifer', 'expression-stats-by-type.rds'))

s$type = gsub('lymphocytes', 'Lymphocytes', s$type)
s$type = gsub('GSE60424.CD4', 'CD4 T-Cells', s$type)
s$type = gsub('GSE60424.CD8', 'CD8 T-Cells', s$type)
s$type = gsub('GSE60424.monocytes', 'Monocytes', s$type)
s$type = gsub('GSE60424.neutrophils', 'Neutrophils', s$type)

out = s %>%
      select(type, gene, gtex.whole_blood.mean, test.mean, cohen.d) %>%
      dplyr::rename(`Tumor/Cell Type` = type, `Gene ID` = gene,
                    `Mean Expression in Whole Blood*` = gtex.whole_blood.mean,
                    `Mean Expression in Tumor/Cell Type*` = test.mean,
                    `Effect Size (Cohen's D)` = cohen.d)

addWorksheet(wb, 'Table S9')
writeData(wb, 'Table S9', out, startRow = 3, keepNA = TRUE)
#----------------------------------------------#
s = readRDS(here('output', '04-decifer', 'main',
                 'cohort_rel_cov_tpm_merged_cor_complete.rds'))

s$type = gsub('lymphocytes', 'Lymphocytes', s$type)
s$type = gsub('GSE60424.CD4', 'CD4 T-Cells', s$type)
s$type = gsub('GSE60424.CD8', 'CD8 T-Cells', s$type)
s$type = gsub('GSE60424.monocytes', 'Monocytes', s$type)
s$type = gsub('GSE60424.neutrophils', 'Neutrophils', s$type)

s = s %>%
    filter(cohort %in% c('brain', 'brain_validation_gbm', 'liver')) %>%
    select(cohort, type, quantile, r) %>%
    dplyr::rename(`Cohort ID` = cohort, `Tumor/Cell Type` = type,
                  `TFBS Quantile` = quantile, `Pearson's R` = r)

s$`Cohort ID` = gsub('brain_validation_gbm', 'Brain Tumor Validation Set (GBM Samples)', s$`Cohort ID`)
s$`Cohort ID` = gsub('brain', 'Brain Tumor Discovery Set', s$`Cohort ID`)
s$`Cohort ID` = gsub('liver', 'Liver Cancer Cohort*', s$`Cohort ID`)

p1 = s %>%
     filter(`Cohort ID` != "Liver Cancer Cohort*")

p2 = s %>%
     filter(`Cohort ID` == "Liver Cancer Cohort*") %>%
     filter(`Tumor/Cell Type` %in% c('LIHC', 'COAD', 'GBM', 'LGG', 'LUAD'))

s = rbind(p1, p2)

addWorksheet(wb, 'Table S10')
writeData(wb, 'Table S10', s, startRow = 3, keepNA = TRUE)
#----------------------------------------------#
saveWorkbook(wb,
             file = here('output', '04-decifer', 'export', 'decifer-supplementary-tables.xlsx'),
             overwrite = TRUE)

