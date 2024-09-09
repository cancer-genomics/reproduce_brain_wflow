library(tidyverse)
library(data.table)
library(openxlsx)
library(here)

cv = fread(here("output", "02-artemis-delfi", "Cross_Validation_scores.csv"))
heldout = fread(here("output", "02-artemis-delfi", "Test_set_scores.csv"))

cv_scores = cv %>%
            select(-V1) %>%
            filter(model %in% c("ARTEMIS_single_DELFI_SSLs_Ensemble", "zscores_ssl")) %>%
            select(id, score, model) %>%
            pivot_wider(id_cols = c('id'), names_from ='model', values_from = 'score') %>%
            dplyr::rename(artemis.delfi = ARTEMIS_single_DELFI_SSLs_Ensemble, copynumber = zscores_ssl) %>%
            mutate(partition = 'training') %>%
            data.frame()

heldout_scores = heldout %>%
  select(-V1) %>%
  filter(model %in% c("ARTEMIS_single_DELFI_SSLs_Ensemble", "zscores_ssl")) %>%
  select(id, score, model) %>%
  pivot_wider(id_cols = c('id'), names_from ='model', values_from = 'score') %>%
  dplyr::rename(artemis.delfi = ARTEMIS_single_DELFI_SSLs_Ensemble, copynumber = zscores_ssl) %>%
  mutate(partition = 'heldout') %>%
  data.frame()

scores = rbind(cv_scores, heldout_scores)

scores$id[scores$id == "CGCNS23P"] = 'CGCNS23P_1'

save(scores, file = here("output", "02-artemis-delfi", "scores.rda"))
