library(tidyverse)
library(data.table)
library(caret)
library(recipes)
library(pROC)
library(gbm)
library(here)


setwd(here("output", "02-artemis-delfi" ))

#These are the models we actually want
model1<-readRDS("ARTEMIS_Ensemble.rds")
model3<-readRDS("ARTEMIS_single_DELFI_SSLs_Ensemble.rds")


# SSL elements generated:
s1<-readRDS("Cov_GBM.rds")
s2<-readRDS("Ratios_ssl.rds")
s3<-readRDS("zscores_ssl.rds")
s4<-readRDS("Epi_ssl.rds")
s5<-readRDS("LINE_ssl.rds")
s6<-readRDS("LTR_ssl.rds")
s7<-readRDS("SINE_ssl.rds")
s8<-readRDS("Sat_ssl.rds")
s9<-readRDS("RNA_TE_ssl.rds")



get_coefs <- function(model) {
  orig_coefs <- coef(model$finalModel, s = model$bestTune$lambda) * (-1)
  pr <- prep(model$recipe)
  model_input <- suppressWarnings(bake(pr, new_data = model$trainingData))

  feature_means <- model_input  %>%
      select(-c(id, type)) %>%
      colMeans()
  feature_sds <- model_input %>%
      select(-c(id, type)) %>%
      as.data.frame() %>%
      summarise_all(sd)
  feature_coefs <- data.frame(features = names(feature_sds),
                            sd = as.numeric(feature_sds))
  feature_coefs <- merge(feature_coefs,
                     data.frame(features = rownames(orig_coefs),
                                orig_coefs = as.numeric(orig_coefs)),
                     by = 'features', all.x = TRUE)
  feature_coefs$scaled_coefs <- feature_coefs$orig_coefs * feature_coefs$sd
  coefs <- feature_coefs %>% filter(scaled_coefs != 0)
  return(coefs)
}


artemis <- get_coefs(model1) %>% select(features,scaled_coefs) %>% mutate(model="artemis")
joint <- get_coefs(model3) %>% select(features,scaled_coefs) %>% mutate(model="joint")

#----------------------#
#The epi coefficients should be multiplied by the epi coefficient in the artemis model
e <- get_coefs(s4) %>% select(features,scaled_coefs) %>% mutate(model="Epigenetics")
e$scaled_coefs <- e$scaled_coefs * (artemis %>% filter(features=="Epi_ssl"))$scaled_coefs

refs<-fread(here("data", "02-artemis-delfi", "Epi_Reference_Bins.csv"))

refs<-refs[map >= 0.90 & gc >= 0.30]
refs <- refs[,chr:=factor(chr, paste0("chr", c(1:22, "X")))]
setkey(refs, chr, start, end)
refs[,bin:=1:.N]
refs$bin<-paste0(refs$ref,"_",refs$bin)

e<-inner_join(e,refs %>% select(chr,start,end,ref,bin),by=c("features"="bin"))
e<-e %>% mutate(ref=if_else(ref=="H3K27me3-human_GM12878_ENCFF001SUI","H3K27me3",ref))
e<-e %>% mutate(ref=if_else(ref=="H3K36me3-human_GM12878_ENCFF001SUJ","H3K36me3",ref))
e<-e %>% mutate(ref=if_else(ref=="H3K9me3-human_GM12878_ENCFF001SUP","H3K9me3",ref))
e<-e %>% mutate(ref=if_else(ref=="states_1_5","Activation",ref))
e<-e %>% mutate(ref=if_else(ref=="states_10_13","3' Transcription",ref))
e<-e %>% mutate(ref=if_else(ref=="states_7_9","Repression",ref))
e$features<-paste0(e$ref,"_",e$chr,":",e$start,"-",e$end)
e<-e %>% select(features,scaled_coefs,model)

#----------------------#
#The Sat coefficients should be multiplied by the sat coefficient in the artemis model
s<-get_coefs(s8)%>% select(features,scaled_coefs) %>% mutate(model="Satellite")
s$scaled_coefs<-s$scaled_coefs * (artemis %>% filter(features=="Sat_ssl"))$scaled_coefs
#----------------------#
# The LINE coefficients need to be multiplied by the LINE coefficient in artemis model

l <- get_coefs(s5) %>% select(features,scaled_coefs) %>% mutate(model="LINE")
l$scaled_coefs<-l$scaled_coefs * (artemis %>% filter(features=="LINE_ssl"))$scaled_coefs
#----------------------#
a<-rbind(e,l,s) %>% mutate(model2="ARTEMIS")
a$scaled_coefs<-a$scaled_coefs * (joint %>% filter(features=="Artemis_Score"))$scaled_coefs
#----------------------#
r <- varImp(s2$finalModel) %>%
     rownames_to_column('features') %>%
     dplyr::rename(scaled_coefs = Overall) %>%
     mutate(model = 'Fragmentation Ratios') %>%
     mutate(model2 = '-')

rssl = (joint %>% filter(features=="Ratios_ssl"))$scaled_coefs
if (length(rssl) == 0){
  r$scaled_coefs<-NA
# joint model does not contain ratio features --> all coefficients are practically 0
}else{
  r$scaled_coefs<-r$scaled_coefs * rssl
}



v <- varImp(s1$finalModel) %>%
     rownames_to_column('features') %>%
     dplyr::rename(scaled_coefs = Overall) %>%
     mutate(model = 'Coverage') %>%
     mutate(model2 = '-')
v$scaled_coefs<-v$scaled_coefs * (joint %>% filter(features=="Cov_GBM"))$scaled_coefs

z<-get_coefs(s3)%>% select(features,scaled_coefs) %>% mutate(model="Aneuploidy Z-scores")
z$scaled_coefs<-z$scaled_coefs * (joint %>% filter(features=="zscores_ssl"))$scaled_coefs

z$model2 = '-'
all = rbind(z, v, a)


# table(all$model)
# Aneuploidy Z-scores  Coverage         Epigenetics             LINE           Satellite
# 38                    473                  91                  50                  32


joint
# features scaled_coefs model
# 1 Artemis_Score    1.4463142 joint
# 2       Cov_GBM    0.2759176 joint
# 3   zscores_ssl    0.1987160 joint


artemis
#      features scaled_coefs   model
# 1    Epi_ssl    1.4940520 artemis
# 2    LTR_ssl    0.6759107 artemis
# 3    Sat_ssl    0.4884257 artemis
# 4   SINE_ssl    0.4585032 artemis
# 5 RNA_TE_ssl   -0.3619743 artemis
# 6   LINE_ssl    0.3108969 artemis


output = rbind(joint, artemis) %>%
         dplyr::relocate('model', .before = 'features') %>%
         mutate(scaled_coefs = round(scaled_coefs, 2))


write.table(output, file = here('output', '02-artemis-delfi', 'model_coefficients.txt'),
            quote = FALSE, sep = '\t', row.names = FALSE)
