library(tidyverse)
library(data.table)
library(openxlsx)
library(here)

S1 = read.xlsx(here("data", "Supplementary Tables 20240628.xlsx"),
               sheet = 1, startRow = 3)

S2 = read.xlsx(here("data", "Supplementary Tables 20240628.xlsx"),
               sheet = 2, startRow = 3)


colnames(S1)[3] = 'Patient.ID'
S1 = S1 %>% select(-2) %>% filter(grepl('CG', Complete.Sample.ID))


base = S1 %>%
  select(Complete.Sample.ID, Clinical.Condition, Age, Gender, Pathology.Classification, Detailed.Pathology, MRI.Enhancement, Ki67, `Sample.Source*`,
        `Tumor.Volume.(cm3)`) %>%
  dplyr::rename(alternate_id = Complete.Sample.ID,
                age = Age,
                gender = Gender,
                pathology._simplified = Pathology.Classification,
                tumor_subset_pathology = Detailed.Pathology,
                mri.enhancement = MRI.Enhancement,
                ki67.simplified = Ki67,
                source = `Sample.Source*`,
                volume = `Tumor.Volume.(cm3)`) %>%
  mutate(tumor.lab.id  = alternate_id,
         type = ifelse(Clinical.Condition == "Healthy", 'healthy', 'cancer')) %>%
  select(alternate_id, type, age, gender, pathology._simplified, tumor_subset_pathology, mri.enhancement, ki67.simplified, source, volume)

base$patient_id = base %>%
                  mutate(p = gsub('P1', 'P', alternate_id)) %>%
                  separate(p, into = c('root', 'suffix'), '_') %>%
                  mutate(r = gsub('P$', '', root)) %>%
                  pull(r)


S2 = S2 %>%
     select(Alternate_ID, `Discovery.Cohort.(Training)`, `Discovery.Cohort.(Held-out)`, `Validation.Cohort`) %>%
     dplyr::rename(alternate_id = Alternate_ID) %>%
     mutate(training = ifelse(`Discovery.Cohort.(Training)` == 'YES', TRUE, FALSE),
         heldout = ifelse(`Discovery.Cohort.(Held-out)` == 'YES', TRUE, FALSE),
         validation = ifelse(`Validation.Cohort` == 'YES', TRUE, FALSE)) %>%
     select(alternate_id, training, heldout, validation)


data = merge(base, S2, by = 'alternate_id', all.x =TRUE) %>%
       arrange(alternate_id) %>%
       mutate(index = seq_along(alternate_id), drop = FALSE, exclude= FALSE) %>%
       dplyr::relocate(index, .before = 'alternate_id')

metadata = data %>%
         select(index, alternate_id, patient_id, age, gender, type, training, validation, heldout,exclude, drop,
         pathology._simplified,tumor_subset_pathology, mri.enhancement, ki67.simplified, volume, source)

save(metadata, file = here("output", "01-rbrain", "metadata.rda"))
