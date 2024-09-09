library(tidyverse)
library(data.table)
library(openxlsx)
library(here)

ddpcr = read.xlsx(here("data", "Supplementary Tables 20240628.xlsx"),
               sheet = 3, startRow = 3)

ts = read.xlsx(here("data", "Supplementary Tables 20240628.xlsx"),
               sheet = 4, startRow = 3)


ddpcr = ddpcr %>%
        mutate(method = 'ddPCR')  %>%
        filter(Sample.ID != 'CGCNS49P1') %>%
        mutate(gene = 'IDH1') %>%
        dplyr::rename(id = Sample.ID)
ddpcr$maf = ddpcr$`Mutant.Allele.Fraction.(%)`
ddpcr$maf[ddpcr$`Mutant.Allele.Fraction.(%)` == "<LOD"] = 1e-6
ddpcr$maf = as.numeric(ddpcr$maf)

ts_pos = ts %>%
     filter(`Mutation.Designation` == "Tumor Mutation") %>%
     select(Sample, Gene.Symbol, `%.Mutant.Reads`) %>%
     dplyr::rename(id = Sample, gene = Gene.Symbol, maf = `%.Mutant.Reads`) %>%
     mutate(method = 'Targeted Sequencing')

pos_ids = ts_pos$id
tested_ids = ts %>% filter(`Mutation.Designation` != "Exclude") %>% pull(Sample)
neg_ids = setdiff(tested_ids, pos_ids)

ts_neg = data.frame(id = neg_ids) %>%
     mutate(gene = "-", maf = 1e-6, method = "Targeted Sequencing")

ts = rbind(ts_neg, ts_pos)

plasma_maf = rbind(ts %>% select(id, gene, maf, method),
             ddpcr %>% select(id, gene, maf, method))

save(plasma_maf, file = here("output", "01-rbrain", "plasma_maf.rda"))

