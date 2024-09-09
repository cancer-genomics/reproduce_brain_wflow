library(tidyverse)
library(data.table)
library(openxlsx)
library(here)

tab = read.xlsx(here("data", "Supplementary Tables 20240628.xlsx"),
                  sheet = 6, startRow = 3)

cbc = tab %>%
      dplyr::rename(id = Sample.ID,
                    pathology = Diagnosis.Group,
                    immune.cells = Immune.Cell.Population,
                    fraction = `Cell.Fraction.(%)`)

save(cbc, file = here("output", "01-rbrain", "cbc.rda"))
