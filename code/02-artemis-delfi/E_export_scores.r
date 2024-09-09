library(tidyverse)
library(data.table)
library(openxlsx)

load(here('output', '02-artemis-delfi', 'scores.rda'))


path = here('..', 'Manuscript files', 'Current version', 'For publication', 'Supplementary Tables 20240710.xlsx')
base = read.xlsx(path, sheet = 2, startRow = 3)


pd = data.frame(id = base$Alternate_ID, index = seq_along(base$Alternate_ID))
pd = merge(pd, scores %>% select(id, artemis.delfi), by = 'id', all.x = TRUE) %>%
    arrange(index) %>%
    mutate(index = NULL)


write.table(pd, file = here('output', '02-artemis-delfi', 'scores_export.txt'),
            quote= FALSE, sep ='\t', row.names = FALSE)
