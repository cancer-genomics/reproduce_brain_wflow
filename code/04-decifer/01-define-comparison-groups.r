.libPaths('/dcs04/scharpf/data/pipeline-hub/pipeline-lib/R-libs/3.17.bioc-release')

library(tidyverse)
library(data.table)
library(openxlsx)

load('output/01-rbrain/metadata.rda')

metadata$alternate_id = gsub('CGCNS23P_1', 'CGCNS23P', metadata$alternate_id)


#-------------------------------------------------------#
# define the comparison groups for TFBS-expr analysis

comps = list()

#-------------------------------------------------------#

recurrent.ids = metadata %>%
                filter(training == TRUE) %>%
		filter(pathology._simplified == "Rec Grade IV" | tumor_subset_pathology == "Post-treatment") %>%
		pull(alternate_id) %>%
		sort()

full.exclude = recurrent.ids

# [1] all samples in novaseq training set
meta <- metadata %>%
        filter(training == TRUE) %>%
	select(alternate_id, type) %>%
	filter(! alternate_id %in% full.exclude)

comps[['brain']] = list(healthy = meta %>% filter(type == 'healthy') %>% pull(alternate_id),
		           cancer = meta %>% filter(type == 'cancer') %>% pull(alternate_id))
			   
#-------------------------------------------------------#
# [2] high grade samples from novaseq training set
h.rows = metadata %>%
          filter(training == TRUE) %>%
	  filter(type == 'healthy')

c.rows = metadata %>%
         filter(training == TRUE) %>%
	 filter(pathology._simplified %in% c('Grade III', 'Grade IV'))  %>%
	 filter(! alternate_id %in% full.exclude)

comps[['brain_high_grade']] = list(healthy = h.rows$alternate_id,
		           cancer = c.rows$alternate_id)

#-------------------------------------------------------#
# [2.5] gbm samples from novaseq training set
h.rows = metadata %>%
          filter(training == TRUE) %>%
	  filter(type == 'healthy')

c.rows = metadata %>%
         filter(training == TRUE) %>%
	 filter(pathology._simplified == 'Grade IV')  %>%
	 filter(! alternate_id %in% full.exclude)

comps[['brain_gbm']] = list(healthy = h.rows$alternate_id,
		           cancer = c.rows$alternate_id)

#-------------------------------------------------------#
# [3] low grade samples from novaseq training set
h.rows = metadata %>%
          filter(training == TRUE) %>%
	  filter(type == 'healthy')

c.rows = metadata %>%
         filter(training == TRUE) %>%
	 filter(pathology._simplified %in% c('Grade I', 'Grade II')) %>%
	 filter(! alternate_id %in% full.exclude)

comps[['brain_low_grade']] = list(healthy = h.rows$alternate_id,
		           cancer = c.rows$alternate_id)

#-------------------------------------------------------#

# [7] novaseq validation set gbms

h.rows = metadata %>%
          filter(training == TRUE) %>%
	  filter(type == 'healthy')

cd = c('Glioblastoma')

c.rows = metadata %>%
         filter(validation == TRUE) %>%
	 filter(tumor_subset_pathology %in% cd) %>%
	 filter(! alternate_id %in% full.exclude)

	 
comps[['brain_validation_gbm']] = list(healthy = h.rows$alternate_id,
		           cancer = c.rows$alternate_id)

#-------------------------------------------------------#
liver = read.xlsx('data/meta/Liver_Clinical_Metadata_spreadsheet_8_11.xlsx', 1)

liver <- liver %>% 
        filter(forveryclean == 1) %>% 
        mutate(type = ifelse(HCCStatus == 'Yes', 'cancer', 'healthy')) %>% 
        select(id, type)


comps[['liver']] = list(healthy = liver %>% filter(type == "healthy") %>% pull(id),
		        cancer = liver %>% filter(type == "cancer") %>% pull(id))
			

#-------------------------------------------------------#
print(lapply(comps, function(x) lapply(x, length) %>% unlist()))
saveRDS(comps, file = 'output/04-decifer/main/comparison_sets.rds')
