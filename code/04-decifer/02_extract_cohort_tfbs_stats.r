.libPaths('/dcs04/scharpf/data/pipeline-hub/pipeline-lib/R-libs/3.17.bioc-release')

# setwd('/dcs07/scharpf/data/nniknafs/delfi-brain/data/manuscript/analysis/tfbs-expr/submission-July-2024')
source('scripts/04-decifer/00_functions.r')

meta = readRDS('data/04-decifer/TFBS_meta.rds')
meta$experiment = gsub(' ', '', meta$experiment)
colnames(meta)[4] = 'experiment.id'
meta$quantile = rank(meta$binding.sites) / length(meta$binding.sites)

# step 1: read and uniformly format the input data
comparisons = readRDS('output/04-decifer/main/comparison_sets.rds')

brain_comps = comparisons[grepl('brain', names(comparisons))]
path = 'data/04-decifer/brain_novaseq.cohort_rel_cov.txt'
part1 = extract_cohort_data(path, brain_comps)

liver_comps = comparisons[grepl('liver', names(comparisons))]
path = 'data/04-decifer/liver_cohort_rel_cov2.txt'
part2 = extract_cohort_data(path, liver_comps)

master = c(part1, part2)

saveRDS(master, file = 'output/04-decifer/main/cohort_data_tfbs_rel_cov.rds')
 
# step 2: find the differences in rel coverage between healthy and cancer
stats = pblapply(master, extract_rel_cov_stats) %>% data.table::rbindlist(., idcol = 'cohort')

# step 3: annotate with metadata
ann = merge(meta, stats, by.y = 'expr', by.x = 'experiment', all.x = TRUE) %>% arrange(cohort, abs(cohen.d))
saveRDS(ann, file = 'output/04-decifer/main/cohort_stats_tfbs_rel_cov.rds')
