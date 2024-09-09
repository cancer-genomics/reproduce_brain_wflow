### DECIFER Analysis

The following steps can be used to perform DECIFER analysis.

```
Rscript scripts/04-decifer/01-defined-comparison-groups.r

Rscript scripts/04-decifer/02_extract_cohort_tfbs_stats.r

Rscript scripts/04-decifer/04_correlate_expression_tfbs_rel_cov.r

# the following steps should be run on a high performance compute cluster
# assuming SLURM scheduler. Each instance of this array job 
# calls the script 05_brain_null.r with a specified random seed.
sbatch --array=1-1000 scripts/04-decifer/brain_null_sample.sh

Rscript scripts/04-decifer/06_brain_null_summarize.r

# Summarizing the analysis results in an excel workbook
Rscript scripts/04-decifer/07_export_tables.r
```
