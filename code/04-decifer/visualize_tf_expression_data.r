# module load conda_R/4.2.x
.libPaths('/dcl01/scharpf/data/pipeline-hub/pipeline-lib/R-libs/3.14-bioc-release/')

library(tidyverse)
library(pbapply)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(scales)
library(getopt)
#----------------------------------------------#
opts = list()
opts$quantile.threshold = 0.1

print(opts)
quantile.threshold = opts$quantile.threshold

clust.method = 'spearman'
output.figure.file = file.path('figures', paste0('TCGA-GTEX-expression-heatmap-tpm-log2-var-',quantile.threshold, '-cluster-', clust.method, '.pdf'))
#----------------------------------------------#
# annotate gene ids
mapping = readRDS('data/meta/gene-id-mapping.rds')

map.array = matrix(mapping$gene_name, ncol = 1)
rownames(map.array) <- mapping$gene_id
#----------------------------------------------#
tfbs.meta <- readRDS('data/meta/TFBS_meta.rds')
tfbs.meta$quantile  = rank(tfbs.meta$binding.sites)/ nrow(tfbs.meta)

tf.gene.ids <- tfbs.meta %>% pull(gene) %>% unique() %>% as.character()
map.array = map.array[map.array[,1] %in% tf.gene.ids,]
#----------------------------------------------#
types = c('gtex.whole_blood', 'AML', 'DLBC', 'GSE60424.CD4', 'GSE60424.CD8',  'GSE60424.monocytes', 'GSE60424.neutrophils',
       	    		      	     'LGG', 'GBM', 'Adult tumor core:Astrocyte', 'BLCA', 'BRCA', 'COAD' ,'KIRC', 'LIHC', 'LUAD', 'LUSC', 'OV' ,'PRAD' ,'STAD')
				     
# read in tpm values
expr <- readRDS('/dcs05/scharpf/data/nniknafs/delfi-brain/tfbs-expr/bulk-immune-rerun/merged/data/merged_tpm.rds')
meta <- expr$meta
mat <- expr$matrix

meta= meta %>% filter(source %in% types)

write.table(meta, 'export/expression_data_set_metadata.txt', quote = FALSE, sep ='\t', row.names = TRUE)
#----------------------------------------------#
# rows = meta %>% filter(study %in% c('tcga', 'gtex'))
# mat = mat[names(map.array),rownames(rows)]

rows = meta %>% filter(source %in% c('gtex.whole_blood', 'AML', 'DLBC', 'GSE60424.CD4', 'GSE60424.CD8',  'GSE60424.monocytes', 'GSE60424.neutrophils',
       	    		      	     'LGG', 'GBM', 'Adult tumor core:Astrocyte', 'BLCA', 'BRCA', 'COAD' ,'KIRC', 'LIHC', 'LUAD', 'LUSC', 'OV' ,'PRAD' ,'STAD'))

rows$source = gsub('GSE60424.CD8', 'CD8 T-Cells', rows$source)
rows$source = gsub('GSE60424.CD4', 'CD4 T-Cells', rows$source)
rows$source = gsub('GSE60424.monocytes', 'Monocytes', rows$source)
rows$source = gsub('GSE60424.neutrophils', 'Neutrophils', rows$source)
rows$source = gsub('Adult tumor core:Astrocyte', 'Tumor Core Astrocytes', rows$source)
rows$source = gsub('gtex.whole_blood', 'Whole Blood', rows$source)

mat = mat[names(map.array),rownames(rows)]
remove(list = c('expr'))
#----------------------------------------------#

# transform tpm values to log scale
mat = log2(mat + 0.001)

# calculate standard deviation of gene expression 
# across cohort of samples, and select those within a variance tier or above
row.sd = apply(mat, 1, sd)
r = rank(row.sd) / length(row.sd)


# select TFS within the top 20% of variance across the entire TCGA+GTEX data set
ind = as.numeric(which(r >= 1 - quantile.threshold))
keep.tfs = names(r[ind])



# select n samples per tumor type
set.seed(1)
sr = split(rows, as.character(rows$source))
n.per.type = 20
keep.samples = pblapply(sr, 
                        function(x){
                            if (nrow(x) < n.per.type){
                               return(rownames(x))}else{
                                sample(rownames(x), size= n.per.type, replace= FALSE)}})

selected.types = c('Whole Blood', 'AML', 'DLBC','CD4 T-Cells', 'CD8 T-Cells',  'Monocytes', 'Neutrophils',
       	    		      	     'LGG', 'GBM', 'Tumor Core Astrocytes', 'BLCA', 'BRCA', 'COAD' ,'KIRC', 'LIHC', 'LUAD', 'LUSC', 'OV' ,'PRAD' ,'STAD')

keep.samples <- keep.samples[selected.types]
keep.samples <-  keep.samples %>% unlist() %>% unique() %>% sort()


pd <- mat[keep.tfs, keep.samples]
# scale rows of matrix
pd <- t(scale(t(pd)))

annot.df = rows[keep.samples,]
annot.df$study = NULL

annot.df$source  = factor(as.character(annot.df$source), levels = selected.types)


annot.df = annot.df %>% arrange(source)
pd = pd[,rownames(annot.df)]



# color scale for heatmap body
myCol <- colorRampPalette(c('dodgerblue', 'black', 'yellow'))(100)
myBreaks <- seq(-3, 3, length.out = 100)


cols = c("#980000", "#cc0000","#e06666","#e69138","#f6b26b",
 "#f9cb9c","#fce5cd","#93c47d","#6aa84f","#38761d",
 "#ffe599","#f1c232","#a4c2f4","#6d9eeb","#1155cc",
 "#ead1dc","#c27ba0","#8e7cc3","#351c75","#cccccc")

names(cols) = c("Whole Blood", "AML", "DLBC", "CD8 T-Cells", "CD4 T-Cells", "Monocytes", "Neutrophils", "LGG", "GBM", "Tumor Core Astrocytes", "BLCA", "KIRC", "COAD", "STAD", "LIHC", "LUAD", "LUSC", "BRCA", "OV", "PRAD")



colAnn <- HeatmapAnnotation(df= annot.df, which = 'col', na_col = 'white', col = list(source = cols), annotation_legend_param = list(source = list(title = 'Tumor/Cell Type')))



# all samples
set.seed(100)
pdf(output.figure.file, width = 12, height = 10)
Heatmap(pd, 
        col = colorRamp2(myBreaks, myCol),
        cluster_rows = TRUE, 
        cluster_columns = TRUE, 
        show_row_dend = TRUE, 
        show_column_dend = TRUE,
	show_column_names = FALSE,
	show_row_names = FALSE,
	top_annotation = colAnn,
        use_raster = FALSE,
	clustering_distance_rows = clust.method,
	clustering_distance_columns = clust.method, heatmap_legend_param = list(title = 'Expression\nLog2(TPM + 0.001)'))

dev.off()


out.data = list(keep.tfs = keep.tfs, keep.samples = keep.samples, annot = annot.df  %>% rownames_to_column('id') %>% filter(id %in% keep.samples),
                expression.mat = pd)
saveRDS(out.data, file = 'export/expression_heatmap_figure_data.rds')
