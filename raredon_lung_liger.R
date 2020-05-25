library(Seurat)
library(sctransform)
library(harmony)
library(liger)
library(reticulate)
library(umap)
library(dplyr)
library(ggplot2)
library(patchwork)

### load lung data ###
lung <- readRDS("/home/s1987963/MacrophageAtlas/raredon_lung.rds")


### normalisation ###
lung <- NormalizeData(lung)

### feature selection ###
lung <- FindVariableFeatures(lung)

### scale data ###
lung.liger <- ScaleData(lung, split.by = "patient.ID", do.center = FALSE)

### joint matrix factorization ###
## parameters ##
# k in the range 20 - 40 works well for most datasets
# lambda: regularization parameter, default is 5 but can lowered down to 1 in scenarios where the dataset differences are expected to be relatively small
# thresh: convergence threshold, default is 1e-6, lower values cause the algorithm to run longer
# max.iter: maximum number of iterations, default is 30
lung.liger <- RunOptimizeALS(lung.liger, k = 20, lambda = 5, split.by = "patient.ID")

### quantile normalisation ###
lung.liger <- RunQuantileAlignSNF(lung.liger, split.by = "patient.ID")

### UMAP visualisation ###
lung.liger <- RunUMAP(lung, dims = 1:ncol(lung.liger[["iNMF"]]), reduction = "iNMF")
liger.plots <- plotByDatasetAndCluster(lung.liger, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = T)
liger.plots[[1]] + liger.plots[[2]]

### save data ###
saveRDS(lung.liger, file = "/home/s1987963/MacrophageAtlas/raredon_lung_liger.rds")