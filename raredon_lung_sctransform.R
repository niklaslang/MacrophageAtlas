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
lung <- readRDS("/home/s1987963/MacrophageAtlas/raredon_lung_4.rds")

### split lung data ###
lung.list <- SplitObject(lung, split.by = "patient.ID")

### perform SCtransform ###
options(future.globals.maxSize = 4000 * 1024^2)
for (i in 1:length(lung.list)) {
  lung.list[[i]] <- SCTransform(lung.list[[i]], verbose = FALSE)
}

### select features for downstream integration ###
lung.features <- SelectIntegrationFeatures(object.list = lung.list, nfeatures = 3000)

### run PrepSCTIntegration ###
lung.list <- PrepSCTIntegration(object.list = lung.list, anchor.features = lung.features, 
                                    verbose = FALSE)

### identify anchors ###
lung.anchors <- FindIntegrationAnchors(object.list = lung.list, normalization.method = "SCT", 
                                           anchor.features = lung.features, verbose = FALSE)

### integrate datasets ###
lung.sctransform <- IntegrateData(anchorset = lung.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)

### PCA ### 
lung.sctransform <- RunPCA(lung.sctransform, verbose = FALSE)
# elbow plot
ElbowPlot(lung.sctransform, ndims = 50)

### UMAP visualisation ###
# first 10 PCs
lung.sctransform <- RunUMAP(lung.sctransform, dims = 1:10)
DimPlot(lung.sctransform, group.by = "patient.ID", pt.size = 0.1)
# first 15 PCs
lung.sctransform <- RunUMAP(lung.sctransform, dims = 1:15)
DimPlot(lung.sctransform, group.by = "patient.ID", pt.size = 0.1)




