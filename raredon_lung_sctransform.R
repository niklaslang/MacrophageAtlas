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

### split lung data ###
lung.list <- SplitObject(lung, split.by = "patient.ID")

### perform SCtransform ###
options(future.globals.maxSize = 8000 * 1024^2)
for (i in 1:length(lung.list)) {
  lung.list[[i]] <- SCTransform(lung.list[[i]], verbose = TRUE)
}

### select features for downstream integration ###
lung.features <- SelectIntegrationFeatures(object.list = lung.list, nfeatures = 3000)

### run PrepSCTIntegration ###
lung.list <- PrepSCTIntegration(object.list = lung.list, anchor.features = lung.features, verbose = TRUE)

### identify anchors ###
lung.anchors <- FindIntegrationAnchors(object.list = lung.list, anchor.features = lung.features,
                                       normalization.method = "SCT", k.filter = 60, verbose = TRUE)

### integrate datasets ###
lung.sctransform <- IntegrateData(anchorset = lung.anchors, normalization.method = "SCT", verbose = TRUE)

### PCA ### 
lung.sctransform <- RunPCA(lung.sctransform, verbose = TRUE)

# elbow plot
ElbowPlot(lung.sctransform, ndims = 50)

### UMAP visualisation ###
dims <- c(6,8,10,12,15)
for(d in dims){
  lung.sctransform <- RunUMAP(lung.sctransform, dims = 1:d)
  umap.batch.plot <- DimPlot(lung, reduction = "umap", label = F)
  eval(parse(text=paste0("UMAP.batch_dim", d, " <- umap.batch.plot")))
}

### save data ###
saveRDS(lung.sctransform, file = "/home/s1987963/MacrophageAtlas/raredon_lung_sctransform.rds")


