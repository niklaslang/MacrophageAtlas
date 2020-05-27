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

### UMAP visualisation of batch effect ###
dims <- c(5,7,9,11,14,17)
for(d in dims){
  lung.sctransform <- RunUMAP(lung.sctransform, dims = 1:d, seed.use=1)
  umap.batch.plot <- DimPlot(lung.sctransform, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
  eval(parse(text=paste0("UMAP.batch_dim", d, " <- umap.batch.plot")))
}

### clustering ###
# evaluate different numbers of PCs and resolutions
dims <- c(7,9,11,14,17)
res <- seq(0.2,1.6,0.1)
for(d in dims){
  lung.sctransform <- RunUMAP(lung.sctransform, dims=1:d, seed.use=1)
  for (r in res) {
    lung.sctransform <- FindNeighbors(lung.sctransform, dims = 1:d)
    lung.sctransform <- FindClusters(lung.sctransform, resolution = r)
    umap.plot <- DimPlot(lung.sctransform, reduction = "umap", pt.size = 0.1)
    batch.plot <- DimPlot(lung.sctransform, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
    eval(parse(text=paste0("UMAP_dim", d, "_res", r, " <- umap.plot + batch.plot")))
  }
}

### preliminary clustering ###
## first 11 PCs ##
## resolution 1.0 ##
lung.sctransform <- RunUMAP(lung.sctransform, dims=1:11, seed.use=1)
lung.sctransform <- FindNeighbors(lung.sctransform, dims = 1:11)
lung.sctransform <- FindClusters(lung.sctransform, resolution = 1.0)

### explore clustering at patient level ###
lung.sctransform$cluster <- Idents(lung.sctransform)
DimPlot(lung.sctransform, group.by = "cluster", split.by = "patient.ID", ncol = 5)

### save data ###
saveRDS(lung.sctransform, file = "/home/s1987963/ds_group/Niklas/raredon_lung/raredon_lung_sctransform.rds")

### marker gene visualization ###
# all MNP markers
MNP.genes <- c("CSF1R", "LYZ", "HLA-DRA", "ITGAX", "ITGAM", "C1QB", "MRC1","CD14", "MNDA",
               "S100A8","S100A9", "CD1C","XCR1", "CD86" )

# macrophage markers
macrophage.genes <- c("CSF1R", "LYZ", "HLA-DRA", "ITGAX", "ITGAM", "C1QB","MRC1","CD68")
# monocyte markers
monocyte.genes <- c("CD14", "MNDA", "S100A8","S100A9")
# dendritic cell markers
dc.genes <- c("CD1C","XCR1", "CD86")
# lineage markers
lineage.genes <- c("EPCAM", "CD3D", "CDH5", "PECAM1", "PDGFRB", "PDGFRA", "CD34", "GZMA", "CD79A", "CD79B")

# set default essay to RNA counts
DefaultAssay(lung.sctransform) <- "RNA"
# Normalize RNA data for visualization purposes
lung.sctransform <- NormalizeData(lung.sctransform, verbose = TRUE)
# feature plot with macrophage markers
FeaturePlot(lung.sctransform, features = macrophage.genes, pt.size = 0.2)
# feature plot with monocyte markers
FeaturePlot(lung.sctransform, features = monocyte.genes, pt.size = 0.2)
# feature plot with DC markers
FeaturePlot(lung.sctransform, features = dc.genes, pt.size = 0.5, ncol = 3)
# feature plot with more general lineage markers
FeaturePlot(lung.sctransform, features = lineage.genes, pt.size = 0.5, ncol =3)








