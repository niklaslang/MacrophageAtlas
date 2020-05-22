library(Seurat)
library(umap)
library(reticulate)
library(dplyr)
library(ggplot2)
library(patchwork)

### load lung.all data ###
lung <- readRDS("/home/s1987963/MacrophageAtlas/raredon_lung_all.rds")

### normalisation ###
lung <- NormalizeData(lung, normalization.method = "LogNormalize", scale.factor = 10000)

### feature selection ###
# selecting the most highly variable genes
lung <- FindVariableFeatures(lung, selection.method = "vst", nfeatures = 2000)

# 10 most highly variable genes
top10 <- head(VariableFeatures(lung), 10)

# plot variable features with and without labels
HVGs.plot <- VariableFeaturePlot(lung)
HVGs.plot <- LabelPoints(plot = HVGs.plot, points = top10, repel = TRUE)
HVGs.plot

### scale data ###
lung <- ScaleData(lung, vars.to.regress = "nFeature_RNA")

### dimensionality reduction: PCA ###
lung <- RunPCA(lung, features = VariableFeatures(object = lung))

## explore PCA results ##
# genes making up the first 15 PCs
print(lung[["pca"]], dims = 1:15, nfeatures = 5)

# plot genes making PC1 and PC2
VizDimLoadings(lung, dims = 1:2, reduction = "pca")

# plot PC1 against PC2
DimPlot(lung, reduction = "pca")

# PC heatmap
for(i in 1:4){
  lower_limit <- 1+(6*(i-1))
  upper_limit <- 6*i
  DimHeatmap(lung, dims = lower_limit:upper_limit, nfeatures = 20, cells = 500, balanced = TRUE)
}

## elbow plot ##
ElbowPlot(lung, ndims = 50)

### save data ###
saveRDS(lung, file = "/home/s1987963/MacrophageAtlas/raredon_lung_pca.rds")

### clustering ###
# evaluate different numbers of PCs and resolutions
dims <- c(7,9,10,11,13,14)
res <- seq(0.2,1.2,0.1)

for(d in dims){
  lung <- RunUMAP(lung, dims=1:d, seed.use=1)
  for (r in res) {
    lung <- FindNeighbors(lung, dims = 1:d)
    lung <- FindClusters(lung, resolution = r)
    umap.plot <- DimPlot(lung, reduction = "umap", label = F)
    batch.plot <- DimPlot(lung, reduction = "umap", group.by = "orig.ident")
    eval(parse(text=paste0("UMAP_dim", d, "_res", r, " <- umap.plot + batch.plot")))
  }
}

### save data ###
saveRDS(lung, file = "/home/s1987963/MacrophageAtlas/raredon_lung_all.rds")

### read data ###
lung <- readRDS("/home/s1987963/MacrophageAtlas/raredon_lung_all.rds")

### find cluster marker genes ###
# find markers for every cluster compared to all remaining cells, report only the positive ones
lung.markers <- FindAllMarkers(lung, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
lung.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

## marker genes ##
# all MNP markers
MNP.genes <- c("CSF1R", "LYZ", "HLA-DRA", "ITGAX", "ITGAM", "C1QB", "MRC1","CD14", "MNDA",
               "S100A8","S100A9", "CD1C","XCR1", "CD86" )
# macrophage markers
macrophage.genes <- c("CSF1R", "LYZ", "HLA-DRA", "ITGAX", "ITGAM", "C1QB", "MRC1")
# monocyte markers
monocyte.genes <- c("CD14", "MNDA", "S100A8","S100A9")
# dendritic cell markers
dc.genes <- c("CD1C","XCR1", "CD86")

# dot plot
DotPlot(lung, features = MNP.genes) + RotatedAxis()

# feature plot
FeaturePlot(lung, features = MNP.genes)

### save file ###
saveRDS(lung, file = "/home/s1987963/MacrophageAtlas/raredon_lung_final.rds")

#### processing lung.11 data ###
### load data ###
lung.11 <- readRDS("/home/s1987963/MacrophageAtlas/raredon_lung_11.rds")

### normalisation ###
lung.11 <- NormalizeData(lung.11, normalization.method = "LogNormalize", scale.factor = 10000)

### feature selection ###
# selecting the most highly variable genes
lung.11 <- FindVariableFeatures(lung.11, selection.method = "vst", nfeatures = 2000)

### scale data ###
lung.11 <- ScaleData(lung.11, vars.to.regress = "nFeature_RNA")

### dimensionality reduction: PCA ###
lung.11 <- RunPCA(lung.11, features = VariableFeatures(object = lung.11))

## explore PCA results ##
# genes making up the first 15 PCs
print(lung.11[["pca"]], dims = 1:15, nfeatures = 5)

# plot genes making PC1 and PC2
VizDimLoadings(lung.11, dims = 1:2, reduction = "pca")

# plot PC1 against PC2
DimPlot(lung.11, reduction = "pca")

# PC heatmap
for(i in 1:4){
  lower_limit <- 1+(6*(i-1))
  upper_limit <- 6*i
  DimHeatmap(lung.11, dims = lower_limit:upper_limit, nfeatures = 20, cells = 500, balanced = TRUE)
}

## elbow plot ##
ElbowPlot(lung.11, ndims = 50)

### save data ###
saveRDS(lung.11, file = "/home/s1987963/MacrophageAtlas/raredon_lung_11_pca.rds")

### clustering ###
## evaluate different numbers of PCs and resolutions ##
# create batch plots
dims <- c(5,6,9,10,11,16,19)
for(d in dims){
  lung.11 <- RunUMAP(lung.11, dims=1:d, seed.use=1)
  batch.plot <- DimPlot(lung.11, reduction = "umap", group.by = "orig.ident", pt.size = 0.1)
  eval(parse(text=paste0("UMAP_dim", d, "batch <- batch.plot")))
}

# compare clustering and batch effect
dims <- c(5,6)
res <- c(.3,1.4,1.5,1.6,1.7,1.8)
#res <- seq(0.2,2.0,0.2)
for(d in dims){
  lung.11 <- RunUMAP(lung.11, dims=1:d, seed.use=1)
  for (r in res) {
    lung.11 <- FindNeighbors(lung.11, dims = 1:d)
    lung.11 <- FindClusters(lung.11, resolution = r)
    umap.plot <- DimPlot(lung.11, reduction = "umap", label = F)
    batch.plot <- DimPlot(lung.11, reduction = "umap", group.by = "orig.ident")
    eval(parse(text=paste0("UMAP_dim", d, "_res", r, " <- umap.plot + batch.plot")))
  }
}

#### processing lung.9 data ###
### load data ###
lung.9 <- readRDS("/home/s1987963/MacrophageAtlas/raredon_lung_9.rds")

### normalisation ###
lung.9 <- NormalizeData(lung.9, normalization.method = "LogNormalize", scale.factor = 10000)

### feature selection ###
# selecting the most highly variable genes
lung.9 <- FindVariableFeatures(lung.9, selection.method = "vst", nfeatures = 2000)

### scale data ###
lung.9 <- ScaleData(lung.9)

### dimensionality reduction: PCA ###
lung.9 <- RunPCA(lung.9, features = VariableFeatures(object = lung.9))

## explore PCA results ##
# genes making up the first 15 PCs
print(lung.9[["pca"]], dims = 1:15, nfeatures = 5)

# plot genes making PC1 and PC2
VizDimLoadings(lung.9, dims = 1:2, reduction = "pca")

# plot PC1 against PC2
DimPlot(lung.9, reduction = "pca")

# PC heatmap
for(i in 1:4){
  lower_limit <- 1+(6*(i-1))
  upper_limit <- 6*i
  DimHeatmap(lung.9, dims = lower_limit:upper_limit, nfeatures = 20, cells = 500, balanced = TRUE)
}

## elbow plot ##
ElbowPlot(lung.9, ndims = 50)

### save data ###
saveRDS(lung.9, file = "/home/s1987963/MacrophageAtlas/raredon_lung_9_pca.rds")

### clustering ###
# evaluate different numbers of PCs and resolutions
dims <- c(6,9,10,11,15,20)
res <- seq(0.2,1.2,0.1)

#

# 
for(d in dims){
  lung.9 <- RunUMAP(lung.9, dims=1:d, seed.use=1)
  for (r in res) {
    lung.9 <- FindNeighbors(lung.9, dims = 1:d)
    lung.9 <- FindClusters(lung.9, resolution = r)
    umap.plot <- DimPlot(lung.9, reduction = "umap", label = F)
    batch.plot <- DimPlot(lung.9, reduction = "umap", group.by = "orig.ident")
    eval(parse(text=paste0("UMAP_dim", d, "_res", r, " <- umap.plot + batch.plot")))
  }
}
