library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(umap)
library(reticulate)

### load lung data ###
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
lung <- ScaleData(lung, features = rownames(lung), vars.to.regress = "percent.mt")

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

### clustering ###
# dims = 15, vary resolution
lung <- FindNeighbors(lung, dims = 1:15)
lung <- FindClusters(lung, resolution = 0.5)

### visualisation ###
## UMAP ##
lung <- RunUMAP(lung, dims = 1:15)
DimPlot(lung, reduction = "umap")
DimPlot(lung, reduction = "umap", group.by = "orig.ident")

## tSNE ##
lung <- RunTSNE(lung, dims = 1:15)
DimPlot(lung, reduction = "tsne")
DimPlot(lung, reduction = "tsne", group.by = "orig.ident")

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

