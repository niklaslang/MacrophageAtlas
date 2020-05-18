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
lung <- ScaleData(lung, vars.to.regress = "percent.mt")

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
UMAPPlot(lung, reduction = "umap")
UMAPPlot(lung, reduction = "umap", group.by = "orig.ident")

## tSNE ##
lung <- RunTSNE(lung, dims = 1:15)
TSNEPlot(lung, reduction = "tsne")
TSNEPlot(lung, reduction = "tsne", group.by = "orig.ident")

### save data ###
saveRDS(lung, file = "/home/s1987963/MacrophageAtlas/raredon_lung_all.rds")

### read data ###
lung <- readRDS("/home/s1987963/MacrophageAtlas/raredon_lung_all.rds")

### find cluster marker genes ###
# find markers for every cluster compared to all remaining cells, report only the positive ones
lung.markers <- FindAllMarkers(lung, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
lung.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

# plot macrophage marker genes
MNP.markers <- list("CSF1R", "LYZ",  "CD68","HLA-DRA", "ITGAX", "ITGAM", "C1QB", "MRC1", # macrophage marker
                    "CD14", "MNDA", # monocyte + macrophage marker
                    "S100A8","S100A9", # monocyte marker
                    "LILRA4","CD1C","XCR1", "CD86" # DC marker
                    )

### save data ###
saveRDS(lung, file = "/home/s1987963/MacrophageAtlas/raredon_lung_all.rds")

# dot plot
DotPlot(lung, features = MNP.markers, split.by = "seurat_clusters") + RotatedAxis()
marker.vln.1 <- VlnPlot(lung, features = macrophage.markers[1])
marker.vln.2 <- VlnPlot(lung, features = macrophage.markers[2])
marker.vln.3 <- VlnPlot(lung, features = macrophage.markers[3])
marker.vln.4 <- VlnPlot(lung, features = macrophage.markers[4])
marker.vln.5 <- VlnPlot(lung, features = macrophage.markers[5])
marker.vln.6 <- VlnPlot(lung, features = macrophage.markers[6])
marker.vln.7 <- VlnPlot(lung, features = macrophage.markers[7])
marker.vln.8 <- VlnPlot(lung, features = macrophage.markers[8])
marker.vln.9 <- VlnPlot(lung, features = macrophage.markers[9])
marker.vln.10 <- VlnPlot(lung, features = macrophage.markers[10])
marker.vln.11 <- VlnPlot(lung, features = macrophage.markers[11])
marker.vln.12 <- VlnPlot(lung, features = macrophage.markers[12])
marker.vln.13 <- VlnPlot(lung, features = macrophage.markers[13])
marker.vln.14 <- VlnPlot(lung, features = macrophage.markers[14])
marker.vln.15 <- VlnPlot(lung, features = macrophage.markers[15])
marker.vln.16 <- VlnPlot(lung, features = macrophage.markers[16])

# plot violins
marker.vln.1 + marker.vln.2
marker.vln.3 + marker.vln.4

# feature plot
FeaturePlot(lung, features = macrophage.markers[c(-3,-6)])

# heat map
top10 <- lung.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(lung, features = top10$gene) + NoLegend()

### save file ###
saveRDS(lung, file = "/home/s1987963/MacrophageAtlas/raredon_lung_final.rds")

