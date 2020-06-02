library(Seurat)
library(reticulate)
library(umap)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

### load data ###
merge.path <- "/home/s1987963/ds_group/Niklas/raredon_lung/merge/"
logtransform.path <- "/home/s1987963/ds_group/Niklas/raredon_lung/merge/logtransform/"
logtransform.regress.nFeatures.path <- "/home/s1987963/ds_group/Niklas/raredon_lung/merge/logtransform.regress.nFeatures/"
lung <- readRDS("/home/s1987963/ds_group/Niklas/raredon_lung/raredon_lung.rds")

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
all.genes <- all.genes <- rownames(lung)
lung.logtransform <- ScaleData(lung)
lung.logtransform.regress.nFeatures <- ScaleData(lung, vars.to.regress = "nFeature_RNA")

### dimensionality reduction: PCA ###
lung.logtransform <- RunPCA(lung.logtransform, features = VariableFeatures(object = lung))
lung.logtransform.regress.nFeatures <- RunPCA(lung.logtransform.regress.nFeatures, features = VariableFeatures(object = lung))

## explore PCA results ##
# genes making up the first 15 PCs
print(lung.logtransform[["pca"]], dims = 1:15, nfeatures = 5)
print(lung.logtransform.regress.nFeatures[["pca"]], dims = 1:15, nfeatures = 5)

# plot genes making PC1 and PC2
VizDimLoadings(lung.logtransform, dims = 1:2, reduction = "pca")
VizDimLoadings(lung.logtransform.regress.nFeatures, dims = 1:2, reduction = "pca")

# plot PC1 against PC2
DimPlot(lung.logtransform, reduction = "pca")
DimPlot(lung.logtransform.regress.nFeatures, reduction = "pca")

# PC heatmap
for(i in 1:4){
  lower_limit <- 1+(6*(i-1))
  upper_limit <- 6*i
  DimHeatmap(lung.logtransform, dims = lower_limit:upper_limit, nfeatures = 20, cells = 500, balanced = TRUE)
  DimHeatmap(lung.logtransform.regress.nFeatures, dims = lower_limit:upper_limit, nfeatures = 20, cells = 500, balanced = TRUE)
}

## elbow plot ##
logtransform.elbow.plot <- ElbowPlot(lung.logtransform, ndims = 50)
png(paste0(merge.path,"logtransform.elbow.plot.png"), width=1000,height=600,units="px")
print(logtransform.elbow.plot)
dev.off()

logtransform.regress.nFeatures.elbow.plot <- ElbowPlot(lung.logtransform.regress.nFeatures, ndims = 50)
png(paste0(merge.path,"logtransform.regress.nFeatures.elbow.plot.png"), width=1000,height=600,units="px")
print(logtransform.regress.nFeatures.elbow.plot)
dev.off()

### LOGTRANSFORM.nREGRESS.NFEATURES ###
### clustering ###
# evaluate different numbers of PCs and resolutions
dims <- c(7,9,10,11,14)
res <- seq(0.5,1.5,0.1)

for(d in dims){
  lung.logtransform.regress.nFeatures <- RunUMAP(lung.logtransform.regress.nFeatures, dims=1:d, seed.use=1)
  for (r in res) {
    lung.logtransform.regress.nFeatures <- FindNeighbors(lung.logtransform.regress.nFeatures, dims = 1:d)
    lung.logtransform.regress.nFeatures <- FindClusters(lung.logtransform.regress.nFeatures, resolution = r)
    umap.plot <- DimPlot(lung.logtransform.regress.nFeatures, reduction = "umap", label = F, pt.size = 0.1)
    batch.plot <- DimPlot(lung.logtransform.regress.nFeatures, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
    eval.plot <- umap.plot + batch.plot
    png(paste0(logtransform.regress.nFeatures.path, "UMAP_dim", d, "_res", r, ".png"), width=1400, height=600, units="px")
    print(eval.plot)
    dev.off()
  }
}

## best (preliminary) clustering ##
# first 10 PCs, resolution 0.9
lung.logtransform.regress.nFeatures <- FindNeighbors(lung.logtransform.regress.nFeatures, dims = 1:10)
lung.logtransform.regress.nFeatures <- FindClusters(lung.logtransform.regress.nFeatures, resolution = 0.9)
# add UMAP
lung.logtransform.regress.nFeatures <- RunUMAP(lung.logtransform.regress.nFeatures, dims=1:10, seed.use=1)
# add tSNE
lung.logtransform.regress.nFeatures <- RunTSNE(lung.logtransform.regress.nFeatures, dims=1:10, seed.use=1)

## explore clustering at patient level ##
patient.clustering <- DimPlot(lung.logtransform.regress.nFeatures, group.by = "seurat_clusters", split.by = "patient.ID", ncol = 7)
png(paste0(logtransform.regress.nFeatures.path,"patient.clustering.png"), width=1800,height=600,units="px")
print(patient.clustering)
dev.off()

### save data ###
saveRDS(lung.logtransform.regress.nFeatures, file = "/home/s1987963/ds_group/Niklas/raredon_lung/merge/raredon_lung_logtransform.regress.nFeatures.rds")

## QC metrics at cluster level ##
cluster.qc.heatmap <- FeaturePlot(lung.logtransform.regress.nFeatures, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(logtransform.regress.nFeatures.path,"cluster.qc.heatmap.png"), width=1500,height=500,units="px")
print(cluster.qc.heatmap)
dev.off()

## marker gene visualization ##
# macrophage markers
macrophage.genes <- c("CSF1R", "LYZ", "HLA-DRA", "ITGAX", "ITGAM", "C1QB","MRC1", "MARCO") #MSR1
# monocyte markers
monocyte.genes <- c("CD14", "MNDA", "S100A8","S100A9")
# dendritic cell markers
dc.genes <- c("CD1C","XCR1", "CD86")
# lineage markers
lineage.genes <- c("EPCAM", "CD3D", "CDH5", "PECAM1", "PDGFRB", "PDGFRA", "CD34", "GZMA", "CD79A", "CD79B")

# set default essay to RNA counts
DefaultAssay(lung.logtransform.regress.nFeatures) <- "RNA"
# Normalize RNA data for visualization purposes
lung.logtransform.regress.nFeatures <- NormalizeData(lung.logtransform.regress.nFeatures, verbose = TRUE)

# feature plot with macrophage markers
macrophage.markers <- FeaturePlot(lung.logtransform.regress.nFeatures, features = macrophage.genes, pt.size = 0.2, ncol = 4) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(logtransform.regress.nFeatures.path,"macrophage.markers.png"), width=1600,height=800,units="px")
print(macrophage.markers)
dev.off()

# feature plot with monocyte markers
monocyte.markers <- FeaturePlot(lung.logtransform.regress.nFeatures, features = monocyte.genes, pt.size = 0.2, ncol = 2) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(logtransform.regress.nFeatures.path,"monocyte.markers.png"), width=800,height=800,units="px")
print(monocyte.markers)
dev.off()

# feature plot with DC markers
dc.markers <- FeaturePlot(lung.logtransform.regress.nFeatures, features = dc.genes, pt.size = 0.2, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(logtransform.regress.nFeatures.path,"dc.markers.png"), width=1200,height=400,units="px")
print(dc.markers)
dev.off()

# feature plot with more general lineage markers
lineage.markers <- FeaturePlot(lung.logtransform.regress.nFeatures, features = lineage.genes, pt.size = 0.2, ncol =3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(logtransform.regress.nFeatures.path,"lineage.markers.png"), width=1200,height=1200,units="px")
print(lineage.markers)
dev.off()

### find cluster marker genes ###
# find markers for every cluster compared to all remaining cells, report only the positive ones
lung.markers <- FindAllMarkers(lung.logtransform.regress.nFeatures, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
lung.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

### LOGTRANSFORM ###
### clustering ###
# evaluate different numbers of PCs and resolutions
dims <- c(7,10,14,22,30)
res <- seq(0.1,1.5,0.1)

for(d in dims){
  lung.logtransform <- RunUMAP(lung.logtransform, dims=1:d, seed.use=1)
  for (r in res) {
    lung.logtransform <- FindNeighbors(lung.logtransform, dims = 1:d)
    lung.logtransform <- FindClusters(lung.logtransform, resolution = r)
    umap.plot <- DimPlot(lung.logtransform, reduction = "umap", label = F, pt.size = 0.1)
    batch.plot <- DimPlot(lung.logtransform, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
    eval.plot <- umap.plot + batch.plot
    png(paste0(logtransform.path, "UMAP_dim", d, "_res", r, ".png"), width=1400, height=600, units="px")
    print(eval.plot)
    dev.off()
  }
}

## best (preliminary) clustering ##
# first 7 PCs, resolution 1.0
lung.logtransform <- FindNeighbors(lung.logtransform, dims = 1:7)
lung.logtransform <- FindClusters(lung.logtransform, resolution = 1.0)
# add UMAP
lung.logtransform <- RunUMAP(lung.logtransform, dims=1:7, seed.use=1)
# add tSNE
lung.logtransform <- RunTSNE(lung.logtransform, dims=1:7, seed.use=1)

## explore clustering at patient level ##
patient.clustering <- DimPlot(lung.logtransform, group.by = "seurat_clusters", split.by = "patient.ID", ncol = 7)
png(paste0(logtransform.path,"patient.clustering.png"), width=1800,height=600,units="px")
print(patient.clustering)
dev.off()

## triple UMAP ##
umap.1 <- DimPlot(lung.logtransform, reduction = "umap", label = F, pt.size = 0.1)
umap.2 <- DimPlot(lung.logtransform, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
umap.3 <- DimPlot(lung.logtransform, reduction = "umap", group.by = "scrublet_auto", pt.size = 0.1)
triple.umap <- umap.1 + umap.2 + umap.3
png(paste0(logtransform.path,"triple.umap.png"), width=1800,height=600,units="px")
print(triple.umap)
dev.off()

### save data ###
saveRDS(lung.logtransform, file = "/home/s1987963/ds_group/Niklas/raredon_lung/merge/raredon_lung_logtransform.rds")

## QC metrics at cluster level ##
cluster.qc.heatmap <- FeaturePlot(lung.logtransform, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(logtransform.path,"cluster.qc.heatmap.png"), width=1500,height=500,units="px")
print(cluster.qc.heatmap)
dev.off()

## marker gene visualization ##
# macrophage markers
macrophage.genes <- c("CSF1R", "LYZ", "HLA-DRA", "ITGAX", "ITGAM", "C1QB","MRC1", "MARCO") #MSR1
# monocyte markers
monocyte.genes <- c("CD14", "MNDA", "S100A8","S100A9")
# dendritic cell markers
dc.genes <- c("CD1C","XCR1", "CD86")
# lineage markers
lineage.genes <- c("EPCAM", #epithelial cells
                   "CDH5", "PECAM1", #endothelial cells
                   "CD3D", "GZMA", #T-cells,
                   "PDGFRB", "PDGFRA", #mesenchymal cells
                   "CD34",  
                   "CD79A", "CD79B" #B-cells
                   )

# set default essay to RNA counts
DefaultAssay(lung.logtransform) <- "RNA"
# Normalize RNA data for visualization purposes
lung.logtransform <- NormalizeData(lung.logtransform, verbose = TRUE)

# feature plot with macrophage markers
macrophage.markers <- FeaturePlot(lung.logtransform, features = macrophage.genes, pt.size = 0.2, ncol = 4) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(logtransform.path,"macrophage.markers.png"), width=1600,height=800,units="px")
print(macrophage.markers)
dev.off()

# feature plot with monocyte markers
monocyte.markers <- FeaturePlot(lung.logtransform, features = monocyte.genes, pt.size = 0.2, ncol = 2) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(logtransform.path,"monocyte.markers.png"), width=800,height=800,units="px")
print(monocyte.markers)
dev.off()

# feature plot with DC markers
dc.markers <- FeaturePlot(lung.logtransform, features = dc.genes, pt.size = 0.2, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(logtransform.path,"dc.markers.png"), width=1200,height=400,units="px")
print(dc.markers)
dev.off()

# feature plot with more general lineage markers
lineage.markers <- FeaturePlot(lung.logtransform, features = lineage.genes, pt.size = 0.2, ncol =3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(logtransform.path,"lineage.markers.png"), width=1200,height=1200,units="px")
print(lineage.markers)
dev.off()

### find cluster marker genes ###
# find markers for every cluster compared to all remaining cells, report only the positive ones
lung.markers <- FindAllMarkers(lung.logtransform, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
lung.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

