library(Seurat)
library(reticulate)
library(umap)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

### path variables ###
lung.path <- "/home/s1987963/ds_group/Niklas/reyfman_lung/reyfman_lung_healthy.rds"
# logtransform.path <- "/home/s1987963/ds_group/Niklas/reyfman_lung/logtransform/uncorrected/" # uncorrected/unadjusted
logtransform.path <- "/home/s1987963/ds_group/Niklas/reyfman_lung/logtransform/regress.nFeatures/" # corrected for nFeature
# logtransform.path <- "/home/s1987963/ds_group/Niklas/reyfman_lung/logtransform/regress.nCounts.mito/" # corrected for nCounts and mitochondrial fraction
  
### read healthy lung data ###
lung <- readRDS(lung.path)

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
# lung.logtransform <- ScaleData(lung) # uncorrected
lung.logtransform <- ScaleData(lung, vars.to.regress = "nFeature_RNA") # adjusted for nFeatures
# lung.logtransform <- ScaleData(lung, vars.to.regress = c("nCount_RNA", "percent.mt")) # adjusted for nCounts + percent.mt

### dimensionality reduction: PCA ###
lung.logtransform <- RunPCA(lung.logtransform, features = VariableFeatures(object = lung))

## elbow plot ##
logtransform.elbow.plot <- ElbowPlot(lung.logtransform, ndims = 50)
png(paste0(logtransform.path,"logtransform.elbow.plot.png"), width=1000,height=600,units="px")
print(logtransform.elbow.plot)
dev.off()

### visualize batch effect - pre-selection of dimensions -  ##
dims <- c(4,6,9,15,17)
for(d in dims){
  lung.logtransform <- RunUMAP(lung.logtransform, dims = 1:d, seed.use=1)
  umap.batch.plot <- DimPlot(lung.logtransform, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
  png(paste0(logtransform.path, "UMAP_dim", d, "batch.png"), width=1400, height=600, units="px")
  print(umap.batch.plot)
  dev.off()
}

### clustering ###
## evaluate different numbers of PCs and resolutions ##
### clustering ###
# evaluate different numbers of PCs and resolutions
dims <- c(6,9) #adjusted for nFeatures
# dims <- c(4,6,9,11,15,17,25) #adjusted for nCounts + percent.mito
#dims <- c(6,8,9,11,15,20,30) # uncorrected/adjusted
res <- seq(0.3,1.5,0.1)
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

### definite clustering ###
## best (preliminary) clustering ##
# uncorrected/adjusted: 8 first PCs, resolution 0.3
# adjusted for nFeatures: 6 first PCs, resolution 0.4
# adjusted for nCounts + percent.mito: 6 first PCs, resolution 0.4
lung.logtransform <- FindNeighbors(lung.logtransform, dims = 1:6)
lung.logtransform <- FindClusters(lung.logtransform, resolution = 0.4)
# add UMAP
lung.logtransform <- RunUMAP(lung.logtransform, dims=1:6, seed.use=1)
# add tSNE
lung.logtransform <- RunTSNE(lung.logtransform, dims=1:6, seed.use=1)

### save data ###
saveRDS(lung.logtransform, file = paste0(logtransform.path, "reyfman_lung_logtransform.rds"))

## explore clustering at patient level ##
patient.clustering <- DimPlot(lung.logtransform, group.by = "seurat_clusters", split.by = "patient.ID", ncol = 4)
png(paste0(logtransform.path,"patient.clustering.png"), width=1600,height=800,units="px")
print(patient.clustering)
dev.off()

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
dc.genes <- c("CD1C","XCR1", "CD86", "CCL17", "S100B", "RGS1")
# lineage markers
lineage.genes <- c("EPCAM", #epithelial cells
                   "CDH5", "PECAM1", "VWF", "KDR", #endothelial cells
                   "PDGFRA", "PDGFRB", "ACTA2", "MYH11",  #mesenchymal cells - maybe "CD34"?
                   "PTPRC", #immune cells
                   "TPSB2", #mast cells - maybe "IL1RL1"?
                   "CD3D", "GZMA", #T-cells
                   "CD79A", "CD79B" #B-cells
)

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
png(paste0(logtransform.path,"dc.markers.png"), width=1200,height=800,units="px")
print(dc.markers)
dev.off()

# feature plot with more general lineage markers
lineage.markers <- FeaturePlot(lung.logtransform, features = lineage.genes, pt.size = 0.2, ncol = 5) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(logtransform.path,"lineage.markers.png"), width=1800,height=1200,units="px")
print(lineage.markers)
dev.off()
