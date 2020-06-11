library(Seurat)
library(reticulate)
library(umap)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

### path.variables ###
liver.path <- "/home/s1987963/ds_group/Niklas/macparland_liver/macparland_liver_filtered.rds"
logtransform.path <- "/home/s1987963/ds_group/Niklas/macparland_liver/logtransform/uncorrected/"
# logtransform.path <- "/home/s1987963/ds_group/Niklas/macparland_liver/logtransform/regress.nFeatures/" # corrected for nFeatures
# logtransform.path <- "/home/s1987963/ds_group/Niklas/macparland_liver/logtransform/regress.nCounts.mito/" # corrected for nCounts and mitochondrial fraction

### load data ###
liver <- readRDS(liver.path)

### normalisation ###
liver.logtransform <- NormalizeData(liver, normalization.method = "LogNormalize", scale.factor = 10000)

### feature selection ###
liver.logtransform <- FindVariableFeatures(liver.logtransform, selection.method = "vst", nfeatures = 2000)

# 10 most highly variable genes
top10 <- head(VariableFeatures(liver.logtransform), 10)

# plot variable features with and without labels
HVGs.plot <- VariableFeaturePlot(liver.logtransform)
HVGs.plot <- LabelPoints(plot = HVGs.plot, points = top10, repel = TRUE)
HVGs.plot

### scale data ###
liver.logtransform <- ScaleData(liver.logtransform)
# liver.logtransform <- ScaleData(liver.logtransform, vars.to.regress = "nFeature_RNA") # adjusted for nFeatures
# liver.logtransform <- ScaleData(liver.logtransform, vars.to.regress = c("nCount_RNA", "percent.mt")) # adjusted for nCounts + percent.mt

### dimensionality reduction: PCA ###
liver.logtransform <- RunPCA(liver.logtransform, features = VariableFeatures(object = liver.logtransform))

## elbow plot ##
logtransform.elbow.plot <- ElbowPlot(liver.logtransform, ndims = 50)
png(paste0(logtransform.path,"logtransform.elbow.plot.png"), width=1000,height=600,units="px")
print(logtransform.elbow.plot)
dev.off()

### visualize batch effect - pre-selection of dimensions ##
dims <- c(8,12,18,25,30)
for(d in dims){
  liver.logtransform <- RunUMAP(liver.logtransform, dims = 1:d, seed.use=1)
  umap.batch.plot <- DimPlot(liver.logtransform, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
  png(paste0(logtransform.path, "UMAP_dim", d, ".batch.png"), width=1400, height=800, units="px")
  print(umap.batch.plot)
  dev.off()
}

### clustering ###
## evaluate different numbers of PCs and resolutions ##
### clustering ###
# evaluate different numbers of PCs and resolutions
dims <- c(8,12,18) # uncorrected/adjusted
# dims <- c(4,6,9,11,15,17,25) #adjusted for nCounts + percent.mito
#dims <- c(6,8,9,11,15,20,30) # adjusted for nFeatures
res <- seq(0.3,1.5,0.1)
for(d in dims){
  liver.logtransform <- RunUMAP(liver.logtransform, dims=1:d, seed.use=1)
  batch.plot <- DimPlot(liver.logtransform, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
  for (r in res) {
    liver.logtransform <- FindNeighbors(liver.logtransform, dims = 1:d)
    liver.logtransform <- FindClusters(liver.logtransform, resolution = r)
    umap.plot <- DimPlot(liver.logtransform, reduction = "umap", label = F, pt.size = 0.1)
    eval.plot <- umap.plot + batch.plot
    png(paste0(logtransform.path, "UMAP_dim", d, "_res", r, ".png"), width=1400, height=600, units="px")
    print(eval.plot)
    dev.off()
  }
}

### definite clustering ###
## best (preliminary) clustering ##
# uncorrected/adjusted: 8 first PCs, resolution XY
# adjusted for nFeatures: XY first PCs, resolution XY
# adjusted for nCounts + percent.mito: XY first PCs, resolution XY
liver.logtransform <- FindNeighbors(liver.logtransform, dims = 1:8)
liver.logtransform <- FindClusters(liver.logtransform, resolution = 1.1)
# add UMAP
liver.logtransform <- RunUMAP(liver.logtransform, dims=1:8, seed.use=1)
# add tSNE
liver.logtransform <- RunTSNE(liver.logtransform, dims=1:8, seed.use=1)

### save data ###
saveRDS(liver.logtransform, paste0(logtransform.path, "macparland_liver_logtransform.rds"))

## explore clustering at patient level ##
patient.clustering <- DimPlot(liver.logtransform, group.by = "seurat_clusters", split.by = "patient.ID", ncol = 9)
png(paste0(logtransform.path,"patient.clustering.png"), width=1600,height=800,units="px")
print(patient.clustering)
dev.off()

## QC metrics at cluster level ##
cluster.qc.heatmap <- FeaturePlot(liver.logtransform, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(logtransform.path,"cluster.qc.heatmap.png"), width=1500,height=500,units="px")
print(cluster.qc.heatmap)
dev.off()

## marker gene visualization ##
# macrophage markers
macrophage.genes <- c("CSF1R", "CD68", "LYZ", "HLA-DRA", "ITGAX", "ITGAM", "C1QB","MRC1", "MARCO", "MSR1")
# monocyte markers
monocyte.genes <- c("CD14", "MNDA", "S100A8","S100A9")
# dendritic cell markers
dc.genes <- c("CD1C","XCR1", "CD86", "CCL17", "S100B", "RGS1")
# liver markers
liver.genes <- c("ALB", "AFP", # hepatocytes
                 "CALCRL","CD32B", # LSECs
                 "KRT19","EPCAM", #cholangiocytes
                 "ACTA2","COL1A1" # Hepatic Stellate Cells
                 )
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
macrophage.markers <- FeaturePlot(liver.logtransform, features = macrophage.genes, pt.size = 0.2, ncol = 5) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(logtransform.path,"macrophage.markers.png"), width=1600,height=800,units="px")
print(macrophage.markers)
dev.off()

# feature plot with monocyte markers
monocyte.markers <- FeaturePlot(liver.logtransform, features = monocyte.genes, pt.size = 0.2, ncol = 2) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(logtransform.path,"monocyte.markers.png"), width=800,height=800,units="px")
print(monocyte.markers)
dev.off()

# feature plot with DC markers
dc.markers <- FeaturePlot(liver.logtransform, features = dc.genes, pt.size = 0.2, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(logtransform.path,"dc.markers.png"), width=1200,height=800,units="px")
print(dc.markers)
dev.off()

# feature plot with liver markers
liver.markers <- FeaturePlot(liver.logtransform, features = liver.genes, pt.size = 0.2, ncol = 4) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(logtransform.path,"liver.markers.png"), width=1600,height=800,units="px")
print(liver.markers)
dev.off()

# feature plot with lineage markers
lineage.markers <- FeaturePlot(liver.logtransform, features = lineage.genes, pt.size = 0.2, ncol = 5) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(logtransform.path,"lineage.markers.png"), width=1500,height=500,units="px")
print(lineage.markers)
dev.off()