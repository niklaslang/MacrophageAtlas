library(Seurat)
library(sctransform)
library(reticulate)
library(umap)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

### path variables ###
liver.path <- "/home/s1987963/ds_group/Niklas/macparland_liver/macparland_liver_filtered.rds"
# sctransform.path <- "/home/s1987963/ds_group/Niklas/macparland_liver/sctransform/"
sctransform.path <- "/home/s1987963/ds_group/Niklas/macparland_liver/sctransform/alternative_clustering/"

### load liver data ###
liver <- readRDS(liver.path)

### SCtransform ###
options(future.globals.maxSize = 8000 * 1024^2)
liver.sctransform <- SCTransform(liver, verbose = TRUE)

### PCA ### 
liver.sctransform <- RunPCA(liver.sctransform, verbose = TRUE)

# elbow plot
sctransform.elbow.plot <- ElbowPlot(liver.sctransform, ndims = 50)
png(paste0(sctransform.path,"sctransform.elbow.plot.png"), width=1000,height=600,units="px")
print(sctransform.elbow.plot)
dev.off()

### visualize batch effect - pre-selection of dimensions ##
dims <- c(7,10,15,17,22) # uncorrected
for(d in dims){
  liver.sctransform <- RunUMAP(liver.sctransform, dims = 1:d, seed.use=1)
  umap.batch.plot <- DimPlot(liver.sctransform, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
  png(paste0(sctransform.path, "UMAP_dim", d, ".batch.png"), width=1400, height=800, units="px")
  print(umap.batch.plot)
  dev.off()
}

### clustering ###
## evaluate different numbers of PCs and resolutions ##
dims <- c(10,15,17)
res <- seq(0.3,1.5,0.1)
for(d in dims){
  liver.sctransform <- RunUMAP(liver.sctransform, dims=1:d, seed.use=1)
  batch.plot <- DimPlot(liver.sctransform, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
  for (r in res) {
    liver.sctransform <- FindNeighbors(liver.sctransform, dims = 1:d)
    liver.sctransform <- FindClusters(liver.sctransform, resolution = r)
    umap.plot <- DimPlot(liver.sctransform, reduction = "umap", label = F, pt.size = 0.1)
    eval.plot <- umap.plot + batch.plot
    png(paste0(sctransform.path, "UMAP_dim", d, "_res", r, ".png"), width=1400, height=600, units="px")
    print(eval.plot)
    dev.off()
  }
}

## preliminary clustering ##
# uncorrected/unadjusted: 17 first PCs, resolution 0.5
# uncorrected/unadjusted: 15 first PCs, resolution 0.5 alternative approach
# adjusted for nFeatures: XY first PCs, resolution XY
# adjusted for nCounts + percent.mito: XY first PCs, resolution
liver.sctransform <- FindNeighbors(liver.sctransform, dims = 1:15)
liver.sctransform <- FindClusters(liver.sctransform, resolution = 0.5)
# run UMAP
liver.sctransform <- RunUMAP(liver.sctransform, dims=1:15, seed.use=1)

### save data ###
saveRDS(liver.sctransform, file = paste0(sctransform.path, "macparland_liver_sctransform.rds"))

# set default essay to RNA counts
DefaultAssay(liver.sctransform) <- "RNA"
# Normalize RNA data for visualization purposes
liver.sctransform <- NormalizeData(liver.sctransform, verbose = TRUE)

## explore clustering at patient level ##
patient.clustering <- DimPlot(liver.sctransform, group.by = "seurat_clusters", split.by = "patient.ID", ncol = 9)
png(paste0(sctransform.path,"patient.clustering.png"), width=1600,height=400,units="px")
print(patient.clustering)
dev.off()

## QC metrics at cluster level ##
cluster.qc.heatmap <- FeaturePlot(liver.sctransform, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(sctransform.path,"cluster.qc.heatmap.png"), width=1500,height=500,units="px")
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
                 "CALCRL", # LSECs
                 "KRT19","EPCAM","FXYD2", #cholangiocytes
                 "ACTA2","COL1A1" # Hepatic Stellate Cells
)
# lineage markers
lineage.genes <- c("EPCAM", #epithelial cells
                   "HBB", #erythroid cells
                   "PTPRC", #immune cells
                   "CD27","IGHG1", # plasma cells
                   "GZMK","KLRF1", #NK cells
                   "TPSB2", #mast cells - maybe "IL1RL1"?
                   "CD3D", "GZMA", #T-cells
                   "CD79A", "CD79B" #B-cells
)

# feature plot with macrophage markers
macrophage.markers <- FeaturePlot(liver.sctransform, features = macrophage.genes, pt.size = 0.2, ncol = 5) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(sctransform.path,"macrophage.markers.png"), width=1600,height=800,units="px")
print(macrophage.markers)
dev.off()

# feature plot with monocyte markers
monocyte.markers <- FeaturePlot(liver.sctransform, features = monocyte.genes, pt.size = 0.2, ncol = 2) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(sctransform.path,"monocyte.markers.png"), width=800,height=800,units="px")
print(monocyte.markers)
dev.off()

# feature plot with DC markers
dc.markers <- FeaturePlot(liver.sctransform, features = dc.genes, pt.size = 0.2, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(sctransform.path,"dc.markers.png"), width=1200,height=800,units="px")
print(dc.markers)
dev.off()

# feature plot with liver markers
liver.markers <- FeaturePlot(liver.sctransform, features = liver.genes, pt.size = 0.2, ncol = 4) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(sctransform.path,"liver.markers.png"), width=1600,height=800,units="px")
print(liver.markers)
dev.off()

# feature plot with lineage markers
lineage.markers <- FeaturePlot(liver.sctransform, features = lineage.genes, pt.size = 0.2, ncol = 4) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(sctransform.path,"lineage.markers.png"), width=1600,height=1200,units="px")
print(lineage.markers)
dev.off()