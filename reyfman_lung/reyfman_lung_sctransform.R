library(Seurat)
library(sctransform)
library(reticulate)
library(umap)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

### path variables ###
sctransform.path <- "/home/s1987963/ds_group/Niklas/reyfman_lung/sctransform/"

### load lung data ###
lung <- readRDS("/home/s1987963/ds_group/Niklas/reyfman_lung/reyfman_lung_healthy.rds")

### SCtransform ###
options(future.globals.maxSize = 8000 * 1024^2)
lung.sctransform <- SCTransform(lung, verbose = TRUE)

### PCA ### 
lung.sctransform <- RunPCA(lung.sctransform, verbose = TRUE)

# elbow plot
sctransform.elbow.plot <- ElbowPlot(lung.sctransform, ndims = 50)
png(paste0(sctransform.path,"sctransform.elbow.plot.png"), width=1000,height=600,units="px")
print(sctransform.elbow.plot)
dev.off()

### clustering ###
dims <- c(6,8,9,11,13)
res <- seq(0.2,1.5,0.1)
for(d in dims){
  lung.sctransform <- RunUMAP(lung.sctransform, dims=1:d, seed.use=1)
  for (r in res) {
    lung.sctransform <- FindNeighbors(lung.sctransform, dims = 1:d)
    lung.sctransform <- FindClusters(lung.sctransform, resolution = r)
    umap.plot <- DimPlot(lung.sctransform, reduction = "umap", label = F, pt.size = 0.1)
    batch.plot <- DimPlot(lung.sctransform, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
    eval.plot <- umap.plot + batch.plot
    png(paste0(sctransform.path, "UMAP_dim", d, "_res", r, ".png"), width=1400, height=600, units="px")
    print(eval.plot)
    dev.off()
  }
}

## preliminary clustering ##
# first 8 PCs, resolution 0.6
lung.sctransform <- FindNeighbors(lung.sctransform, dims = 1:8)
lung.sctransform <- FindClusters(lung.sctransform, resolution = 0.6)
# run UMAP
lung.sctransform <- RunUMAP(lung.sctransform, dims=1:8, seed.use=1)
# run tSNE
lung.sctransform <- RunTSNE(lung.sctransform, dims=1:8, seed.use=1)

### save data ###
saveRDS(lung.sctransform, file = paste0(sctransform.path, "raredon_lung_sctransform.rds"))

## explore clustering at patient level ##
patient.clustering <- DimPlot(lung.sctransform, group.by = "seurat_clusters", split.by = "patient.ID", ncol = 4)
png(paste0(sctransform.path,"patient.clustering.png"), width=1600,height=800,units="px")
print(patient.clustering)
dev.off()

## triple UMAP ##
#umap.1 <- DimPlot(lung.sctransform, reduction = "umap", label = F, pt.size = 0.1)
#umap.2 <- DimPlot(lung.sctransform, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
#umap.3 <- DimPlot(lung.sctransform, reduction = "umap", group.by = "scrublet_auto", pt.size = 0.5)
#triple.umap <- umap.1 + umap.2 + umap.3
#png(paste0(sctransform.path,"triple.umap.png"), width=1800,height=600,units="px")
#print(triple.umap)
#dev.off()

## QC metrics at cluster level ##
cluster.qc.heatmap <- FeaturePlot(lung.sctransform, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(sctransform.path,"cluster.qc.heatmap.png"), width=1500,height=500,units="px")
print(cluster.qc.heatmap)
dev.off()

## marker gene visualization ##
# macrophage markers
macrophage.genes <- c("CSF1R", "LYZ", "HLA-DRA", "ITGAX", "ITGAM", "C1QB","MRC1", "MARCO") #MSR1 maybe?
# monocyte markers
monocyte.genes <- c("CD14", "MNDA", "S100A8","S100A9")
# dendritic cell markers
dc.genes <- c("CD1C","XCR1", "CD86")
# lineage markers
lineage.genes <- c("EPCAM", #epithelial cells
                   "CDH5", "PECAM1", #endothelial cells
                   "PDGFRB", "PDGFRA", "CD34", #mesenchymal cells
                   "PTPRC", #immune cells
                   "CD3D", "GZMA", #T-cells
                   "CD79A", "CD79B" #B-cells
)

# set default essay to RNA counts
DefaultAssay(lung.sctransform) <- "RNA"
# Normalize RNA data for visualization purposes
lung.sctransform <- NormalizeData(lung.sctransform, verbose = TRUE)

# feature plot with macrophage markers
macrophage.markers <- FeaturePlot(lung.sctransform, features = macrophage.genes, pt.size = 0.2, ncol = 4) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(sctransform.path,"macrophage.markers.png"), width=1600,height=800,units="px")
print(macrophage.markers)
dev.off()

# feature plot with monocyte markers
monocyte.markers <- FeaturePlot(lung.sctransform, features = monocyte.genes, pt.size = 0.2, ncol = 2) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(sctransform.path,"monocyte.markers.png"), width=800,height=800,units="px")
print(monocyte.markers)
dev.off()

# feature plot with DC markers
dc.markers <- FeaturePlot(lung.sctransform, features = dc.genes, pt.size = 0.2, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(sctransform.path,"dc.markers.png"), width=1200,height=400,units="px")
print(dc.markers)
dev.off()

# feature plot with more general lineage markers
lineage.markers <- FeaturePlot(lung.sctransform, features = lineage.genes, pt.size = 0.2, ncol = 5) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(sctransform.path,"lineage.markers.png"), width=2000,height=800,units="px")
print(lineage.markers)
dev.off()
