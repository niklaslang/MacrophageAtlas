library(Seurat)
library(sctransform)
library(reticulate)
library(umap)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

### load lung data ###
sctransform.path <- "/home/s1987963/ds_group/Niklas/raredon_lung/sctransform/"
sctransform.uncorrected.path <- "/home/s1987963/ds_group/Niklas/raredon_lung/sctransform/uncorrected/"
sctransform.regress.nfeatures.path <- "/home/s1987963/ds_group/Niklas/raredon_lung/sctransform/regress.nFeatures/"
lung <- readRDS(paste0("/home/s1987963/ds_group/Niklas/raredon_lung/raredon_lung.rds"))

### SCtransform ###
options(future.globals.maxSize = 8000 * 1024^2)
lung.sctransform <- SCTransform(lung, verbose = TRUE)
lung.sctransform.regress.nfeatures <- SCTransform(lung,vars.to.regress = "nFeature_RNA", verbose = TRUE)

### PCA ### 
lung.sctransform <- RunPCA(lung.sctransform, verbose = TRUE)
lung.sctransform.regress.nfeatures <- RunPCA(lung.sctransform.regress.nfeatures, verbose = TRUE)

# elbow plot
sctransform.elbow.plot <- ElbowPlot(lung.sctransform, ndims = 50)
png(paste0(sctransform.path,"sctransform.elbow.plot.png"), width=1000,height=600,units="px")
print(sctransform.elbow.plot)
dev.off()

sctransform.regress.nfeatures.elbow.plot <- ElbowPlot(lung.sctransform.regress.nfeatures, ndims = 50)
png(paste0(sctransform.path,"sctransform.regress.nfeatures.elbow.plot.png"), width=1000,height=600,units="px")
print(sctransform.regress.nfeatures.elbow.plot)
dev.off()

### SCTRANSFORM.uncorrected ###
### clustering ###
## evaluate different numbers of PCs and resolutions ##
dims <- c(8,10,15,17,25)
res <- seq(0.2,1.5,0.1)
for(d in dims){
  lung.sctransform <- RunUMAP(lung.sctransform, dims=1:d, seed.use=1)
  for (r in res) {
    lung.sctransform <- FindNeighbors(lung.sctransform, dims = 1:d)
    lung.sctransform <- FindClusters(lung.sctransform, resolution = r)
    umap.plot <- DimPlot(lung.sctransform, reduction = "umap", label = F, pt.size = 0.1)
    batch.plot <- DimPlot(lung.sctransform, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
    eval.plot <- umap.plot + batch.plot
    png(paste0(sctransform.uncorrected.path, "UMAP_dim", d, "_res", r, ".png"), width=1400, height=600, units="px")
    print(eval.plot)
    dev.off()
  }
}

## preliminary clustering ##
# first 10 PCs, resolution 0.7
lung.sctransform <- FindNeighbors(lung.sctransform, dims = 1:10)
lung.sctransform <- FindClusters(lung.sctransform, resolution = 0.7)
# run UMAP
lung.sctransform <- RunUMAP(lung.sctransform, dims=1:10, seed.use=1)
# run tSNE
lung.sctransform <- RunTSNE(lung.sctransform, dims=1:10, seed.use=1)

### save data ###
saveRDS(lung.sctransform, file = paste0(sctransform.uncorrected.path, "raredon_lung_sctransform.uncorrected.rds"))

## explore clustering at patient level ##
patient.clustering <- DimPlot(lung.sctransform, group.by = "seurat_clusters", split.by = "patient.ID", ncol = 7)
png(paste0(sctransform.uncorrected.path,"patient.clustering.png"), width=1800,height=600,units="px")
print(patient.clustering)
dev.off()

## triple UMAP ##
umap.1 <- DimPlot(lung.sctransform, reduction = "umap", label = F, pt.size = 0.1)
umap.2 <- DimPlot(lung.sctransform, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
umap.3 <- DimPlot(lung.sctransform, reduction = "umap", group.by = "scrublet_auto", pt.size = 0.5)
triple.umap <- umap.1 + umap.2 + umap.3
png(paste0(sctransform.uncorrected.path,"triple.umap.png"), width=1800,height=600,units="px")
print(triple.umap)
dev.off()

## QC metrics at cluster level ##
cluster.qc.heatmap <- FeaturePlot(lung.sctransform, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(sctransform.uncorrected.path,"cluster.qc.heatmap.png"), width=1500,height=500,units="px")
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
png(paste0(sctransform.uncorrected.path,"macrophage.markers.png"), width=1600,height=800,units="px")
print(macrophage.markers)
dev.off()

# feature plot with monocyte markers
monocyte.markers <- FeaturePlot(lung.sctransform, features = monocyte.genes, pt.size = 0.2, ncol = 2) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(sctransform.uncorrected.path,"monocyte.markers.png"), width=800,height=800,units="px")
print(monocyte.markers)
dev.off()

# feature plot with DC markers
dc.markers <- FeaturePlot(lung.sctransform, features = dc.genes, pt.size = 0.2, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(sctransform.uncorrected.path,"dc.markers.png"), width=1200,height=400,units="px")
print(dc.markers)
dev.off()

# feature plot with more general lineage markers
lineage.markers <- FeaturePlot(lung.sctransform, features = lineage.genes, pt.size = 0.2, ncol = 5) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(sctransform.uncorrected.path,"lineage.markers.png"), width=2000,height=800,units="px")
print(lineage.markers)
dev.off()

### SCTRANSFORM.regress.nFeatures ###
### clustering ###
## evaluate different numbers of PCs and resolutions ##
dims <- c(8,10,15,23)
res <- seq(0.5,1.5,0.1)
for(d in dims){
  lung.sctransform.regress.nfeatures <- RunUMAP(lung.sctransform.regress.nfeatures, dims=1:d, seed.use=1)
  for (r in res) {
    lung.sctransform.regress.nfeatures <- FindNeighbors(lung.sctransform.regress.nfeatures, dims = 1:d)
    lung.sctransform.regress.nfeatures <- FindClusters(lung.sctransform.regress.nfeatures, resolution = r)
    umap.plot <- DimPlot(lung.sctransform.regress.nfeatures, reduction = "umap", label = F, pt.size = 0.1)
    batch.plot <- DimPlot(lung.sctransform.regress.nfeatures, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
    eval.plot <- umap.plot + batch.plot
    png(paste0(sctransform.regress.nfeatures.path, "UMAP_dim", d, "_res", r, ".png"), width=1400, height=600, units="px")
    print(eval.plot)
    dev.off()
  }
}

## preliminary clustering ##
# first 10 PCs, resolution 1.0
lung.sctransform.regress.nfeatures <- FindNeighbors(lung.sctransform.regress.nfeatures, dims = 1:10)
lung.sctransform.regress.nfeatures <- FindClusters(lung.sctransform.regress.nfeatures, resolution = 1.0)
# run UMAP
lung.sctransform.regress.nfeatures <- RunUMAP(lung.sctransform.regress.nfeatures, dims=1:10, seed.use=1)
# run tSNE
lung.sctransform.regress.nfeatures <- RunTSNE(lung.sctransform.regress.nfeatures, dims=1:10, seed.use=1)

### save data ###
saveRDS(lung.sctransform.regress.nfeatures, file = paste0(sctransform.regress.nfeatures.path, "raredon_lung_sctransform.regress.nfeatures.rds"))

## triple UMAP ##
umap.1 <- DimPlot(lung.sctransform.regress.nfeatures, reduction = "umap", label = F, pt.size = 0.1)
umap.2 <- DimPlot(lung.sctransform.regress.nfeatures, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
umap.3 <- DimPlot(lung.sctransform.regress.nfeatures, reduction = "umap", group.by = "scrublet_auto", pt.size = 0.5)
triple.umap <- umap.1 + umap.2 + umap.3
png(paste0(sctransform.regress.nfeatures.path,"triple.umap.png"), width=1800,height=600,units="px")
print(triple.umap)
dev.off()

## explore clustering at patient level ##
patient.clustering <- DimPlot(lung.sctransform.regress.nfeatures, group.by = "seurat_clusters", split.by = "patient.ID", ncol = 7)
png(paste0(sctransform.regress.nfeatures.path,"patient.clustering.png"), width=1800,height=600,units="px")
print(patient.clustering)
dev.off()

## QC metrics at cluster level ##
cluster.qc.heatmap <- FeaturePlot(lung.sctransform.regress.nfeatures, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(sctransform.regress.nfeatures.path,"cluster.qc.heatmap.png"), width=1500,height=500,units="px")
print(cluster.qc.heatmap)
dev.off()

# set default essay to RNA counts
DefaultAssay(lung.sctransform.regress.nfeatures) <- "RNA"
# Normalize RNA data for visualization purposes
lung.sctransform.regress.nfeatures <- NormalizeData(lung.sctransform.regress.nfeatures, verbose = TRUE)

# feature plot with macrophage markers
macrophage.markers <- FeaturePlot(lung.sctransform.regress.nfeatures, features = macrophage.genes, pt.size = 0.2, ncol = 4) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(sctransform.regress.nfeatures.path,"macrophage.markers.png"), width=1600,height=800,units="px")
print(macrophage.markers)
dev.off()

# feature plot with monocyte markers
monocyte.markers <- FeaturePlot(lung.sctransform.regress.nfeatures, features = monocyte.genes, pt.size = 0.2, ncol = 2) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(sctransform.regress.nfeatures.path,"monocyte.markers.png"), width=800,height=800,units="px")
print(monocyte.markers)
dev.off()

# feature plot with DC markers
dc.markers <- FeaturePlot(lung.sctransform.regress.nfeatures, features = dc.genes, pt.size = 0.2, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(sctransform.regress.nfeatures.path,"dc.markers.png"), width=1200,height=400,units="px")
print(dc.markers)
dev.off()

# feature plot with more general lineage markers
lineage.markers <- FeaturePlot(lung.sctransform.regress.nfeatures, features = lineage.genes, pt.size = 0.2, ncol = 5) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(sctransform.regress.nfeatures.path,"lineage.markers.png"), width=2000,height=800,units="px")
print(lineage.markers)
dev.off()