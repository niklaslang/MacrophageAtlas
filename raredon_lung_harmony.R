library(Seurat)
library(harmony)
library(reticulate)
library(umap)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(patchwork)

### path variables ###
harmony.path <- "/home/s1987963/ds_group/Niklas/raredon_lung/harmony/"
# harmony.path <- "/Users/Niklas/Documents/MASTER/harmony/raredon_lung/"

### load lung data ###
lung <- readRDS("/home/s1987963/ds_group/Niklas/raredon_lung/raredon_lung.rds")
lung.harmony <- readRDS("/home/s1987963/ds_group/Niklas/raredon_lung/harmony/raredon_lung_harmony.rds")

### normalisation ###
lung <- NormalizeData(lung, normalization.method = "LogNormalize", scale.factor = 10000)

### feature selection ###
# selecting the most highly variable genes
lung <- FindVariableFeatures(lung, selection.method = "vst", nfeatures = 2000)

### scale data ###
lung <- ScaleData(lung)
# lung <- ScaleData(lung,vars.to.regress = "nFeature_RNA")

### dimensionality reduction: PCA ###
lung <- RunPCA(lung, features = VariableFeatures(object = lung))

## PCA elbow plot ##
pca.elbow.plot <- ElbowPlot(lung, ndims = 50, reduction = "pca")
png(paste0(harmony.path,"pca.elbow.plot.png"), width=1000,height=600,units="px")
print(pca.elbow.plot)
dev.off()

## visualise PCs ##
DimPlot(object = lung, reduction = "pca", pt.size = .1, group.by = "patient.ID")
VlnPlot(object = lung, features = "PC_1", group.by = "patient.ID", pt.size = .1)

## run harmony ##
lung.harmony <- RunHarmony(lung, "patient.ID", plot_convergence = TRUE)

## explore harmony results ##
DimPlot(object = lung.harmony, reduction = "harmony", pt.size = .1, group.by = "patient.ID")
VlnPlot(object = lung.harmony, features = "harmony_1", group.by = "patient.ID", pt.size = .1)

## harmony elbow plot ##
harmony.elbow.plot <- ElbowPlot(lung.harmony, ndims = 50, reduction = "harmony")
png(paste0(harmony.path,"harmony.elbow.plot.png"), width=1000,height=600,units="px")
print(harmony.elbow.plot)
dev.off()


### clustering ###
dims <- c(6,8,10,14,16,19)
# dims <- c(5,8,10,13,16,19) # adjusted for nFeature
## pre-selection of dimensions - visualize batch effect ##
for(d in dims){
  lung.harmony <- RunUMAP(lung.harmony, dims = 1:d, seed.use=1)
  umap.batch.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
  png(paste0(harmony.path, "UMAP_dim", d, "batch.png"), width=1400, height=600, units="px")
  print(umap.batch.plot)
  dev.off()
}

### clustering ###
## evaluate different numbers of PCs and resolutions ##
dims <- c(6,8,10) # uncorrected
res <- seq(0.5,1.5,0.1)
for(d in dims){
  lung.harmony <- RunUMAP(lung.harmony, reduction = "harmony", dims=1:d, seed.use=1)
  for (r in res) {
    lung.harmony <- FindNeighbors(lung.harmony, reduction = "harmony", dims = 1:d)
    lung.harmony <- FindClusters(lung.harmony, reduction = "harmony", resolution = r)
    umap.plot <- DimPlot(lung.harmony, reduction = "umap", label = F, pt.size = 0.1)
    batch.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
    eval.plot <- umap.plot + batch.plot
    png(paste0(harmony.path, "UMAP_dim", d, "_res", r, ".png"), width=1400, height=600, units="px")
    print(eval.plot)
    dev.off()
  }
}

## preliminary clustering ##
# uncorrected:first 10 PCs, resolution 1.4
# regress.out: first 8 PCs, resolutioin 1.0
lung.harmony <- FindNeighbors(lung.harmony, reduction = "harmony", dims = 1:10)
lung.harmony <- FindClusters(lung.harmony, reduction = "harmony", resolution = 1.4)
# run UMAP
lung.harmony <- RunUMAP(lung.harmony, reduction = "harmony", dims=1:10, seed.use=1)
# run tSNE
lung.harmony <- RunTSNE(lung.harmony, reduction = "harmony", dims=1:10, seed.use=1)

### save data ###
saveRDS(lung.harmony, file = paste0(harmony.path, "raredon_lung_harmony.rds"))

## explore clustering at patient level ##
patient.clustering <- DimPlot(lung.harmony, group.by = "seurat_clusters", split.by = "patient.ID", ncol = 7)
png(paste0(harmony.path,"patient.clustering.png"), width=1800,height=600,units="px")
print(patient.clustering)
dev.off()

## triple UMAP ##
umap.1 <- DimPlot(lung.harmony, reduction = "umap", label = F, pt.size = 0.1)
umap.2 <- DimPlot(lung.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
umap.3 <- DimPlot(lung.harmony, reduction = "umap", group.by = "scrublet_auto", pt.size = 0.5)
triple.umap <- umap.1 + umap.2 + umap.3
png(paste0(harmony.path,"triple.umap.png"), width=1800,height=600,units="px")
print(triple.umap)
dev.off()

## QC metrics at cluster level ##
cluster.qc.heatmap <- FeaturePlot(lung.harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"cluster.qc.heatmap.png"), width=1500,height=500,units="px")
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
                   "TPSB2", "IL1RL1",#mast cells - maybe "IL1RL1"?
                   "CD3D", "GZMA", #T-cells
                   "CD79A", "CD79B" #B-cells
)

# feature plot with macrophage markers
macrophage.markers <- FeaturePlot(lung.harmony, features = macrophage.genes, pt.size = 0.2, ncol = 4) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"macrophage.markers.png"), width=1600,height=800,units="px")
print(macrophage.markers)
dev.off()

# feature plot with monocyte markers
monocyte.markers <- FeaturePlot(lung.harmony, features = monocyte.genes, pt.size = 0.2, ncol = 2) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"monocyte.markers.png"), width=800,height=800,units="px")
print(monocyte.markers)
dev.off()

# feature plot with DC markers
dc.markers <- FeaturePlot(lung.harmony, features = dc.genes, pt.size = 0.2, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"dc.markers.png"), width=1200,height=800,units="px")
print(dc.markers)
dev.off()

# feature plot with more general lineage markers
lineage.markers <- FeaturePlot(lung.harmony, features = lineage.genes, pt.size = 0.2, ncol = 5) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"lineage.markers.png"), width=1800,height=1200,units="px")
print(lineage.markers)
dev.off()
