library(Seurat)
library(reticulate)
library(umap)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

### path.variables ###
blood.path <- "/home/s1987963/ds_group/Niklas/reyes_blood/reyes_blood_healthy.rds"
logtransform.path <- "/home/s1987963/ds_group/Niklas/reyes_blood/logtransform/uncorrected/"
# logtransform.path <- "/home/s1987963/ds_group/Niklas/reyes_blood/logtransform/uncorrected/" # corrected for nFeatures
# logtransform.path <- "/home/s1987963/ds_group/Niklas/reyes_blood/logtransform/regress.nCounts.mito/" # corrected for nCounts and mitochondrial fraction

### load data ###
blood <- readRDS(blood.path)

### normalisation ###
blood.logtransform <- NormalizeData(blood, normalization.method = "LogNormalize", scale.factor = 10000)

### feature selection ###
blood.logtransform <- FindVariableFeatures(blood.logtransform, selection.method = "vst", nfeatures = 2000)

# 10 most highly variable genes
top10 <- head(VariableFeatures(blood.logtransform), 10)

# plot variable features with and without labels
HVGs.plot <- VariableFeaturePlot(blood.logtransform)
HVGs.plot <- LabelPoints(plot = HVGs.plot, points = top10, repel = TRUE)
HVGs.plot

### scale data ###
blood.logtransform <- ScaleData(blood.logtransform)
# blood.logtransform <- ScaleData(blood.logtransform, vars.to.regress = "nFeature_RNA") # adjusted for nFeatures
# blood.logtransform <- ScaleData(blood.logtransform, vars.to.regress = c("nCount_RNA", "percent.mt")) # adjusted for nCounts + percent.mt

### dimensionality reduction: PCA ###
blood.logtransform <- RunPCA(blood.logtransform, features = VariableFeatures(object = blood.logtransform))

## elbow plot ##
logtransform.elbow.plot <- ElbowPlot(blood.logtransform, ndims = 50)
png(paste0(logtransform.path,"logtransform.elbow.plot.png"), width=1000,height=600,units="px")
print(logtransform.elbow.plot)
dev.off()

### visualize batch effect - pre-selection of dimensions ##
dims <- c(6,9,12,14,17,21,30,50)
for(d in dims){
  blood.logtransform <- RunUMAP(blood.logtransform, dims = 1:d, seed.use=1)
  umap.batch.plot <- DimPlot(blood.logtransform, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
  png(paste0(logtransform.path, "UMAP_dim", d, ".batch.png"), width=1400, height=600, units="px")
  print(umap.batch.plot)
  dev.off()
}

### clustering ###
## evaluate different numbers of PCs and resolutions ##
### clustering ###
# evaluate different numbers of PCs and resolutions
dims <- c(6,9,12,14,30,50) #adjusted for nFeatures
# dims <- c(4,6,9,11,15,17,25) #adjusted for nCounts + percent.mito
#dims <- c(6,8,9,11,15,20,30) # uncorrected/adjusted
res <- seq(0.1,1,0.1)
for(d in dims){
  blood.logtransform <- RunUMAP(blood.logtransform, dims=1:d, seed.use=1)
  for (r in res) {
    blood.logtransform <- FindNeighbors(blood.logtransform, dims = 1:d)
    blood.logtransform <- FindClusters(blood.logtransform, resolution = r)
    umap.plot <- DimPlot(blood.logtransform, reduction = "umap", label = F, pt.size = 0.1)
    batch.plot <- DimPlot(blood.logtransform, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
    eval.plot <- umap.plot + batch.plot
    png(paste0(logtransform.path, "UMAP_dim", d, "_res", r, ".png"), width=1400, height=600, units="px")
    print(eval.plot)
    dev.off()
  }
}

### definite clustering ###
## best (preliminary) clustering ##
# uncorrected/adjusted: 9 first PCs, resolution 0.1
# adjusted for nFeatures: XY first PCs, resolution XY
# adjusted for nCounts + percent.mito: XY first PCs, resolution XY
blood.logtransform <- FindNeighbors(blood.logtransform, dims = 1:9)
blood.logtransform <- FindClusters(blood.logtransform, resolution = 0.1)
# add UMAP
blood.logtransform <- RunUMAP(blood.logtransform, dims=1:9, seed.use=1)
# add tSNE
blood.logtransform <- RunTSNE(blood.logtransform, dims=1:9, seed.use=1)

### save data ###
blood.logtransform$organ <- "blood"
saveRDS(blood.logtransform, paste0(logtransform.path, "reyes_blood_logtransform.rds"))

## explore clustering at patient level ##
patient.clustering <- DimPlot(blood.logtransform, group.by = "seurat_clusters", split.by = "patient.ID", ncol = 9)
png(paste0(logtransform.path,"patient.clustering.png"), width=1600,height=800,units="px")
print(patient.clustering)
dev.off()

## explore original clustering at##
celltype.clustering <- DimPlot(blood.logtransform, group.by = "cell_type")
png(paste0(logtransform.path,"celltype.clustering.png"), width=1600,height=800,units="px")
print(celltype.clustering)
dev.off()

## QC metrics at cluster level ##
cluster.qc.heatmap <- FeaturePlot(blood.logtransform, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
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
# lineage markers
lineage.genes <- c("PTPRC", #immune cells
                   "CD3D", "GZMA", #T-cells
                   "CD79A", "CD79B" #B-cells
)

# feature plot with macrophage markers
macrophage.markers <- FeaturePlot(blood.logtransform, features = macrophage.genes, pt.size = 0.2, ncol = 5) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(logtransform.path,"macrophage.markers.png"), width=1600,height=800,units="px")
print(macrophage.markers)
dev.off()

# feature plot with monocyte markers
monocyte.markers <- FeaturePlot(blood.logtransform, features = monocyte.genes, pt.size = 0.2, ncol = 2) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(logtransform.path,"monocyte.markers.png"), width=800,height=800,units="px")
print(monocyte.markers)
dev.off()

# feature plot with DC markers
dc.markers <- FeaturePlot(blood.logtransform, features = dc.genes, pt.size = 0.2, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(logtransform.path,"dc.markers.png"), width=1200,height=800,units="px")
print(dc.markers)
dev.off()

# feature plot with more general lineage markers
lineage.markers <- FeaturePlot(blood.logtransform, features = lineage.genes, pt.size = 0.2, ncol = 5) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(logtransform.path,"lineage.markers.png"), width=1500,height=500,units="px")
print(lineage.markers)
dev.off()