library(Seurat)
library(harmony)
library(reticulate)
library(umap)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
options(future.globals.maxSize = 8000 * 1024^2)

### path.variables ###
liver.path <- "/home/s1987963/ds_group/Niklas/macparland_liver/macparland_liver_filtered.rds"
harmony.path <- "/home/s1987963/ds_group/Niklas/macparland_liver/harmony/"

### load data ###
liver <- readRDS(liver.path)

### normalisation ###
liver <- NormalizeData(liver, normalization.method = "LogNormalize", scale.factor = 10000)

### feature selection ###
# selecting the most highly variable genes
liver <- FindVariableFeatures(liver, selection.method = "vst", nfeatures = 2000)

### scale data ###
liver <- ScaleData(liver) # uncorrected
# liver <- ScaleData(liver, vars.to.regress = "nFeature_RNA") # adjusted for nFeatures
# liver <- ScaleData(liver, vars.to.regress = c("nCount_RNA", "percent.mt")) # adjusted for nCounts + percent.mt

### dimensionality reduction: PCA ###
liver <- RunPCA(liver, features = VariableFeatures(object = liver))

## PCA elbow plot ##
pca.elbow.plot <- ElbowPlot(liver, ndims = 50, reduction = "pca")
png(paste0(harmony.path,"pca.elbow.plot.png"), width=1000,height=600,units="px")
print(pca.elbow.plot)
dev.off()

## visualize PCs ##
DimPlot(object = liver, reduction = "pca", pt.size = .1, group.by = "patient.ID")
VlnPlot(object = liver, features = "PC_1", group.by = "patient.ID", pt.size = .1)

## run harmony ##
liver.harmony <- liver %>% RunHarmony("patient.ID", plot_convergence = TRUE)

## explore harmony results ##
DimPlot(object = liver.harmony, reduction = "harmony", pt.size = .1, group.by = "patient.ID")
VlnPlot(object = liver.harmony, features = "harmony_1", group.by = "patient.ID", pt.size = .1)

## harmony elbow plot ##
harmony.elbow.plot <- ElbowPlot(liver.harmony, ndims = 50, reduction = "harmony")
png(paste0(harmony.path,"harmony.elbow.plot.png"), width=1000,height=600,units="px")
print(harmony.elbow.plot)
dev.off()

### visualize batch effect - pre-selection of dimensions ##
dims <- c(6,8,9,12,14,17,24)
for(d in dims){
  liver.harmony <- RunUMAP(liver.harmony, reduction = "harmony", dims = 1:d, seed.use=1)
  umap.batch.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
  png(paste0(harmony.path, "UMAP_dim", d, ".batch.png"), width=1400, height=800, units="px")
  print(umap.batch.plot)
  dev.off()
}

## clustering ##
#dims <- c() #adjusted for nFeatures
# dims <- c() #adjusted for nCounts + percent.mito
# dims <- c(6,8,9,12,14,17) # uncorrected/adjusted
dims <- c(9,12,14,17)
res <- seq(0.5,1.5,0.1)
for(d in dims){
  liver.harmony <- RunUMAP(liver.harmony, reduction = "harmony", dims=1:d, seed.use=1)
  
  # plot batch effect
  batch.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
  for (r in res) {
    liver.harmony <- FindNeighbors(liver.harmony, reduction = "harmony", dims = 1:d)
    liver.harmony <- FindClusters(liver.harmony, reduction = "harmony", resolution = r)
    umap.plot <- DimPlot(liver.harmony, reduction = "umap", label = F, pt.size = 0.1)
    
    # create eval plot
    eval.plot <- umap.plot + batch.plot
    png(paste0(harmony.path, "UMAP_dim", d, "_res", r, ".png"), width=1500, height=600, units="px")
    print(eval.plot)
    dev.off()
  }
}

## preliminary clustering ##
# uncorrected/adjusted: 17 first PCs, resolution 1
# adjusted for nFeatures: XY first PCs, resolution XY
# adjusted for nCounts + percent.mito: XY first PCs, resolution
liver.harmony <- FindNeighbors(liver.harmony, reduction = "harmony", dims = 1:17)
liver.harmony <- FindClusters(liver.harmony, reduction = "harmony", resolution = 1)
# run UMAP
liver.harmony <- RunUMAP(liver.harmony, reduction = "harmony", dims=1:17, seed.use=1)

### save data ###
saveRDS(liver.harmony, file = paste0(harmony.path, "macparland_liver_harmony.rds"))

## explore clustering at patient level ##
patient.clustering <- DimPlot(liver.harmony, group.by = "seurat_clusters", split.by = "patient.ID", ncol = 9)
png(paste0(harmony.path,"patient.clustering.png"), width=1600,height=400,units="px")
print(patient.clustering)
dev.off()

## QC metrics at cluster level ##
cluster.qc.heatmap <- FeaturePlot(liver.harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"cluster.qc.heatmap.png"), width=1500,height=500,units="px")
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
macrophage.markers <- FeaturePlot(liver.harmony, features = macrophage.genes, pt.size = 0.2, ncol = 5) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"macrophage.markers.png"), width=1800,height=800,units="px")
print(macrophage.markers)
dev.off()

# feature plot with monocyte markers
monocyte.markers <- FeaturePlot(liver.harmony, features = monocyte.genes, pt.size = 0.2, ncol = 2) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"monocyte.markers.png"), width=800,height=800,units="px")
print(monocyte.markers)
dev.off()

# feature plot with DC markers
dc.markers <- FeaturePlot(liver.harmony, features = dc.genes, pt.size = 0.2, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"dc.markers.png"), width=1200,height=800,units="px")
print(dc.markers)
dev.off()

# feature plot with liver markers
liver.markers <- FeaturePlot(liver.harmony, features = liver.genes, pt.size = 0.2, ncol = 4) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"liver.markers.png"), width=1600,height=800,units="px")
print(liver.markers)
dev.off()

# feature plot with lineage markers
lineage.markers <- FeaturePlot(liver.harmony, features = lineage.genes, pt.size = 0.2, ncol = 4) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"lineage.markers.png"), width=1600,height=1200,units="px")
print(lineage.markers)
dev.off()
