library(Seurat)
library(harmony)
library(reticulate)
library(umap)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(patchwork)
options(future.globals.maxSize = 8000 * 1024^2)

### path variables ###
blood.path <- "/home/s1987963/ds_group/Niklas/reyes_blood/reyes_blood_healthy.rds"
harmony.path <- "/home/s1987963/ds_group/Niklas/reyes_blood/harmony/uncorrected/"

### load blood data ###
blood <- readRDS(blood.path)

### normalisation ###
blood <- NormalizeData(blood, normalization.method = "LogNormalize", scale.factor = 10000)

### feature selection ###
# selecting the most highly variable genes
blood <- FindVariableFeatures(blood, selection.method = "vst", nfeatures = 2000)

### scale data ###
blood <- ScaleData(blood) # uncorrected
# blood <- ScaleData(blood, vars.to.regress = "nFeature_RNA") # adjusted for nFeatures
# blood <- ScaleData(blood, vars.to.regress = c("nCount_RNA", "percent.mt")) # adjusted for nCounts + percent.mt

### dimensionality reduction: PCA ###
blood <- RunPCA(blood, features = VariableFeatures(object = blood))

## PCA elbow plot ##
pca.elbow.plot <- ElbowPlot(blood, ndims = 50, reduction = "pca")
png(paste0(harmony.path,"pca.elbow.plot.png"), width=1000,height=600,units="px")
print(pca.elbow.plot)
dev.off()

## visualise PCs ##
DimPlot(object = blood, reduction = "pca", pt.size = .1, group.by = "patient.ID")
VlnPlot(object = blood, features = "PC_1", group.by = "patient.ID", pt.size = .1)

## run harmony ##
blood.harmony <- blood %>% RunHarmony("patient.ID", plot_convergence = TRUE)

## explore harmony results ##
DimPlot(object = blood.harmony, reduction = "harmony", pt.size = .1, group.by = "patient.ID")
VlnPlot(object = blood.harmony, features = "harmony_1", group.by = "patient.ID", pt.size = .1)

## harmony elbow plot ##
harmony.elbow.plot <- ElbowPlot(blood.harmony, ndims = 50, reduction = "harmony")
png(paste0(harmony.path,"harmony.elbow.plot.png"), width=1000,height=600,units="px")
print(harmony.elbow.plot)
dev.off()

### visualize batch effect - pre-selection of dimensions ##
dims <- c(8,12,14,18,29,50)
for(d in dims){
  blood.harmony <- RunUMAP(blood.harmony, reduction = "harmony", dims = 1:d, seed.use=1)
  umap.batch.plot <- DimPlot(blood.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
  png(paste0(harmony.path, "UMAP_dim", d, ".batch.png"), width=1400, height=600, units="px")
  print(umap.batch.plot)
  dev.off()
}

## clustering ##
#dims <- c() #adjusted for nFeatures
# dims <- c() #adjusted for nCounts + percent.mito
# dims <- c(4,6,10,13,18) # uncorrected/adjusted
dims <- c(12,18,29,50)
res <- seq(0.1,0.5,0.1)
for(d in dims){
  blood.harmony <- RunUMAP(blood.harmony, reduction = "harmony", dims=1:d, seed.use=1)
  
  # plot batch effect
  batch.plot <- DimPlot(blood.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
  
  # plot original idents
  ident.plot <- DimPlot(blood.harmony, reduction = "umap", group.by = "cell_type", pt.size = 0.1)
  
  for (r in res) {
    blood.harmony <- FindNeighbors(blood.harmony, reduction = "harmony", dims = 1:d)
    blood.harmony <- FindClusters(blood.harmony, reduction = "harmony", resolution = r)
    umap.plot <- DimPlot(blood.harmony, reduction = "umap", label = F, pt.size = 0.1)
    
    # create triple plot
    eval.plot <- umap.plot + batch.plot + ident.plot
    png(paste0(harmony.path, "UMAP_dim", d, "_res", r, ".png"), width=1500, height=600, units="px")
    print(eval.plot)
    dev.off()
  }
}

### save data ###
saveRDS(blood.harmony, file = paste0(harmony.path, "reyes_blood_harmony.rds"))

## preliminary clustering ##
# uncorrected/adjusted: 29 first PCs, resolution 0.1
# adjusted for nFeatures: XY first PCs, resolution XY
# adjusted for nCounts + percent.mito: XY first PCs, resolution
blood.harmony <- FindNeighbors(blood.harmony, reduction = "harmony", dims = 1:29)
blood.harmony <- FindClusters(blood.harmony, reduction = "harmony", resolution = 0.1)
# run UMAP
blood.harmony <- RunUMAP(blood.harmony, reduction = "harmony", dims=1:29, seed.use=1)

### save data ###
blood.harmony$organ <- "blood"
saveRDS(blood.harmony, file = paste0(harmony.path, "reyes_blood_harmony.rds"))

## explore clustering at patient level ##
patient.clustering <- DimPlot(blood.harmony, group.by = "seurat_clusters", split.by = "patient.ID", ncol = 9)
png(paste0(harmony.path,"patient.clustering.png"), width=1600,height=800,units="px")
print(patient.clustering)
dev.off()

## explore original clustering ##
celltype.clustering <- DimPlot(blood.harmony, group.by = "cell_type")
png(paste0(harmony.path,"celltype.clustering.png"), width=1600,height=800,units="px")
print(celltype.clustering)
dev.off()

## QC metrics at cluster level ##
cluster.qc.heatmap <- FeaturePlot(blood.harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
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
# lineage markers
lineage.genes <- c("PTPRC", #immune cells
                   "CD3D", "GZMA", #T-cells
                   "CD79A", "CD79B" #B-cells
)

# feature plot with macrophage markers
macrophage.markers <- FeaturePlot(blood.harmony, features = macrophage.genes, pt.size = 0.2, ncol = 5) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"macrophage.markers.png"), width=1800,height=800,units="px")
print(macrophage.markers)
dev.off()

# feature plot with monocyte markers
monocyte.markers <- FeaturePlot(blood.harmony, features = monocyte.genes, pt.size = 0.2, ncol = 2) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"monocyte.markers.png"), width=800,height=800,units="px")
print(monocyte.markers)
dev.off()

# feature plot with DC markers
dc.markers <- FeaturePlot(blood.harmony, features = dc.genes, pt.size = 0.2, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"dc.markers.png"), width=1200,height=800,units="px")
print(dc.markers)
dev.off()

# feature plot with more general lineage markers
lineage.markers <- FeaturePlot(blood.harmony, features = lineage.genes, pt.size = 0.2, ncol = 5) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"lineage.markers.png"), width=1500,height=500,units="px")
print(lineage.markers)
dev.off()