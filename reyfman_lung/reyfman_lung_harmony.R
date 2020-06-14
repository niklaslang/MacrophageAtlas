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
lung.path <- "/home/s1987963/ds_group/Niklas/reyfman_lung/reyfman_lung_healthy.rds"
harmony.path <- "/home/s1987963/ds_group/Niklas/reyfman_lung/harmony/dim10_annotation/"

### load lung data ###
lung <- readRDS(lung.path)

### normalisation ###
lung <- NormalizeData(lung, normalization.method = "LogNormalize", scale.factor = 10000)

### feature selection ###
# selecting the most highly variable genes
lung <- FindVariableFeatures(lung, selection.method = "vst", nfeatures = 2000)

### scale data ###
lung <- ScaleData(lung) # uncorrected
# lung <- ScaleData(lung, vars.to.regress = "nFeature_RNA") # adjusted for nFeatures
# lung <- ScaleData(lung, vars.to.regress = c("nCount_RNA", "percent.mt")) # adjusted for nCounts + percent.mt

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
lung.harmony <- lung %>% RunHarmony("patient.ID", plot_convergence = TRUE)

## explore harmony results ##
DimPlot(object = lung.harmony, reduction = "harmony", pt.size = .1, group.by = "patient.ID")
VlnPlot(object = lung.harmony, features = "harmony_1", group.by = "patient.ID", pt.size = .1)

## harmony elbow plot ##
harmony.elbow.plot <- ElbowPlot(lung.harmony, ndims = 50, reduction = "harmony")
png(paste0(harmony.path,"harmony.elbow.plot.png"), width=1000,height=600,units="px")
print(harmony.elbow.plot)
dev.off()

## clustering ##
#dims <- c() #adjusted for nFeatures
# dims <- c() #adjusted for nCounts + percent.mito
# dims <- c(4,6,10,13,18) # uncorrected/adjusted
dims <- c(4,6,10,13,18)
res <- seq(0.3,1.5,0.1)
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
# uncorrected/adjusted: 10 first PCs, resolution 0.4
# adjusted for nFeatures: XY first PCs, resolution XY
# adjusted for nCounts + percent.mito: XY first PCs, resolution
lung.harmony <- FindNeighbors(lung.harmony, reduction = "harmony", dims = 1:10)
lung.harmony <- FindClusters(lung.harmony, reduction = "harmony", resolution = 0.4)
# run UMAP
lung.harmony <- RunUMAP(lung.harmony, reduction = "harmony", dims=1:10, seed.use=1)
# run tSNE
lung.harmony <- RunTSNE(lung.harmony, reduction = "harmony", dims=1:10, seed.use=1)

### save data ###
saveRDS(lung.harmony, file = paste0(harmony.path, "raredon_lung_harmony.rds"))

## explore clustering at patient level ##
patient.clustering <- DimPlot(lung.harmony, group.by = "seurat_clusters", split.by = "patient.ID", pt.size = 0.2, ncol = 4)
png(paste0(harmony.path,"patient.clustering.png"), width=1600,height=800,units="px")
print(patient.clustering)
dev.off()

## QC metrics at cluster level ##
cluster.qc.heatmap <- FeaturePlot(lung.harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"cluster.qc.heatmap.png"), width=1500,height=500,units="px")
print(cluster.qc.heatmap)
dev.off()

## marker gene visualization ##
#MNP genes
MNP.genes <- c("CD14", "FCGR3A","CSF1R", "CD68", "LYZ", 
               "CCR2", "CX3CR1", "S100A8", "MARCO", "MNDA", 
               "TIMD4", "CD163", "C1QB", "F13A1",	"MCEMP1",	
               "TREM2",	"CD9",	"VCAN")
# cDC 1 genes
cDC1.genes <- c("XCR1", "CLEC10A", "FCER1A", "CLEC9A")

# cDC 2 genes 
cDC2.genes <- c("CD1C","HLA-DRA",	"HLA-DRB1",	"CCR7",	"ITGAX",
                "CD1E")

# pDC genes 
pDC.genes <- c("LILRA4", "IRF8",	"IRF7",	"CLEC4C")

# T cells
Tcell.genes <- c("CD3D",	"CD4",	"CD8A", "IL7R", "CD2", 
                 "CCR7",	"LTB", "TRAC")

# NK cells
NK.genes <- c("GZMA",	"NKG7",	"XCL1",	"FGFBP2",	"GNLY",	
              "GZMK",	"KLRF1",	"KLRC1")

# B cells
Bcell.genes <- c("CD79A", "CD79B",	"MS4A1")

# plasma cells
plasmacell.genes <- c("IGHG1", "MZB1",	"IGKC",	"JCHAIN",	"JSRP1")

# mast cells
mastcell.genes <- c("TPSB2",	"TPSAB1",	"CPA3",	"MS4A2")

# erys
ery.genes <- c("HBA1", "HBB")

# epithelial cells
epithelial.genes <- c("EPCAM",	"KRT19", "ALB", "TF",	"HNF4A",	"SFTPC")

# endothelial cells
endothelial.genes <- c("KDR",	"CD34", "VWF", "CLDN5", "PECAM1", 
                       "ICAM2",	"PDPN",	"CLEC14A",	"MMRN1")

# mesenchymal cells
mesenchymal.genes <- c("COL1A1",	"COL3A1", "ACTA2", "MYH11",	"PDGFRA", 
                       "PDGFRB", "RGS5",	"MSLN",	"DCN")

# proliferating cells
proliferation.genes <- c("MKI67",	"TOP2A")

# feature plot with MNP markers
MNP.markers <- FeaturePlot(lung.harmony, features = MNP.genes, pt.size = 0.5, ncol = 6) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"MNP.markers.png"), width=1800,height=900,units="px")
print(MNP.markers)
dev.off()

# feature plot with cDC 1 markers
cDC1.markers <- FeaturePlot(lung.harmony, features = cDC1.genes, pt.size = 0.5, ncol = 4) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"cDC1.markers.png"), width=1600,height=400,units="px")
print(cDC1.markers)
dev.off()

# feature plot with cDC 2 markers
cDC2.markers <- FeaturePlot(lung.harmony, features = cDC2.genes, pt.size = 0.5, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"cDC2.markers.png"), width=1200,height=800,units="px")
print(cDC2.markers)
dev.off()

# feature plot with pDC markers
pDC.markers <- FeaturePlot(lung.harmony, features = pDC.genes, pt.size = 0.5, ncol = 4) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"pDC.markers.png"), width=1600,height=400,units="px")
print(pDC.markers)
dev.off()

# feature plot with T cell markers
Tcell.markers <- FeaturePlot(lung.harmony, features = Tcell.genes, pt.size = 0.5, ncol = 4) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"Tcell.markers.png"), width=1600,height=800,units="px")
print(Tcell.markers)
dev.off()

# feature plot with NK cell markers
NK.markers <- FeaturePlot(lung.harmony, features = NK.genes, pt.size = 0.5, ncol = 4) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"NK.markers.png"), width=1600,height=800,units="px")
print(NK.markers)
dev.off()

# feature plot with B cell markers
Bcell.markers <- FeaturePlot(lung.harmony, features = Bcell.genes, pt.size = 0.5, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"Bcell.markers.png"), width=1200,height=400,units="px")
print(Bcell.markers)
dev.off()

# feature plot with plasma cell markers
plasmacell.markers <- FeaturePlot(lung.harmony, features = plasmacell.genes, pt.size = 0.5, ncol = 5) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"plasmacell.markers.png"), width=1500,height=300,units="px")
print(plasmacell.markers)
dev.off()

# feature plot with mast cell markers
mastcell.markers <- FeaturePlot(lung.harmony, features = mastcell.genes, pt.size = 0.5, ncol = 4) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"mastcell.markers.png"), width=1600,height=400,units="px")
print(mastcell.markers)
dev.off()

# feature plot with epithelial cell markers
epithelial.markers <- FeaturePlot(lung.harmony, features = epithelial.genes, pt.size = 0.5, ncol = 5) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"epithelial.markers.png"), width=1500,height=400,units="px")
print(epithelial.markers)
dev.off()

# feature plot with endothelial cell markers
endothelial.markers <- FeaturePlot(lung.harmony, features = endothelial.genes, pt.size = 0.5, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"endothelial.markers.png"), width=1200,height=1200,units="px")
print(endothelial.markers)
dev.off()

# feature plot with mesenchymal cell markers
mesecnhymal.markers <- FeaturePlot(lung.harmony, features = mesenchymal.genes, pt.size = 0.5, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"mesenchymal.markers.png"), width=1200,height=1200,units="px")
print(mesenchymal.markers)
dev.off()

# feature plot with erythroid cell markers
ery.markers <- FeaturePlot(lung.harmony, features = ery.genes, pt.size = 0.5, ncol = 2) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"ery.markers.png"), width=800,height=400,units="px")
print(ery.markers)
dev.off()

# proliferating cells
proliferation.markers <- FeaturePlot(lung.harmony, features = proliferation.genes, pt.size = 0.5, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"proliferation.markers.png"), width=800,height=400,units="px")
print(proliferation.markers)
dev.off()