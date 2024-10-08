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
harmony.path <- "/home/s1987963/ds_group/Niklas/macparland_liver/harmony/dim17_annotation/"

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
# dims <- c() # uncorrected/adjusted
dims <- c(9,12,14,17)
res <- seq(0.5,1.5,0.1)
for(d in dims){
  liver.harmony <- RunUMAP(liver.harmony, reduction = "harmony", dims=1:d, seed.use=1)
  liver.harmony <- FindNeighbors(liver.harmony, reduction = "harmony", dims = 1:d)
  
  # plot batch effect
  batch.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
  
  for (r in res) {
    # perform clustering
    liver.harmony <- FindClusters(liver.harmony, reduction = "harmony", resolution = r)
    
    # plot clusters
    umap.plot <- DimPlot(liver.harmony, reduction = "umap", label = F, pt.size = 0.1)
    
    # plot clusters and batch effect side by side
    eval.plot <- umap.plot + batch.plot
    
    # save output
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
liver.harmony <- FindClusters(liver.harmony, reduction = "harmony", resolution = 1.0)
# run UMAP
liver.harmony <- RunUMAP(liver.harmony, reduction = "harmony", dims=1:17, seed.use=1)

### save data ###
saveRDS(liver.harmony, file = paste0(harmony.path, "macparland_liver_harmony.rds"))

## explore clustering at patient level ##
patient.clustering <- DimPlot(liver.harmony, group.by = "seurat_clusters", split.by = "patient.ID", pt.size = 0.2, ncol = 5)
png(paste0(harmony.path,"patient.clustering.png"), width=1500,height=400,units="px")
print(patient.clustering)
dev.off()

## QC metrics at cluster level ##
cluster.qc.heatmap <- FeaturePlot(liver.harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
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

# feature plot with immune cell marker
immunecell.markers <- FeaturePlot(liver.harmony, features = c("PTPRC"), pt.size = 0.5, ncol = 1) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"immunecell.markers.png"), width=1000,height=1000,units="px")
print(immunecell.markers)
dev.off()

# feature plot with MNP markers
MNP.markers <- FeaturePlot(liver.harmony, features = MNP.genes, pt.size = 0.5, ncol = 6) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"MNP.markers.png"), width=1800,height=900,units="px")
print(MNP.markers)
dev.off()

# feature plot with cDC 1 markers
cDC1.markers <- FeaturePlot(liver.harmony, features = cDC1.genes, pt.size = 0.5, ncol = 4) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"cDC1.markers.png"), width=1600,height=400,units="px")
print(cDC1.markers)
dev.off()

# feature plot with cDC 2 markers
cDC2.markers <- FeaturePlot(liver.harmony, features = cDC2.genes, pt.size = 0.5, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"cDC2.markers.png"), width=1200,height=800,units="px")
print(cDC2.markers)
dev.off()

# feature plot with pDC markers
pDC.markers <- FeaturePlot(liver.harmony, features = pDC.genes, pt.size = 0.5, ncol = 4) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"pDC.markers.png"), width=1600,height=400,units="px")
print(pDC.markers)
dev.off()

# feature plot with T cell markers
Tcell.markers <- FeaturePlot(liver.harmony, features = Tcell.genes, pt.size = 0.5, ncol = 4) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"Tcell.markers.png"), width=1600,height=800,units="px")
print(Tcell.markers)
dev.off()

# feature plot with NK cell markers
NK.markers <- FeaturePlot(liver.harmony, features = NK.genes, pt.size = 0.5, ncol = 4) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"NK.markers.png"), width=1600,height=800,units="px")
print(NK.markers)
dev.off()

# feature plot with B cell markers
Bcell.markers <- FeaturePlot(liver.harmony, features = Bcell.genes, pt.size = 0.5, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"Bcell.markers.png"), width=1200,height=400,units="px")
print(Bcell.markers)
dev.off()

# feature plot with plasma cell markers
plasmacell.markers <- FeaturePlot(liver.harmony, features = plasmacell.genes, pt.size = 0.5, ncol = 5) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"plasmacell.markers.png"), width=1500,height=300,units="px")
print(plasmacell.markers)
dev.off()

# feature plot with mast cell markers
mastcell.markers <- FeaturePlot(liver.harmony, features = mastcell.genes, pt.size = 0.5, ncol = 4) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"mastcell.markers.png"), width=1600,height=400,units="px")
print(mastcell.markers)
dev.off()

# feature plot with epithelial cell markers
epithelial.markers <- FeaturePlot(liver.harmony, features = epithelial.genes, pt.size = 0.5, ncol = 5) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"epithelial.markers.png"), width=1500,height=400,units="px")
print(epithelial.markers)
dev.off()

# feature plot with endothelial cell markers
endothelial.markers <- FeaturePlot(liver.harmony, features = endothelial.genes, pt.size = 0.5, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"endothelial.markers.png"), width=1200,height=1200,units="px")
print(endothelial.markers)
dev.off()

# feature plot with mesenchymal cell markers
mesenchymal.markers <- FeaturePlot(liver.harmony, features = mesenchymal.genes, pt.size = 0.5, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"mesenchymal.markers.png"), width=1200,height=1200,units="px")
print(mesenchymal.markers)
dev.off()

# feature plot with erythroid cell markers
ery.markers <- FeaturePlot(liver.harmony, features = ery.genes, pt.size = 0.5, ncol = 2) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"ery.markers.png"), width=800,height=400,units="px")
print(ery.markers)
dev.off()

# proliferating cells
proliferation.markers <- FeaturePlot(liver.harmony, features = proliferation.genes, pt.size = 0.5, ncol = 2) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(harmony.path,"proliferation.markers.png"), width=800,height=400,units="px")
print(proliferation.markers)
dev.off()

## compare PTPRC expression across clusters ##
umap.plot <- DimPlot(liver.harmony, reduction = "umap", label = T, label.size = 6, pt.size = 0.1)
ptprc.plot <- VlnPlot(object = liver.harmony, features = c("PTPRC"), group.by = "seurat_clusters", pt.size = 0.1) + NoLegend()
immuneclusters.plot <- umap.plot + immunecell.markers - ptprc.plot + plot_layout(ncol=1, widths=c(2,1))
png(paste0(harmony.path, "immuneclusters.png"), width=1800,height=1200,units="px")
print(immuneclusters.plot)
dev.off()