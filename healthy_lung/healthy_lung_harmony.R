library(Seurat)
library(harmony)
library(future)
library(reticulate)
library(umap)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(patchwork)

### path variables ###
raredon.lung.path <- "/home/s1987963/ds_group/Niklas/raredon_lung/raredon_lung.rds"
reyfman.lung.path <- "/home/s1987963/ds_group/Niklas/reyfman_lung/reyfman_lung_healthy.rds"
habermann.lung.path <- "/home/s1987963/ds_group/Niklas/habermann_lung/habermann_lung_healthy.rds"
lung.path <- "/home/s1987963/ds_group/Niklas/healthy_lung/"
harmony.samples.path <- "/home/s1987963/ds_group/Niklas/healthy_lung/harmonize_samples/"

### read files ###
raredon.lung <- readRDS(raredon.lung.path) # 17867 cells
reyfman.lung <- readRDS(reyfman.lung.path) # 42144 cells
habermann.lung <- readRDS(habermann.lung.path) # 56222 cells

### add meta data ###
raredon.lung$study <- "raredon_lung"
raredon.lung$cohort <- "Chicago"
reyfman.lung$study <- "reyfman_lung"
reyfman.lung$cohort <- "Yale"
habermann.lung$study <- "habermann_lung"
habermann.lung$cohort <-"Nashville"

### normalization ###
raredon.lung <- NormalizeData(raredon.lung, normalization.method = "LogNormalize", scale.factor = 10000)
reyfman.lung <- NormalizeData(reyfman.lung, normalization.method = "LogNormalize", scale.factor = 10000)
habermann.lung <- NormalizeData(habermann.lung, normalization.method = "LogNormalize", scale.factor = 10000)

### merge data ###
lung <- merge(raredon.lung, c(reyfman.lung, habermann.lung)) #total 116233 cells

### feature selection ###
# select the top 4000 HVGs of each study and only consider genes that are HVGs in all studies
raredon.lung <- FindVariableFeatures(raredon.lung, selection.method = "vst", nfeatures = 4000)
reyfman.lung <- FindVariableFeatures(reyfman.lung, selection.method = "vst", nfeatures = 4000)
habermann.lung <- FindVariableFeatures(habermann.lung, selection.method = "vst", nfeatures = 4000)
# option 1: only consider genes that are HVGs in at least studies
lung.HVGs <- c()
lung.HVGs <- Reduce(append, list(VariableFeatures(raredon.lung), VariableFeatures(reyfman.lung), VariableFeatures(habermann.lung)))
lung.HVGs <- unique(lung.HVGs[which(table(lung.HVGs) > 1)])
# option 2: only consider genes that are HVGs in all studies
#lung.HVGs <- Reduce(intersect, list(VariableFeatures(raredon.lung), VariableFeatures(reyfman.lung), VariableFeatures(habermann.lung)))
# set HVGs
VariableFeatures(lung) <- lung.HVGs
#set new path variable
# harmony.samples.path <- "/home/s1987963/ds_group/Niklas/healthy_lung/harmonize_samples/HVGs_3studies/"
harmony.samples.path <- "/home/s1987963/ds_group/Niklas/healthy_lung/harmonize_samples/HVGs_2studies/"

# remove Seurat Objects
rm(raredon.lung)
rm(reyfman.lung)
rm(habermann.lung)

### scale data ###
lung <- ScaleData(lung) # uncorrected
# lung <- ScaleData(lung, vars.to.regress = "nFeature_RNA") # adjusted for nFeatures
# lung <- ScaleData(lung, vars.to.regress = c("nCount_RNA", "percent.mt")) # adjusted for nCounts + percent.mt

### dimensionality reduction: PCA ###
lung <- RunPCA(lung, features = VariableFeatures(object = lung))

# elbow plot
pca.elbow.plot <- ElbowPlot(lung, ndims = 50, reduction = "pca")
png(paste0(harmony.samples.path,"pca.elbow.plot.png"), width=1000,height=600,units="px")
print(pca.elbow.plot)
dev.off()

# visualize PCs
PC1_2.samples.plot <- DimPlot(object = lung, reduction = "pca", pt.size = .1, group.by = "patient.ID")
png(paste0(harmony.samples.path,"PC1_2.samples.png"), width=1000,height=1000,units="px")
print(PC1_2.samples.plot)
dev.off()

PC1_2.studies.plot <- DimPlot(object = lung, reduction = "pca", pt.size = .1, group.by = "study")
png(paste0(harmony.samples.path,"PC1_2.studies.png"), width=1000,height=1000,units="px")
print(PC1_2.studies.plot)
dev.off()

# save genes making up the PCs  
sink(paste0(harmony.samples.path, "PC_genes.txt"))
print(lung[["pca"]], dims = 1:50, nfeatures = 20)
sink()

### integration with harmony ###
lung.harmony <- lung %>% RunHarmony("patient.ID", theta = 2, reduction.save = "harmony_theta2", plot_convergence = TRUE) # harmonize all 12 samples independently

# harmony elbow plot
harmony.elbow.plot <- ElbowPlot(lung.harmony, ndims = 50, reduction = "harmony_theta2")
png(paste0(harmony.samples.path,"harmony_theta2.elbow.plot.png"), width=1000,height=600,units="px")
print(harmony.elbow.plot)
dev.off()

# visualize harmony PCs
harmony.PC1_2.samples.plot <- DimPlot(object = lung.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "patient.ID")
png(paste0(harmony.samples.path,"harmony_theta2.PC1_2.samples.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.samples.plot)
dev.off()

harmony.PC1_2.studies.plot <- DimPlot(object = lung.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "study")
png(paste0(harmony.samples.path,"harmony_theta2.PC1_2.studies.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.studies.plot)
dev.off()

# save genes making up the PCs
sink(paste0(harmony.samples.path, "harmony_PC_genes.txt"))
print(lung.harmony[["harmony_theta2"]], dims = 1:50, nfeatures = 20)
sink()

### save data ###
saveRDS(lung.harmony, paste0(harmony.samples.path, "healthy_lung_harmony_samples.rds"))

### explore different numbers of harmony PCs ###
# MNP genes
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

# visualize batch effect and marker gene expression
dims <- c(50) # harmonize samples
for(d in dims){
  
  # create folder
  dir.create(paste0(harmony.samples.path, "dim", d, "_annotation"))
  dim.path <- paste0(harmony.samples.path, "dim", d, "_annotation/")
  
  # run UMAP
  lung.harmony <- RunUMAP(lung.harmony, reduction = "harmony_theta2", dims = 1:d, seed.use=1)
  
  ## batch plots ##
  # no.1 overview - by patients | by study
  sample.batch.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.01)
  study.batch.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = "study", pt.size = 0.1)
  batch1.plot <- sample.batch.plot + study.batch.plot
  png(paste0(dim.path, "UMAP_dim", d, ".batch.png"), width=1500, height=600, units="px")
  print(batch1.plot)
  dev.off()
  
  # no.2 - coloured by patients
  batch2.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".patient.batch1.png"), width=1600, height=1000, units="px")
  print(batch2.plot)
  dev.off()
  
  # no.3 - split by patients
  batch3.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = NULL, split.by = "patient.ID", pt.size = 0.01, ncol = 8)
  png(paste0(dim.path, "UMAP_dim", d, ".patient.batch2.png"), width=1600, height=1000, units="px")
  print(batch3.plot)
  dev.off()
  
  # no.4 - coloured by study
  batch4.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = "study", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".study.batch1.png"), width=1600, height=1000, units="px")
  print(batch4.plot)
  dev.off()
  
  # no.5 - split by study
  batch5.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = NULL, split.by = "study", pt.size = 0.01, ncol = 3)
  png(paste0(dim.path, "UMAP_dim", d, ".study.batch2.png"), width=1800, height=600, units="px")
  print(batch5.plot)
  dev.off()
  
  ## QC metrics at cluster level ##
  cluster.qc.heatmap <- FeaturePlot(lung.harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "cluster.qc.heatmap.png"), width=1500,height=500,units="px")
  print(cluster.qc.heatmap)
  dev.off()
  
  ## gene expression heatmaps ##
  # feature plot with immune cell marker
  immunecell.markers <- FeaturePlot(lung.harmony, features = c("PTPRC"), pt.size = 0.5, ncol = 1) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "immunecell.markers.png"), width=1000,height=1000,units="px")
  print(immunecell.markers)
  dev.off()
  
  # feature plot with MNP markers
  MNP.markers <- FeaturePlot(lung.harmony, features = MNP.genes, pt.size = 0.5, ncol = 6) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "MNP.markers.png"), width=1800,height=900,units="px")
  print(MNP.markers)
  dev.off()
  
  # feature plot with cDC 1 markers
  cDC1.markers <- FeaturePlot(lung.harmony, features = cDC1.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "cDC1.markers.png"), width=1600,height=400,units="px")
  print(cDC1.markers)
  dev.off()
  
  # feature plot with cDC 2 markers
  cDC2.markers <- FeaturePlot(lung.harmony, features = cDC2.genes, pt.size = 0.5, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path, "UMAP_dim", d, "cDC2.markers.png"), width=1200,height=800,units="px")
  print(cDC2.markers)
  dev.off()
  
  # feature plot with pDC markers
  pDC.markers <- FeaturePlot(lung.harmony, features = pDC.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "pDC.markers.png"), width=1600,height=400,units="px")
  print(pDC.markers)
  dev.off()
  
  # feature plot with T cell markers
  Tcell.markers <- FeaturePlot(lung.harmony, features = Tcell.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "Tcell.markers.png"), width=1600,height=800,units="px")
  print(Tcell.markers)
  dev.off()
  
  # feature plot with NK cell markers
  NK.markers <- FeaturePlot(lung.harmony, features = NK.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "NK.markers.png"), width=1600,height=800,units="px")
  print(NK.markers)
  dev.off()
  
  # feature plot with B cell markers
  Bcell.markers <- FeaturePlot(lung.harmony, features = Bcell.genes, pt.size = 0.5, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "Bcell.markers.png"), width=1200,height=400,units="px")
  print(Bcell.markers)
  dev.off()
  
  # feature plot with plasma cell markers
  plasmacell.markers <- FeaturePlot(lung.harmony, features = plasmacell.genes, pt.size = 0.5, ncol = 5) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "plasmacell.markers.png"), width=1500,height=300,units="px")
  print(plasmacell.markers)
  dev.off()
  
  # feature plot with mast cell markers
  mastcell.markers <- FeaturePlot(lung.harmony, features = mastcell.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"mastcell.markers.png"), width=1600,height=400,units="px")
  print(mastcell.markers)
  dev.off()
  
  # feature plot with epithelial cell markers
  epithelial.markers <- FeaturePlot(lung.harmony, features = epithelial.genes, pt.size = 0.5, ncol = 5) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"epithelial.markers.png"), width=1500,height=400,units="px")
  print(epithelial.markers)
  dev.off()
  
  # feature plot with endothelial cell markers
  endothelial.markers <- FeaturePlot(lung.harmony, features = endothelial.genes, pt.size = 0.5, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"endothelial.markers.png"), width=1200,height=1200,units="px")
  print(endothelial.markers)
  dev.off()
  
  # feature plot with mesenchymal cell markers
  mesenchymal.markers <- FeaturePlot(lung.harmony, features = mesenchymal.genes, pt.size = 0.5, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"mesenchymal.markers.png"), width=1200,height=1200,units="px")
  print(mesenchymal.markers)
  dev.off()
  
  # feature plot with erythroid cell markers
  ery.markers <- FeaturePlot(lung.harmony, features = ery.genes, pt.size = 0.5, ncol = 2) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"ery.markers.png"), width=800,height=400,units="px")
  print(ery.markers)
  dev.off()
  
  # proliferating cells
  proliferation.markers <- FeaturePlot(lung.harmony, features = proliferation.genes, pt.size = 0.5, ncol = 2) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"proliferation.markers.png"), width=800,height=400,units="px")
  print(proliferation.markers)
  dev.off()
  
}

### clustering ###
# only for selected number of dims (for reasons of computational efficiency!)
dims <- c(20,50)
res <- seq(0.1,1.5,0.1)
for(d in dims){
  
  # set path variable
  dim.path <- paste0(harmony.samples.path, "dim", d, "_annotation/")
  
  # visualization with UMAP
  lung.harmony <- RunUMAP(lung.harmony, reduction = "harmony_theta2", dims=1:d, seed.use=1)
  
  # plot batch effect
  batch.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
  
  # create kNN graph
  lung.harmony <- FindNeighbors(lung.harmony, reduction = "harmony_theta2", dims = 1:d)
  
  for (r in res) {
    
    lung.harmony <- FindClusters(lung.harmony, reduction = "harmony_theta2", resolution = r)
    umap.plot <- DimPlot(lung.harmony, reduction = "umap", label = F, pt.size = 0.1)
    
    # create eval plot
    eval.plot <- umap.plot + batch.plot
    png(paste0(dim.path, "UMAP_dim", d, "_res", r, ".png"), width=1500, height=600, units="px")
    print(eval.plot)
    dev.off()
    
    # clustering at patient level
    patient.clustering <- DimPlot(lung.harmony, group.by = "seurat_clusters", split.by = "patient.ID", pt.size = 0.1, ncol = 8)
    png(paste0(harmony.path,"UMAP_dim", d, "_res", r,"patient.clustering.png"), width=1600,height=1000,units="px")
    print(patient.clustering)
    dev.off()
    
    # clustering at study level
    patient.clustering <- DimPlot(lung.harmony, group.by = "seurat_clusters", split.by = "study", pt.size = 0.01, ncol = 2)
    png(paste0(harmony.path,"UMAP_dim", d, "_res", r, "study.clustering.png"), width=1800,height=800,units="px")
    print(study.clustering)
    dev.off()
    
  }
}

