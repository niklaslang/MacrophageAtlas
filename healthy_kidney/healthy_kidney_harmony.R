library(Seurat)
library(harmony)
library(future)
library(reticulate)
library(umap)
library(dplyr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(patchwork)

### path variables ###
menon.kidney.path <- "/home/s1987963/ds_group/Niklas/menon_kidney/menon_kidney_filtered.rds"
liao.kidney.path <- "/home/s1987963/ds_group/Niklas/liao_kidney/liao_kidney_filtered.rds"
wilson.kidney.path <- "/home/s1987963/ds_group/Niklas/wilson_kidney/healthy/wilson_kidney_healthy_filtered.rds"
kidney.path <- "/home/s1987963/ds_group/Niklas/healthy_kidney/"
harmony.samples.path <- "/home/s1987963/ds_group/Niklas/healthy_kidney/harmonize_samples/"

### read files ###
menon.kidney <- readRDS(menon.kidney.path) # 25724 cells
liao.kidney <- readRDS(liao.kidney.path) #  22128 cells
wilson.kidney <- readRDS(wilson.kidney.path) #  13207 cells

### add meta data ###
menon.kidney$technology <- "scRNAseq"
liao.kidney$technology <- "scRNAseq"
wilson.kidney$technology <- "snRNAseq"

### normalization ###
menon.kidney <- NormalizeData(menon.kidney, normalization.method = "LogNormalize", scale.factor = 10000)
liao.kidney <- NormalizeData(liao.kidney, normalization.method = "LogNormalize", scale.factor = 10000)
wilson.kidney <- NormalizeData(wilson.kidney, normalization.method = "LogNormalize", scale.factor = 10000)

### merge data ###
kidney <- merge(menon.kidney, c(liao.kidney, wilson.kidney)) #total 61059 healthy kidney cells

### feature selection ###
# select the top 4000 HVGs of each study and only consider genes that are HVGs in all studies
menon.kidney <- FindVariableFeatures(menon.kidney, selection.method = "vst", nfeatures = 4000)
liao.kidney <- FindVariableFeatures(liao.kidney, selection.method = "vst", nfeatures = 4000)
wilson.kidney <- FindVariableFeatures(wilson.kidney, selection.method = "vst", nfeatures = 4000)
## option 1: only consider genes that are HVGs in at least studies
#kidney.HVGs <- c()
#kidney.HVGs <- Reduce(append, list(VariableFeatures(menon.kidney), VariableFeatures(liao.kidney), VariableFeatures(wilson.kidney)))
#kidney.HVGs <- unique(kidney.HVGs[which(table(kidney.HVGs) > 1)]) # 2672 HVGs
# option 2: only consider genes that are HVGs in all studies
kidney.HVGs <- Reduce(intersect, list(VariableFeatures(menon.kidney), VariableFeatures(liao.kidney), VariableFeatures(wilson.kidney))) # 784 HVGs
# set HVGs
VariableFeatures(kidney) <- kidney.HVGs
#set new path variable
harmony.samples.path <- "/home/s1987963/ds_group/Niklas/healthy_kidney/harmonize_samples/HVGs_3studies/"
#harmony.samples.path <- "/home/s1987963/ds_group/Niklas/healthy_kidney/harmonize_samples/HVGs_2studies/"

# remove Seurat Objects
rm(menon.kidney)
rm(liao.kidney)
rm(wilson.kidney)

### scale data ###
kidney <- ScaleData(kidney) # uncorrected
# kidney <- ScaleData(kidney, vars.to.regress = "nFeature_RNA") # adjusted for nFeatures
# kidney <- ScaleData(kidney, vars.to.regress = c("nCount_RNA", "percent.mt")) # adjusted for nCounts + percent.mt

### dimensionality reduction: PCA ###
kidney <- RunPCA(kidney, features = VariableFeatures(object = kidney))

# elbow plot
pca.elbow.plot <- ElbowPlot(kidney, ndims = 50, reduction = "pca")
png(paste0(harmony.samples.path,"pca.elbow.plot.png"), width=1000,height=600,units="px")
print(pca.elbow.plot)
dev.off()

# visualize PCs
PC1_2.samples.plot <- DimPlot(object = kidney, reduction = "pca", pt.size = .1, group.by = "patient.ID")
png(paste0(harmony.samples.path,"PC1_2.samples.png"), width=1000,height=1000,units="px")
print(PC1_2.samples.plot)
dev.off()

PC1_2.studies.plot <- DimPlot(object = kidney, reduction = "pca", pt.size = .1, group.by = "study")
png(paste0(harmony.samples.path,"PC1_2.studies.png"), width=1000,height=1000,units="px")
print(PC1_2.studies.plot)
dev.off()

# save genes making up the PCs  
sink(paste0(harmony.samples.path, "PC_genes.txt"))
print(kidney[["pca"]], dims = 1:50, nfeatures = 20)
sink()

### integration with harmony ###
kidney.harmony <- kidney %>% RunHarmony("patient.ID", theta = 2, reduction.save = "harmony_theta2", plot_convergence = TRUE) # harmonize all 30 samples independently

# harmony elbow plot
harmony.elbow.plot <- ElbowPlot(kidney.harmony, ndims = 50, reduction = "harmony_theta2")
png(paste0(harmony.samples.path,"harmony_theta2.elbow.plot.png"), width=1000,height=600,units="px")
print(harmony.elbow.plot)
dev.off()

# visualize harmony PCs
harmony.PC1_2.samples.plot <- DimPlot(object = kidney.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "patient.ID")
png(paste0(harmony.samples.path,"harmony_theta2.PC1_2.samples.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.samples.plot)
dev.off()

harmony.PC1_2.studies.plot <- DimPlot(object = kidney.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "study")
png(paste0(harmony.samples.path,"harmony_theta2.PC1_2.studies.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.studies.plot)
dev.off()

# save genes making up the PCs
sink(paste0(harmony.samples.path, "harmony_PC_genes.txt"))
print(kidney.harmony[["harmony_theta2"]], dims = 1:50, nfeatures = 20)
sink()

### save data ###
saveRDS(kidney.harmony, paste0(harmony.samples.path, "healthy_kidney_harmony_samples.rds"))

# visualize batch effect
dims <- c(14,30,36,40,50)
for(d in dims){
  
  # create folder
  dir.create(paste0(harmony.samples.path, "dim", d, "_annotation"))
  dim.path <- paste0(harmony.samples.path, "dim", d, "_annotation/")
  
  # run UMAP
  kidney.harmony <- RunUMAP(kidney.harmony, reduction = "harmony_theta2", dims = 1:d, seed.use=1)
  
  ## batch plots ##
  # no.1 overview - by patients | by study
  sample.batch.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.01)
  study.batch.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = "study", pt.size = 0.1)
  batch1.plot <- sample.batch.plot + study.batch.plot
  png(paste0(dim.path, "UMAP_dim", d, ".batch.png"), width=1500, height=600, units="px")
  print(batch1.plot)
  dev.off()
  
  # no.2 - colored by patients
  batch2.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".patient.batch1.png"), width=1600, height=1000, units="px")
  print(batch2.plot)
  dev.off()
  
  # no.3 - split by patients
  batch3.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = "study", split.by = "patient.ID", pt.size = 0.01, ncol = 6)
  png(paste0(dim.path, "UMAP_dim", d, ".patient.batch2.png"), width=1800, height=1500, units="px")
  print(batch3.plot)
  dev.off()
  
  # no.4 - colored by study
  batch4.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = "study", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".study.batch1.png"), width=1600, height=1000, units="px")
  print(batch4.plot)
  dev.off()
  
  # no.5 - split by study
  batch5.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = "study", split.by = "study", pt.size = 0.01, ncol = 3)
  png(paste0(dim.path, "UMAP_dim", d, ".study.batch2.png"), width=1800, height=600, units="px")
  print(batch5.plot)
  dev.off()
}

### marker genes ###
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

# broad lineage markers
overview.markers <- c("PTPRC", "EPCAM", "PECAM1", "PDGFRA")

# epithelial cells
epithelial.genes <- c("EPCAM",
                      "NPHS1", "NPHS2", "WT1", #podocytes
                      "KIT", # CD IC A ("SCL26A4", # CD IC B)
                      "AQP2", "KRT18", #CD PC
                      "SLC8A1", #CNT
                      "SLC12A3", #DCT
                      "SLC12A1", #TAL
                      "SPP1", #tLOH
                      "LRP2", "SLC34A1" #PCT
                      )

# endothelial cells
endothelial.genes <- c("KDR",	"CD34", "VWF", "CLDN5", "PECAM1", 
                       "ICAM2",	"PDPN",	"CLEC14A",	"MMRN1", "FLT1")

# mesenchymal cells
mesenchymal.genes <- c("COL1A1",	"COL3A1", "ACTA2", "MYH11",	"PDGFRA", 
                       "PDGFRB", "MSLN",	"DCN",
                       "RGS5" #vSMC
                       )

# proliferating cells
proliferation.genes <- c("MKI67",	"TOP2A")

# visualize marker gene expression
dims <- c(40)
for(d in dims){
  
  # create folder
  #dir.create(paste0(harmony.samples.path, "dim", d, "_annotation"))
  dim.path <- paste0(harmony.samples.path, "dim", d, "_annotation/")
  
  # run UMAP
  kidney.harmony <- RunUMAP(kidney.harmony, reduction = "harmony_theta2", dims = 1:d, seed.use=1)
  
  ### batch plots ##
  ## no.1 overview - by patients | by study
  #sample.batch.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.01)
  #study.batch.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = "study", pt.size = 0.1)
  #batch1.plot <- sample.batch.plot + study.batch.plot
  #png(paste0(dim.path, "UMAP_dim", d, ".batch.png"), width=1500, height=600, units="px")
  #print(batch1.plot)
  #dev.off()
  #
  ## no.2 - coloured by patients
  #batch2.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.01)
  #png(paste0(dim.path, "UMAP_dim", d, ".patient.batch1.png"), width=1600, height=1000, units="px")
  #print(batch2.plot)
  #dev.off()
  #
  ## no.3 - split by patients
  #batch3.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = NULL, split.by = "patient.ID", pt.size = 0.01, ncol = 8)
  #png(paste0(dim.path, "UMAP_dim", d, ".patient.batch2.png"), width=1600, height=1000, units="px")
  #print(batch3.plot)
  #dev.off()
  #
  ## no.4 - coloured by study
  #batch4.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = "study", pt.size = 0.01)
  #png(paste0(dim.path, "UMAP_dim", d, ".study.batch1.png"), width=1600, height=1000, units="px")
  #print(batch4.plot)
  #dev.off()
  #
  ## no.5 - split by study
  #batch5.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = NULL, split.by = "study", pt.size = 0.01, ncol = 3)
  #png(paste0(dim.path, "UMAP_dim", d, ".study.batch2.png"), width=1800, height=600, units="px")
  #print(batch5.plot)
  #dev.off()
  
  ## QC metrics at cluster level ##
  cluster.qc.heatmap <- FeaturePlot(kidney.harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "cluster.qc.heatmap.png"), width=1500,height=500,units="px")
  print(cluster.qc.heatmap)
  dev.off()
  
  ## gene expression heatmaps ##
  # overview feature plot
  overview.markers <- FeaturePlot(kidney.harmony, features = overview.markers, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "overview.markers.png"), width=1600,height=500,units="px")
  print(overview.markers)
  dev.off()
  
  # feature plot with immune cell marker
  immunecell.markers <- FeaturePlot(kidney.harmony, features = c("PTPRC"), pt.size = 0.5, ncol = 1) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "immunecell.markers.png"), width=1000,height=1000,units="px")
  print(immunecell.markers)
  dev.off()
  
  # feature plot with MNP markers
  MNP.markers <- FeaturePlot(kidney.harmony, features = MNP.genes, pt.size = 0.5, ncol = 6) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "MNP.markers.png"), width=1800,height=900,units="px")
  print(MNP.markers)
  dev.off()
  
  # feature plot with cDC 1 markers
  cDC1.markers <- FeaturePlot(kidney.harmony, features = cDC1.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "cDC1.markers.png"), width=1600,height=400,units="px")
  print(cDC1.markers)
  dev.off()
  
  # feature plot with cDC 2 markers
  cDC2.markers <- FeaturePlot(kidney.harmony, features = cDC2.genes, pt.size = 0.5, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path, "UMAP_dim", d, "cDC2.markers.png"), width=1200,height=800,units="px")
  print(cDC2.markers)
  dev.off()
  
  # feature plot with pDC markers
  pDC.markers <- FeaturePlot(kidney.harmony, features = pDC.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "pDC.markers.png"), width=1600,height=400,units="px")
  print(pDC.markers)
  dev.off()
  
  # feature plot with T cell markers
  Tcell.markers <- FeaturePlot(kidney.harmony, features = Tcell.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "Tcell.markers.png"), width=1600,height=800,units="px")
  print(Tcell.markers)
  dev.off()
  
  # feature plot with NK cell markers
  NK.markers <- FeaturePlot(kidney.harmony, features = NK.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "NK.markers.png"), width=1600,height=800,units="px")
  print(NK.markers)
  dev.off()
  
  # feature plot with B cell markers
  Bcell.markers <- FeaturePlot(kidney.harmony, features = Bcell.genes, pt.size = 0.5, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "Bcell.markers.png"), width=1200,height=400,units="px")
  print(Bcell.markers)
  dev.off()
  
  # feature plot with plasma cell markers
  plasmacell.markers <- FeaturePlot(kidney.harmony, features = plasmacell.genes, pt.size = 0.5, ncol = 5) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "plasmacell.markers.png"), width=1500,height=300,units="px")
  print(plasmacell.markers)
  dev.off()
  
  # feature plot with mast cell markers
  mastcell.markers <- FeaturePlot(kidney.harmony, features = mastcell.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"mastcell.markers.png"), width=1600,height=400,units="px")
  print(mastcell.markers)
  dev.off()
  
  # feature plot with epithelial cell markers
  epithelial.markers <- FeaturePlot(kidney.harmony, features = epithelial.genes, pt.size = 0.5, ncol = 5) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"epithelial.markers.png"), width=1500,height=900,units="px")
  print(epithelial.markers)
  dev.off()
  
  # feature plot with endothelial cell markers
  endothelial.markers <- FeaturePlot(kidney.harmony, features = endothelial.genes, pt.size = 0.5, ncol = 5) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"endothelial.markers.png"), width=1500,height=600,units="px")
  print(endothelial.markers)
  dev.off()
  
  # feature plot with mesenchymal cell markers
  mesenchymal.markers <- FeaturePlot(kidney.harmony, features = mesenchymal.genes, pt.size = 0.5, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"mesenchymal.markers.png"), width=1200,height=1200,units="px")
  print(mesenchymal.markers)
  dev.off()
  
  # feature plot with erythroid cell markers
  ery.markers <- FeaturePlot(kidney.harmony, features = ery.genes, pt.size = 0.5, ncol = 2) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"ery.markers.png"), width=800,height=400,units="px")
  print(ery.markers)
  dev.off()
  
  # proliferating cells
  proliferation.markers <- FeaturePlot(kidney.harmony, features = proliferation.genes, pt.size = 0.5, ncol = 2) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"proliferation.markers.png"), width=800,height=400,units="px")
  print(proliferation.markers)
  dev.off()
  
}


### clustering ###
# only for selected number of dims (for reasons of computational efficiency!)
dims <- c(36,40,50)
res <- seq(1.5,2.0,0.1)
for(d in dims){
  
  # set path variable
  dim.path <- paste0(harmony.samples.path, "dim", d, "_annotation/")
  
  # visualization with UMAP
  kidney.harmony <- RunUMAP(kidney.harmony, reduction = "harmony_theta2", dims=1:d, seed.use=1)
  
  # plot batch effect
  batch.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = "study", pt.size = 0.1)
  
  # create kNN graph
  kidney.harmony <- FindNeighbors(kidney.harmony, reduction = "harmony_theta2", dims = 1:d)
  
  for (r in res) {
    
    kidney.harmony <- FindClusters(kidney.harmony, reduction = "harmony_theta2", resolution = r)
    umap.plot <- DimPlot(kidney.harmony, reduction = "umap", label = F, pt.size = 0.1)
    
    # create eval plot
    eval.plot <- umap.plot + batch.plot
    png(paste0(dim.path, "UMAP_dim", d, "_res", r, ".png"), width=1500, height=600, units="px")
    print(eval.plot)
    dev.off()
    
    # clustering at patient level
    patient.clustering <- DimPlot(kidney.harmony, group.by = "seurat_clusters", split.by = "patient.ID", pt.size = 0.1, ncol = 6)
    png(paste0(dim.path,"UMAP_dim", d, "_res", r,".patient.clustering.png"), width=1800,height=1500,units="px")
    print(patient.clustering)
    dev.off()
    
    # clustering at study level
    study.clustering <- DimPlot(kidney.harmony, group.by = "seurat_clusters", split.by = "study", pt.size = 0.01, ncol = 3)
    png(paste0(dim.path,"UMAP_dim", d, "_res", r, ".study.clustering.png"), width=1800,height=600,units="px")
    print(study.clustering)
    dev.off()
    
  }
}
