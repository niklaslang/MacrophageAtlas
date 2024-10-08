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
raredon.lung.path <- "/home/s1987963/ds_group/Niklas/raredon_lung/raredon_lung.rds"
reyfman.lung.path <- "/home/s1987963/ds_group/Niklas/reyfman_lung/reyfman_lung_healthy.rds"
habermann.lung.path <- "/home/s1987963/ds_group/Niklas/habermann_lung/healthy/habermann_lung_healthy.rds"
lung.path <- "/home/s1987963/ds_group/Niklas/healthy_lung/"
harmony.samples.path <- "/home/s1987963/ds_group/Niklas/healthy_lung/harmonize_samples/"

### read files ###
raredon.lung <- readRDS(raredon.lung.path) # 17867 cells
reyfman.lung <- readRDS(reyfman.lung.path) # 42144 cells
habermann.lung <- readRDS(habermann.lung.path) # 44873  cells

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
lung <- merge(raredon.lung, c(reyfman.lung, habermann.lung)) #total 104884 cells

### feature selection ###
# select the top 4000 HVGs of each study and only consider genes that are HVGs in all studies
raredon.lung <- FindVariableFeatures(raredon.lung, selection.method = "vst", nfeatures = 4000)
reyfman.lung <- FindVariableFeatures(reyfman.lung, selection.method = "vst", nfeatures = 4000)
habermann.lung <- FindVariableFeatures(habermann.lung, selection.method = "vst", nfeatures = 4000)
## option 1: only consider genes that are HVGs in at least studies
#lung.HVGs <- c()
#lung.HVGs <- Reduce(append, list(VariableFeatures(raredon.lung), VariableFeatures(reyfman.lung), VariableFeatures(habermann.lung)))
#lung.HVGs <- unique(lung.HVGs[which(table(lung.HVGs) > 1)]) # 2745 HVGs
# option 2: only consider genes that are HVGs in all studies
lung.HVGs <- Reduce(intersect, list(VariableFeatures(raredon.lung), VariableFeatures(reyfman.lung), VariableFeatures(habermann.lung))) # 1585 HVGs
# set HVGs
VariableFeatures(lung) <- lung.HVGs
#set new path variable
harmony.samples.path <- "/home/s1987963/ds_group/Niklas/healthy_lung/harmonize_samples/HVGs_3studies/"
#harmony.samples.path <- "/home/s1987963/ds_group/Niklas/healthy_lung/harmonize_samples/HVGs_2studies/"

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
lung.harmony <- lung %>% RunHarmony("patient.ID", theta = 2, reduction.save = "harmony_theta2", plot_convergence = TRUE) # harmonize all 32 samples independently

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

# broad lineage markers
lineage.genes <- c("PTPRC", "EPCAM", "PECAM1", "PDGFRA")

# epithelial cells
epithelial.genes <- c("EPCAM",
                      "VEGFA", "AGER", # AT I
                      "SFTPC", # AT II
                      "SCGB1A1", "MUC5AC", # club cells (clara cells)
                      "FOXJ1", # ciliated
                      "PIFO", "KCTD12", # ciliated 1
                      "APOD", # ciliated 2
                      "KRT5", "TP63", # basal stem cells
                      "MUC5B", "CSF3", "SCGB3A2", # goblets (secretory cells)
                      "CALCA", # neuroendocrine
                      "FOXI1", "CFTR" # ionocytes
                      )

# endothelial cells
endothelial.genes <- c("KDR",	"CD34", "VWF", "CLDN5", "PECAM1", 
                       "ICAM2",	"PDPN",	"CLEC14A",	"MMRN1")

# mesenchymal cells
mesenchymal.genes <- c("COL1A1",	"COL3A1", "ACTA2", "MYH11",	"PDGFRA", 
                       "PDGFRB", "RGS5",	"MSLN",	"DCN", "WT1")

# proliferating cells
proliferation.genes <- c("MKI67",	"TOP2A")

# visualize batch effect
dims <- c(50,40,30,24,14) # harmonize samples
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
  
  # no.2 - colored by patients
  batch2.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".patient.batch1.png"), width=1600, height=1000, units="px")
  print(batch2.plot)
  dev.off()
  
  # no.3 - split by patients
  batch3.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = "study", split.by = "patient.ID", pt.size = 0.01, ncol = 8)
  png(paste0(dim.path, "UMAP_dim", d, ".patient.batch2.png"), width=1600, height=1000, units="px")
  print(batch3.plot)
  dev.off()
  
  # no.4 - colored by study
  batch4.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = "study", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".study.batch1.png"), width=1600, height=1000, units="px")
  print(batch4.plot)
  dev.off()
  
  # no.5 - split by study
  batch5.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = "study", split.by = "study", pt.size = 0.01, ncol = 3)
  png(paste0(dim.path, "UMAP_dim", d, ".study.batch2.png"), width=1800, height=600, units="px")
  print(batch5.plot)
  dev.off()
}

# visualize marker gene expression
dims <- c(24)
for(d in dims){
  
  # create folder
  #dir.create(paste0(harmony.samples.path, "dim", d, "_annotation"))
  dim.path <- paste0(harmony.samples.path, "dim", d, "_annotation/")
  
  # run UMAP
  lung.harmony <- RunUMAP(lung.harmony, reduction = "harmony_theta2", dims = 1:d, seed.use=1)

  ## QC metrics at cluster level ##
  cluster.qc.heatmap <- FeaturePlot(lung.harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "cluster.qc.heatmap.png"), width=1500,height=500,units="px")
  print(cluster.qc.heatmap)
  dev.off()
  
  ## gene expression heatmaps ##
  # overview feature plot
  lineage.overview <- FeaturePlot(lung.harmony, features = lineage.genes, pt.size = 0.5, ncol = 4) & 
   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "overview.markers.png"), width=1600,height=500,units="px")
  print(lineage.overview )
  dev.off()
  
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
  epithelial.markers <- FeaturePlot(lung.harmony, features = epithelial.genes, pt.size = 0.5, ncol = 6) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"epithelial.markers.png"), width=1800,height=1000,units="px")
  print(epithelial.markers)
  dev.off()
  
  # feature plot with endothelial cell markers
  endothelial.markers <- FeaturePlot(lung.harmony, features = endothelial.genes, pt.size = 0.5, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"endothelial.markers.png"), width=1200,height=1200,units="px")
  print(endothelial.markers)
  dev.off()
  
  # feature plot with mesenchymal cell markers
  mesenchymal.markers <- FeaturePlot(lung.harmony, features = mesenchymal.genes, pt.size = 0.5, ncol = 5) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"mesenchymal.markers.png"), width=1500,height=600,units="px")
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
  
  rm(cluster.qc.heatmap, lineage.overview, mesenchymal.markers, endothelial.markers, epithelial.markers, immunecell.markers,
     MNP.markers, cDC1.markers, cDC2.markers, pDC.markers, mastcell.markers, 
     Tcell.markers, NK.markers, Bcell.markers, plasmacell.markers,
     ery.markers, proliferation.markers)
}

### clustering ###
# only for selected number of dims (for reasons of computational efficiency!)
dims <- c(40)
res <- seq(1.1,2.0,0.1)
for(d in dims){
  
  # set path variable
  dim.path <- paste0(harmony.samples.path, "dim", d, "_annotation/")
  
  # visualization with UMAP
  lung.harmony <- RunUMAP(lung.harmony, reduction = "harmony_theta2", dims=1:d, seed.use=1)
  
  # plot batch effect
  batch.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = "study", pt.size = 0.1)
  
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
    png(paste0(dim.path,"UMAP_dim", d, "_res", r,".patient.clustering.png"), width=1600,height=1000,units="px")
    print(patient.clustering)
    dev.off()
    
    # clustering at study level
    study.clustering <- DimPlot(lung.harmony, group.by = "seurat_clusters", split.by = "study", pt.size = 0.01, ncol = 3)
    png(paste0(dim.path,"UMAP_dim", d, "_res", r, ".study.clustering.png"), width=1800,height=600,units="px")
    print(study.clustering)
    dev.off()
    
  }
}

### preliminary clustering ###
lung.harmony <- FindNeighbors(lung.harmony, reduction = "harmony_theta2", dims = 1:40)
lung.harmony <- FindClusters(lung.harmony, reduction = "harmony_theta2", resolution = 1.5)
# run UMAP
lung.harmony <- RunUMAP(lung.harmony, reduction = "harmony_theta2", dims = 1:40, seed.use=1)

# compare PTPRC expression across clusters
umap.plot <- DimPlot(lung.harmony, reduction = "umap", label = T, label.size = 6, pt.size = 0.1)
ptprc.plot <- VlnPlot(object = lung.harmony, features = c("PTPRC"), group.by = "seurat_clusters", pt.size = 0.1) + NoLegend()
immuneclusters.plot <- umap.plot + immunecell.markers - ptprc.plot + plot_layout(ncol=1, widths=c(2,1))
png(paste0(harmony.samples.path, "dim40_annotation/immuneclusters.png"), width=1800,height=1200,units="px")
print(immuneclusters.plot)
dev.off()

### compute cluster marker genes ###
lung.markers <- FindAllMarkers(lung.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
lung.top50.markers <- lung.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
write.csv(lung.markers, file = paste0(harmony.samples.path, "dim40_annotation/ALL_marker_genes.csv"))
write.csv(lung.top50.markers, file = paste0(harmony.samples.path, "dim40_annotation/top50_marker_genes.csv"))

### cell type annotation ###
cluster.annotation <- c("Healthy Lung Alveolar Macrophage 1", "Healthy Lung ATII 1", "Healthy Lung ATII 2", 
                        "Healthy Lung Macrophage 1", "Healthy Lung Macrophage 2", "Healthy Lung Monocyte", 
                        "Healthy Lung T cell", "Helathy Lung Club Cell", "Healthy Lung NK 1", "Healthy Lung Ciliated 1",
                        "Healthy Lung ATI", "Healthy Lung Cap Endo 1", "Healthy Lung cDC", "Healthy Lung Mesenchyme 2",
                        "Healthy Lung NK 2", "Healthy Lung Proliferating 1", "Healthy Lung Goblet", "Healthy Lung LymphEndo", 
                        "Healthy Lung ArtEndo", "Healthy Lung Venous/Peribronchial Endo", "Healthy Lung Cap Endo 2", 
                        "Healthy Lung Macrophage 3", "Healthy Lung Ciliated 2", "Healthy Lung Mast cell", "Healthy Lung B cell", 
                        "Healthy Lung Plasma cell", "Healthy Lung Mesenchyme 1", "Healthy Lung Macrophage 4", "Healthy Lung ATII 3", 
                        "Healthy Lung ATII 4", "Healthy Lung Macrophage 5", "NK/MP doublets", "Healthy Lung Macrophage 6", 
                        "Healthy Lung Basal Cell", "Healthy Lung ATII 5", "Healthy Lung Macrophage 7", "Endo/mac doublet",
                        "NK/Epithelial doublets", "Healthy Lung Alveolar Macrophage 2", "Mes/mac doublet", 
                        "Healthy Lung Alveolar Macrophage 3", "Endo/Ep Doublet", "B cell/Ep Doublet"
                        )

names(cluster.annotation) <- levels(lung.harmony)
lung.harmony <- RenameIdents(lung.harmony, cluster.annotation)

# save annotated UMAP
annotated.umap.plot <- DimPlot(lung.harmony, reduction = "umap", label = T, label.size = 5, pt.size = 0.1)
png(paste0(harmony.samples.path, "dim40_annotation/UMAP_annotated.png"), width=1800,height=1200,units="px")
print(annotated.umap.plot)
dev.off()

### cell lineage annotation ###
cell.data <- data.table(barcode = colnames(lung.harmony),
                        celltype = Idents(lung.harmony))

lineage.annotation <- c("MP", "Epithelia", "Epithelia", "MP","MP", "MP","T cell", "Epithelia",
                        "NK cell", "Epithelia", "Epithelia", "Endothelia","MP", "Mesenchyme","NK cell",
                        "Proliferating", "Epithelia", "Endothelia", "Endothelia", "Endothelia", "Endothelia", 
                        "MP", "Epithelia", "Mast cell", "B cell", "Plasma cell", "Mesenchyme", "MP", "Epithelia", 
                        "Epithelia", "MP", "Doublet", "MP", "Epithelia", "Epithelia", "MP", "Doublet",
                        "Doublet", "MP", "Doublet", "MP", "Doublet", "Doublet"
                        )

lineage.data <- data.table(celltype = cluster.annotation, lineage = lineage.annotation)
meta.data <- merge(cell.data, lineage.data, by = "celltype")
meta.data <- data.frame(meta.data, row.names = meta.data$barcode)
meta.data$barcode <- NULL
meta.data$celltype <- NULL
lung.harmony <- AddMetaData(lung.harmony, meta.data, col.name = "lineage")

# save annotated UMAP
annotated.umap.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = "lineage", label = T, label.size = 5, pt.size = 0.1)
png(paste0(harmony.samples.path, "dim40_annotation/UMAP_annotated.lineage.png"), width=1800,height=1200,units="px")
print(annotated.umap.plot)
dev.off()

### save R session ###
save.image(file = "/home/s1987963/ds_group/Niklas/healthy_lung/harmonize_samples/HVGs_3studies/healthy_lung_harmony.RData")

### save data ###
saveRDS(lung.harmony, "/home/s1987963/ds_group/Niklas/healthy_organs/healthy_lung_annotated.rds")
