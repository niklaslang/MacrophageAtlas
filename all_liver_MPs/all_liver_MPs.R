library(Seurat)
library(harmony)
library(future)
library(reticulate)
library(umap)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(patchwork)
library(data.table)

### path variables ###
MP.path <- "/home/s1987963/ds_group/Niklas/combined_MPs/combined_liver_MPs_annotated.rds"
harmony.samples.path <- "/home/s1987963/ds_group/Niklas/MP_liver/harmonize_samples/HVGs_conditions/"

### load data ###
MP <- readRDS(MP.path) # 16646 cells

#### split healthy data ###
MP.healthy <- subset(MP, subset = condition == "healthy")
MP.fibrotic <- subset(MP, subset = condition == "fibrotic")

### feature selection: HVGs ###
#MP <- FindVariableFeatures(MP, selection.method = "vst", nfeatures = 2000)
MP.healthy <- FindVariableFeatures(MP.healthy, selection.method = "vst", nfeatures = 4000)
MP.fibrotic <- FindVariableFeatures(MP.fibrotic, selection.method = "vst", nfeatures = 4000)
MP.HVGs <- intersect(VariableFeatures(MP.fibrotic), VariableFeatures(MP.healthy)) # 1631
VariableFeatures(MP) <- MP.HVGs

### scale data ###
MP <- ScaleData(MP) # uncorrected

### dimensionality reduction: PCA ###
MP <- RunPCA(MP, features = VariableFeatures(object = MP))

# elbow plot
pca.elbow.plot <- ElbowPlot(MP, ndims = 50, reduction = "pca")
png(paste0(harmony.samples.path,"pca.elbow.plot.png"), width=1000,height=600,units="px")
print(pca.elbow.plot)
dev.off()

# visualize PCs
PC1_2.samples.plot <- DimPlot(object = MP, reduction = "pca", pt.size = .1, group.by = "patient.ID")
png(paste0(harmony.samples.path,"PC1_2.samples.png"), width=1000,height=1000,units="px")
print(PC1_2.samples.plot)
dev.off()

PC1_2.studies.plot <- DimPlot(object = MP, reduction = "pca", pt.size = .1, group.by = "dataset")
png(paste0(harmony.samples.path,"PC1_2.studies.png"), width=1000,height=1000,units="px")
print(PC1_2.studies.plot)
dev.off()

# save genes making up the PCs  
sink(paste0(harmony.samples.path, "PC_genes.txt"))
print(MP[["pca"]], dims = 1:50, nfeatures = 20)
sink()

### integration with harmony ###
# harmonize samples
MP.harmony <- MP %>% RunHarmony("patient.ID", theta = 2, reduction.save = "harmony_theta2", plot_convergence = TRUE) # harmonize all 12 samples independently

## harmony elbow plot ##
harmony.elbow.plot <- ElbowPlot(MP.harmony, ndims = 50, reduction = "harmony_theta2")
png(paste0(harmony.samples.path,"harmony_theta2.elbow.plot.png"), width=1000,height=600,units="px")
print(harmony.elbow.plot)
dev.off()

## explore harmony coordinates ##
## visualize PCs ##
harmony.PC1_2.samples.plot <- DimPlot(object = MP.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "patient.ID")
png(paste0(harmony.samples.path,"harmony_theta2.PC1_2.samples.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.samples.plot)
dev.off()

harmony.PC1_2.studies.plot <- DimPlot(object = MP.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "dataset")
png(paste0(harmony.samples.path,"harmony_theta2.PC1_2.studies.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.studies.plot)
dev.off()

### save genes making up the PCs ### 
sink(paste0(harmony.samples.path, "harmony_PC_genes.txt"))
print(MP.harmony[["harmony_theta2"]], dims = 1:50, nfeatures = 20)
sink()

### save data ###
saveRDS(MP.harmony, paste0(harmony.samples.path, "all_MP_harmony.rds"))

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

# proliferating cells
proliferation.genes <- c("MKI67",	"TOP2A")

# visualize batch effect and marker gene expression
dims <- c(50)
for(d in dims){
  
  # create folder
  dir.create(paste0(harmony.samples.path, "dim", d, "_annotation"))
  dim.path <- paste0(harmony.samples.path, "dim", d, "_annotation/")
  
  # run UMAP
  MP.harmony <- RunUMAP(MP.harmony, reduction = "harmony_theta2", dims = 1:d, seed.use=1)
  
  ## visualize batch effect 
  # no.1 overview
  sample.batch.plot <- DimPlot(MP.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
  dataset.batch.plot <- DimPlot(MP.harmony, reduction = "umap", group.by = "dataset", pt.size = 0.1)
  condition.batch.plot <- DimPlot(MP.harmony, reduction = "umap", group.by = "condition", pt.size = 0.1)
  batch1.plot <- sample.batch.plot + dataset.batch.plot + condition.batch.plot
  png(paste0(dim.path, "UMAP_dim", d, ".batch.png"), width=1800, height=600, units="px")
  print(batch1.plot)
  dev.off()
  
  # no.2 - colored by patients
  batch2.plot <- DimPlot(MP.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".patient.batch1.png"), width=1600, height=1000, units="px")
  print(batch2.plot)
  dev.off()
  
  # no.3 - split by patients
  batch3.plot <- DimPlot(MP.harmony, reduction = "umap", group.by = "study", split.by = "patient.ID", pt.size = 0.01, ncol = 6)
  png(paste0(dim.path, "UMAP_dim", d, ".patient.batch2.png"), width=1800, height=600, units="px")
  print(batch3.plot)
  dev.off()
  
  # no.4 - colored by dataset
  batch4.plot <- DimPlot(MP.harmony, reduction = "umap", group.by = "dataset", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".study.batch1.png"), width=1600, height=1000, units="px")
  print(batch4.plot)
  dev.off()
  
  # no.5 - split by dataset
  batch5.plot <- DimPlot(MP.harmony, reduction = "umap", split.by = "dataset", pt.size = 0.01, ncol = 3)
  png(paste0(dim.path, "UMAP_dim", d, ".study.batch2.png"), width=1800, height=600, units="px")
  print(batch5.plot)
  dev.off()
  
  # no.6 - colored by condition
  batch6.plot <- DimPlot(MP.harmony, reduction = "umap", group.by = "condition", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".condition.batch1.png"), width=1000, height=1000, units="px")
  print(batch6.plot)
  dev.off()
  
  # no.7 - split by condition
  batch7.plot <- DimPlot(MP.harmony, reduction = "umap", group.by = "study", split.by = "condition", pt.size = 0.01, ncol = 2)
  png(paste0(dim.path, "UMAP_dim", d, ".condition.batch2.png"), width=1800, height=900, units="px")
  print(batch7.plot)
  dev.off()
  
  # no.9 - colored by lineage annotation
  batch9.plot <- DimPlot(MP.harmony, reduction = "umap", group.by = "lineage_integrated", label=TRUE, pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".lineage.batch2.png"), width=1600, height=1000, units="px")
  print(batch9.plot)
  dev.off()

  ## compare healthy vs. fibrotic ##
  compare <- DimPlot(MP.harmony, reduction = "umap", split.by = "condition", group.by = "condition", 
                     pt.size = 0.01) + NoLegend() + scale_color_viridis(discrete=TRUE)
  png(paste0(dim.path, "UMAP_dim", d, "healthy_fibrotic_1.png"), width=1800, height=900, units="px")
  print(compare)
  dev.off()
  
  ## compare original annotation ##
  annotation <- DimPlot(MP.harmony, reduction = "umap", split.by = "condition", group.by = "celltype_integrated",
                        label = TRUE, repel = TRUE, pt.size = 0.01) + NoLegend()
  png(paste0(dim.path, "UMAP_dim", d, "healthy_fibrotic_2.png"), width=1800, height=900, units="px")
  print(annotation)
  dev.off()
  
  ## compare original annotation in healthy vs. fibrotic data ##
  plot <- compare / annotation
  png(paste0("UMAP_dim", d, "healthy_fibrotic_3.png"), width=2000, height=2000, units="px")
  print(plot)
  dev.off()
  
  ## QC metrics at cluster level ##
  cluster.qc.heatmap <- FeaturePlot(MP.harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "cluster.qc.heatmap.png"), width=1500,height=500,units="px")
  print(cluster.qc.heatmap)
  dev.off()
  
  # feature plot with MNP markers
  MNP.markers <- FeaturePlot(MP.harmony, features = MNP.genes, pt.size = 0.5, ncol = 6) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "MNP.markers.png"), width=1800,height=900,units="px")
  print(MNP.markers)
  dev.off()
  
  # feature plot with cDC 1 markers
  cDC1.markers <- FeaturePlot(MP.harmony, features = cDC1.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "cDC1.markers.png"), width=1600,height=400,units="px")
  print(cDC1.markers)
  dev.off()
  
  # feature plot with cDC 2 markers
  cDC2.markers <- FeaturePlot(MP.harmony, features = cDC2.genes, pt.size = 0.5, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path, "UMAP_dim", d, "cDC2.markers.png"), width=1200,height=800,units="px")
  print(cDC2.markers)
  dev.off()
  
  # feature plot with pDC markers
  pDC.markers <- FeaturePlot(MP.harmony, features = pDC.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "pDC.markers.png"), width=1600,height=400,units="px")
  print(pDC.markers)
  dev.off()
  
  # proliferating cells
  proliferation.markers <- FeaturePlot(MP.harmony, features = proliferation.genes, pt.size = 0.5, ncol = 2) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"proliferation.markers.png"), width=800,height=400,units="px")
  print(proliferation.markers)
  dev.off()
  
}

### clustering ###
res <- seq(0.5,1.5,0.1)
for(d in dims){
  
  # set path variable
  dim.path <- paste0(harmony.samples.path, "dim", d, "_annotation/")
  
  # calculate UMAP
  MP.harmony <- RunUMAP(MP.harmony, reduction = "harmony_theta2", dims=1:d, seed.use=1)
  
  # construct kNN graph
  MP.harmony <- FindNeighbors(MP.harmony, reduction = "harmony_theta2", dims = 1:d)
  
  # plot batch effect
  batch.plot <- DimPlot(MP.harmony, reduction = "umap", group.by = "condition", pt.size = 0.01)
  
  for (r in res) {
    # compute clusters
    MP.harmony <- FindClusters(MP.harmony, reduction = "harmony_theta2", resolution = r)
    
    # plot clusters
    cluster.plot <- DimPlot(MP.harmony, reduction = "umap", label = TRUE, pt.size = 0.1)
    
    # create eval plot
    eval.plot <- cluster.plot + batch.plot
    png(paste0(dim.path, "UMAP_dim", d, "_res", r, ".png"), width=1500, height=600, units="px")
    print(eval.plot)
    dev.off()
  }
}

### preliminary clustering (HVGs_overall) ###
MP.harmony <- FindNeighbors(MP.harmony, reduction = "harmony_theta2", dims = 1:50)
MP.harmony <- FindClusters(MP.harmony, reduction = "harmony_theta2", resolution = 0.7)
# run UMAP
MP.harmony <- RunUMAP(MP.harmony, reduction = "harmony_theta2", dims=1:50, seed.use=1)

### save data ###
saveRDS(MP.harmony, paste0(harmony.samples.path, "all_MP_harmony.rds"))

## compute cluster marker genes ###
MP.markers <- FindAllMarkers(MP.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
MP.top50.markers <- MP.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
write.csv(MP.markers, file = paste0(harmony.samples.path, "dim50_annotation/ALL_marker_genes.csv"))
write.csv(MP.top50.markers, file = paste0(harmony.samples.path, "dim50_annotation/top50_marker_genes.csv"))
