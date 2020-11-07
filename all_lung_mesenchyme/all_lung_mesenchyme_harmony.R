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
mesenchyme.path <- "/home/s1987963/ds_group/Niklas/combined_mesenchyme/combined_lung_mesenchyme_annotated.rds"
harmony.samples.path <- "/home/s1987963/ds_group/Niklas/mesenchyme_lung/harmonize_samples/HVGs_conditions/"
dir.create(harmony.samples.path, recursive = TRUE)

### load data ###
mesenchyme <- readRDS(mesenchyme.path) # 9574 cells

#### split healthy data ###
mesenchyme.healthy <- subset(mesenchyme, subset = condition == "healthy")
mesenchyme.fibrotic <- subset(mesenchyme, subset = condition == "fibrotic")

### feature selection: HVGs ###
#mesenchyme <- FindVariableFeatures(mesenchyme, selection.method = "vst", nfeatures = 2000)
mesenchyme.healthy <- FindVariableFeatures(mesenchyme.healthy, selection.method = "vst", nfeatures = 4000)
mesenchyme.fibrotic <- FindVariableFeatures(mesenchyme.fibrotic, selection.method = "vst", nfeatures = 4000)
mesenchyme.HVGs <- intersect(VariableFeatures(mesenchyme.fibrotic), VariableFeatures(mesenchyme.healthy)) # 1965
VariableFeatures(mesenchyme) <- mesenchyme.HVGs

### scale data ###
mesenchyme <- ScaleData(mesenchyme) # uncorrected

### dimensionality reduction: PCA ###
mesenchyme <- RunPCA(mesenchyme, features = VariableFeatures(object = mesenchyme))

# elbow plot
pca.elbow.plot <- ElbowPlot(mesenchyme, ndims = 50, reduction = "pca")
png(paste0(harmony.samples.path,"pca.elbow.plot.png"), width=1000,height=600,units="px")
print(pca.elbow.plot)
dev.off()

# visualize PCs
PC1_2.samples.plot <- DimPlot(object = mesenchyme, reduction = "pca", pt.size = .1, group.by = "patient.ID")
png(paste0(harmony.samples.path,"PC1_2.samples.png"), width=1000,height=1000,units="px")
print(PC1_2.samples.plot)
dev.off()

PC1_2.studies.plot <- DimPlot(object = mesenchyme, reduction = "pca", pt.size = .1, group.by = "study")
png(paste0(harmony.samples.path,"PC1_2.studies.png"), width=1000,height=1000,units="px")
print(PC1_2.studies.plot)
dev.off()

# save genes making up the PCs  
sink(paste0(harmony.samples.path, "PC_genes.txt"))
print(mesenchyme[["pca"]], dims = 1:50, nfeatures = 20)
sink()

### integration with harmony ###
# harmonize samples
mesenchyme.harmony <- mesenchyme %>% RunHarmony("patient.ID", theta = 2, reduction.save = "harmony_theta2", plot_convergence = TRUE) # harmonize all 12 samples independently

## harmony elbow plot ##
harmony.elbow.plot <- ElbowPlot(mesenchyme.harmony, ndims = 50, reduction = "harmony_theta2")
png(paste0(harmony.samples.path,"harmony_theta2.elbow.plot.png"), width=1000,height=600,units="px")
print(harmony.elbow.plot)
dev.off()

## explore harmony coordinates ##
## visualize PCs ##
harmony.PC1_2.samples.plot <- DimPlot(object = mesenchyme.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "patient.ID")
png(paste0(harmony.samples.path,"harmony_theta2.PC1_2.samples.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.samples.plot)
dev.off()

harmony.PC1_2.studies.plot <- DimPlot(object = mesenchyme.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "study")
png(paste0(harmony.samples.path,"harmony_theta2.PC1_2.studies.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.studies.plot)
dev.off()

### save genes making up the PCs ### 
sink(paste0(harmony.samples.path, "harmony_PC_genes.txt"))
print(mesenchyme.harmony[["harmony_theta2"]], dims = 1:50, nfeatures = 20)
sink()

### save data ###
saveRDS(mesenchyme.harmony, paste0(harmony.samples.path, "all_mesenchyme_harmony.rds"))

### explore different numbers of harmony PCs ###
# mesenchymal and proliferation marker genes
marker.genes <- c("COL1A1",	"COL3A1", "ACTA2", "MYH11",	"PDGFRA", 
                  "PDGFRB", "RGS5",	"MSLN", "DCN", "CD34",
                  "THY11", "VIPR1", "PTHR1", "MKI67",	"TOP2A",
                  "NGFR", "ADAMTSL2", "MFAP4")

adv.fib.genes <- c("PI16", "MFAP5", "SFRP1", "SFRP2")

inflammatory.fib.genes <- c("CXCL12", "FGF7", "COL3A1", "SCG2")

lipo.fib.genes <- c("ITGA8", "NPNT")
myo.fib.genes <- c("DIO2", "ASPN")
lipo.myo.transition.genes <- c("FBLN2", "PIEZO2")
mesothelia.genes <- c("KRT19", "SYT4")
pericyte.genes <- c("COX4I2", "PDGRB", "ITGA7", "SSTR2")
SMC.genes <- c("ACTA2", "MYH11")

# visualize batch effect and marker gene expression
dims <- c(20,31, 36,39,45,50)
for(d in dims){
  
  # create folder
  dir.create(paste0(harmony.samples.path, "dim", d, "_annotation"))
  dim.path <- paste0(harmony.samples.path, "dim", d, "_annotation/")
  
  # run UMAP
  mesenchyme.harmony <- RunUMAP(mesenchyme.harmony, reduction = "harmony_theta2", dims = 1:d, seed.use=1)
  
  ## visualize batch effect 
  # no.1 overview
  sample.batch.plot <- DimPlot(mesenchyme.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
  dataset.batch.plot <- DimPlot(mesenchyme.harmony, reduction = "umap", group.by = "study", pt.size = 0.1)
  condition.batch.plot <- DimPlot(mesenchyme.harmony, reduction = "umap", group.by = "condition", pt.size = 0.1)
  batch1.plot <- sample.batch.plot + dataset.batch.plot + condition.batch.plot
  png(paste0(dim.path, "UMAP_dim", d, ".batch.png"), width=1800, height=600, units="px")
  print(batch1.plot)
  dev.off()
  
  # no.2 - colored by patients
  batch2.plot <- DimPlot(mesenchyme.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".patient.batch1.png"), width=1600, height=1000, units="px")
  print(batch2.plot)
  dev.off()
  
  # no.4 - colored by dataset
  batch4.plot <- DimPlot(mesenchyme.harmony, reduction = "umap", group.by = "study", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".study.batch1.png"), width=1600, height=1000, units="px")
  print(batch4.plot)
  dev.off()
  
  # no.5 - split by dataset
  batch5.plot <- DimPlot(mesenchyme.harmony, reduction = "umap", split.by = "study", pt.size = 0.01, ncol = 3)
  png(paste0(dim.path, "UMAP_dim", d, ".study.batch2.png"), width=1800, height=600, units="px")
  print(batch5.plot)
  dev.off()
  
  # no.6 - colored by condition
  batch6.plot <- DimPlot(mesenchyme.harmony, reduction = "umap", group.by = "condition", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".condition.batch1.png"), width=1000, height=1000, units="px")
  print(batch6.plot)
  dev.off()
  
  # no.7 - split by condition
  batch7.plot <- DimPlot(mesenchyme.harmony, reduction = "umap", group.by = "study", split.by = "condition", pt.size = 0.01, ncol = 2)
  png(paste0(dim.path, "UMAP_dim", d, ".condition.batch2.png"), width=1800, height=900, units="px")
  print(batch7.plot)
  dev.off()
  
  # no.9 - colored by lineage annotation
  batch9.plot <- DimPlot(mesenchyme.harmony, reduction = "umap", group.by = "lineage_integrated", label=TRUE, pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".lineage.batch2.png"), width=1600, height=1000, units="px")
  print(batch9.plot)
  dev.off()
  
  ## QC metrics at cluster level ##
  cluster.qc.heatmap <- FeaturePlot(mesenchyme.harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "cluster.qc.heatmap.png"), width=1500,height=500,units="px")
  print(cluster.qc.heatmap)
  dev.off()
  
  # feature plot with mesenchymal markers
  mesenchyme.markers <- FeaturePlot(mesenchyme.harmony, features = marker.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "marker_genes.png"), width=1000,height=1000,units="px")
  print(mesenchyme.markers)
  dev.off()
  
  # feature plot with adventitial fibroblast markers
  adv.fib.markers <- FeaturePlot(mesenchyme.harmony, features = adv.fib.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "_adventitial_fibroblasts_markers.png"), width=1200,height=400,units="px")
  print(adv.fib.markers)
  dev.off()
  
  # feature plot with inflammatory fibroblast markers
  inflammatory.fib.markers <- FeaturePlot(mesenchyme.harmony, features = inflammatory.fib.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "_inflammatory_fibroblasts_markers.png"), width=1200,height=400,units="px")
  print(inflammatory.fib.markers)
  dev.off()
  
  # feature plot with lipofibroblast markers
  lipo.fib.markers <- FeaturePlot(mesenchyme.harmony, features = lipo.fib.genes, pt.size = 0.5, ncol = 2) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "_lipo_fibroblasts_markers.png"), width=800,height=400,units="px")
  print(lipo.fib.markers)
  dev.off()
  
  # feature plot with myofibroblast markers
  myo.fib.markers <- FeaturePlot(mesenchyme.harmony, features = myo.fib.genes, pt.size = 0.5, ncol = 2) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "_myo_fibroblasts_markers.png"), width=800,height=400,units="px")
  print(myo.fib.markers)
  dev.off()
  
  # feature plot with lipo-myofibroblast transition markers
  lipo.myo.transition.markers <- FeaturePlot(mesenchyme.harmony, features = lipo.myo.transition.genes, pt.size = 0.5, ncol = 2) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "_lipo_myo_transition_markers.png"), width=800,height=400,units="px")
  print(lipo.myo.transition.markers)
  dev.off()
  
  # feature plot with mesothelia markers
  mesothelia.markers <- FeaturePlot(mesenchyme.harmony, features = mesothelia.genes, pt.size = 0.5, ncol = 2) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "_mesothelia_markers.png"), width=800,height=400,units="px")
  print(mesothelia.markers)
  dev.off()
  
  # feature plot with pericyte markers
  pericyte.markers <- FeaturePlot(mesenchyme.harmony, features = pericyte.genes, pt.size = 0.5, ncol = 2) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "_pericyte_markers.png"), width=1200,height=400,units="px")
  print(pericyte.markers)
  dev.off()
  
  # feature plot with SMC markers
  SMC.markers <- FeaturePlot(mesenchyme.harmony, features = SMC.genes, pt.size = 0.5, ncol = 2) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "_SMC_markers.png"), width=800, height=400,units="px")
  print(SMC.markers)
  dev.off()
  
}

### clustering ###
res <- seq(0.1,1.5,0.1)
for(d in dims){
  
  # set path variable
  dim.path <- paste0(harmony.samples.path, "dim", d, "_annotation/")
  
  # calculate UMAP
  mesenchyme.harmony <- RunUMAP(mesenchyme.harmony, reduction = "harmony_theta2", dims=1:d, seed.use=1)
  
  # construct kNN graph
  mesenchyme.harmony <- FindNeighbors(mesenchyme.harmony, reduction = "harmony_theta2", dims = 1:d)
  
  # plot batch effect
  batch.plot <- DimPlot(mesenchyme.harmony, reduction = "umap", group.by = "condition", pt.size = 0.01)
  
  for (r in res) {
    # compute clusters
    mesenchyme.harmony <- FindClusters(mesenchyme.harmony, reduction = "harmony_theta2", resolution = r)
    
    # plot clusters
    cluster.plot <- DimPlot(mesenchyme.harmony, reduction = "umap", label = TRUE, pt.size = 0.1)
    
    # create eval plot
    eval.plot <- cluster.plot + batch.plot
    png(paste0(dim.path, "UMAP_dim", d, "_res", r, ".png"), width=1500, height=600, units="px")
    print(eval.plot)
    dev.off()
  }
}

### preliminary clustering (HVGs_overall) ###
mesenchyme.harmony <- FindNeighbors(mesenchyme.harmony, reduction = "harmony_theta2", dims = 1:50)
mesenchyme.harmony <- FindClusters(mesenchyme.harmony, reduction = "harmony_theta2", resolution = 0.7)
# run UMAP
mesenchyme.harmony <- RunUMAP(mesenchyme.harmony, reduction = "harmony_theta2", dims=1:50, seed.use=1)

### save data ###
saveRDS(mesenchyme.harmony, paste0(harmony.samples.path, "all_mesenchyme_harmony.rds"))

## compute cluster marker genes ###
mesenchyme.markers <- FindAllMarkers(mesenchyme.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
mesenchyme.top50.markers <- mesenchyme.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
write.csv(mesenchyme.markers, file = paste0(harmony.samples.path, "dim50_annotation/ALL_marker_genes.csv"))
write.csv(mesenchyme.top50.markers, file = paste0(harmony.samples.path, "dim50_annotation/top50_marker_genes.csv"))
