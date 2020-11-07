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
mesenchyme.path <- "/home/s1987963/ds_group/Niklas/combined_mesenchyme/combined_liver_mesenchyme_annotated.rds"
harmony.samples.path <- "/home/s1987963/ds_group/Niklas/mesenchyme_liver/harmonize_samples/HVGs_conditions/"
dir.create(harmony.samples.path, recursive = TRUE)

### load data ###
mesenchyme <- readRDS(mesenchyme.path) # 2728 cells

#### split healthy data ###
mesenchyme.healthy <- subset(mesenchyme, subset = condition == "healthy")
mesenchyme.fibrotic <- subset(mesenchyme, subset = condition == "fibrotic")

### feature selection: HVGs ###
#mesenchyme <- FindVariableFeatures(mesenchyme, selection.method = "vst", nfeatures = 2000)
mesenchyme.healthy <- FindVariableFeatures(mesenchyme.healthy, selection.method = "vst", nfeatures = 4000)
mesenchyme.fibrotic <- FindVariableFeatures(mesenchyme.fibrotic, selection.method = "vst", nfeatures = 4000)
mesenchyme.HVGs <- intersect(VariableFeatures(mesenchyme.fibrotic), VariableFeatures(mesenchyme.healthy)) # 1873
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

PC1_2.studies.plot <- DimPlot(object = mesenchyme, reduction = "pca", pt.size = .1, group.by = "dataset")
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

harmony.PC1_2.studies.plot <- DimPlot(object = mesenchyme.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "dataset")
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

# visualize batch effect and marker gene expression
dims <- c(24,33,41,50)
for(d in dims){
  
  # create folder
  dir.create(paste0(harmony.samples.path, "dim", d, "_annotation"))
  dim.path <- paste0(harmony.samples.path, "dim", d, "_annotation/")
  
  # run UMAP
  mesenchyme.harmony <- RunUMAP(mesenchyme.harmony, reduction = "harmony_theta2", dims = 1:d, seed.use=1)
  
  ## visualize batch effect 
  # no.1 overview
  sample.batch.plot <- DimPlot(mesenchyme.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
  dataset.batch.plot <- DimPlot(mesenchyme.harmony, reduction = "umap", group.by = "dataset", pt.size = 0.1)
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
  
  # no.3 - split by patients
  batch3.plot <- DimPlot(mesenchyme.harmony, reduction = "umap", group.by = "study", split.by = "patient.ID", pt.size = 0.01, ncol = 6)
  png(paste0(dim.path, "UMAP_dim", d, ".patient.batch2.png"), width=1800, height=600, units="px")
  print(batch3.plot)
  dev.off()
  
  # no.4 - colored by dataset
  batch4.plot <- DimPlot(mesenchyme.harmony, reduction = "umap", group.by = "dataset", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".study.batch1.png"), width=1600, height=1000, units="px")
  print(batch4.plot)
  dev.off()
  
  # no.5 - split by dataset
  batch5.plot <- DimPlot(mesenchyme.harmony, reduction = "umap", split.by = "dataset", pt.size = 0.01, ncol = 3)
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
  
  ## compare healthy vs. fibrotic ##
  compare <- DimPlot(mesenchyme.harmony, reduction = "umap", split.by = "condition", group.by = "condition", 
                     pt.size = 0.01) + NoLegend() + scale_color_viridis(discrete=TRUE)
  png(paste0(dim.path, "UMAP_dim", d, "healthy_fibrotic_1.png"), width=1800, height=900, units="px")
  print(compare)
  dev.off()
  
  ## compare original annotation ##
  annotation <- DimPlot(mesenchyme.harmony, reduction = "umap", split.by = "condition", group.by = "celltype_integrated",
                        label = TRUE, repel = TRUE, pt.size = 0.01) + NoLegend()
  png(paste0(dim.path, "UMAP_dim", d, "healthy_fibrotic_2.png"), width=1800, height=900, units="px")
  print(annotation)
  dev.off()
  
  ## compare original annotation in healthy vs. fibrotic data ##
  plot <- compare / annotation
  png(paste0(dim.path, "UMAP_dim", d, "healthy_fibrotic_3.png"), width=2000, height=2000, units="px")
  print(plot)
  dev.off()
  
  ## QC metrics at cluster level ##
  cluster.qc.heatmap <- FeaturePlot(mesenchyme.harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "cluster.qc.heatmap.png"), width=1500,height=500,units="px")
  print(cluster.qc.heatmap)
  dev.off()
  
  # feature plot with MNP markers
  mesenchyme.markers <- FeaturePlot(mesenchyme.harmony, features = marker.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "marker_genes.png"), width=1000,height=1000,units="px")
  print(mesenchyme.markers)
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

### preliminary clustering ###
mesenchyme.harmony <- FindNeighbors(mesenchyme.harmony, reduction = "harmony_theta2", dims = 1:50)
mesenchyme.harmony <- FindClusters(mesenchyme.harmony, reduction = "harmony_theta2", resolution = 0.8)
# run UMAP
mesenchyme.harmony <- RunUMAP(mesenchyme.harmony, reduction = "harmony_theta2", dims=1:50, seed.use=1)

### save data ###
saveRDS(mesenchyme.harmony, paste0(harmony.samples.path, "all_mesenchyme_harmony.rds"))

## compute cluster marker genes ###
mesenchyme.markers <- FindAllMarkers(mesenchyme.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
mesenchyme.top50.markers <- mesenchyme.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
write.csv(mesenchyme.markers, file = paste0(harmony.samples.path, "dim50_annotation/ALL_marker_genes.csv"))
write.csv(mesenchyme.top50.markers, file = paste0(harmony.samples.path, "dim50_annotation/top50_marker_genes.csv"))

### celltype annotation ###
cluster.annotation <- c("VSMC 1", "VSMC 2", "HSC 1", "VSMC 3", "Myofibroblasts", 
                        "VSMC 4", "Mesothelia", "HSC 2",
                        "Mes/Leuk doublet", "Mes/Hep doublet", "Mes/Endo doublet", "Mes/Cholangio")
names(cluster.annotation) <- levels(mesenchyme.harmony)
mesenchyme.harmony <- RenameIdents(mesenchyme.harmony, cluster.annotation)

# add cell types to meta data
cell.data <- data.table(barcode = colnames(mesenchyme.harmony),
                        celltype = Idents(mesenchyme.harmony))
cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$barcode <- NULL
mesenchyme.harmony <- AddMetaData(mesenchyme.harmony, cell.data, col.name = "mesenchyme_celltype")

# order cell type factors
mesenchyme.harmony$mesenchyme_celltype = factor(mesenchyme.harmony$mesenchyme_celltype)
mesenchyme.harmony$mesenchyme_celltype <- factor(mesenchyme.harmony$mesenchyme_celltype, levels(mesenchyme.harmony$mesenchyme_celltype)[c( 3,8, # HSC
                                                                                                                                           7, # Mesothelia
                                                                                                                                           5, # Myofibroblasts
                                                                                                                                           1,2,4,6, # VSMC
                                                                                                                                           10,12,11,9 # Doublets
                                                                                                                                           )])

### cell lineage annotation ###
cell.data <- data.table(barcode = colnames(mesenchyme.harmony),
                        celltype = Idents(mesenchyme.harmony))

lineage.annotation <- c("Mesenchyme", "Mesenchyme", "Mesenchyme", "Mesenchyme", 
                        "Mesenchyme", "Mesenchyme", "Mesenchyme", "Mesenchyme",
                        "Doublet", "Doublet", "Doublet", "Doublet")

lineage.data <- data.table(celltype = cluster.annotation, lineage = lineage.annotation)
meta.data <- merge(cell.data, lineage.data, by = "celltype")
meta.data <- data.frame(meta.data, row.names = meta.data$barcode)
meta.data$barcode <- NULL
meta.data$celltype <- NULL
mesenchyme.harmony <- AddMetaData(mesenchyme.harmony, meta.data, col.name = "lineage_integrated")

# order lineage factors
mesenchyme.harmony$lineage_integrated = factor(mesenchyme.harmony$lineage_integrated)
mesenchyme.harmony$lineage_integrated<- factor(mesenchyme.harmony$lineage_integrated, levels(mesenchyme.harmony$lineage_integrated)[c(2,1)])

### save data ###
saveRDS(mesenchyme.harmony, "/home/s1987963/ds_group/Niklas/combined_mesenchyme/combined_liver_mesenchyme_subclustered_annotated.rds")
