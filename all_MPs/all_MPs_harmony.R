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
library(scales)
library(viridis)

### path variables ###
liver.path <- "/home/s1987963/ds_group/Niklas/combined_MPs/combined_liver_MPs_annotated.rds"
lung.path <- "/home/s1987963/ds_group/Niklas/combined_MPs/combined_lung_MPs_annotated.rds"
kidney.path <- "/home/s1987963/ds_group/Niklas/combined_MPs/combined_kidney_MPs_annotated.rds"
blood.path <- "/home/s1987963/ds_group/Niklas/combined_MPs/combined_blood_MPs_annotated.rds"
MP.path <- "/home/s1987963/ds_group/Niklas/combined_organs_MPs/harmonize_samples_organs/HVGs_2organs/"

### read data ###
liver <- readRDS(liver.path)
lung <- readRDS(lung.path)
kidney <- readRDS(kidney.path)
blood <- readRDS(blood.path)

### correct blood meta data ###
### create meta cell types ###
cell.data <- data.table(barcode = colnames(blood),
                        condition = blood$condition)
cell.data[, condition_new := ifelse(condition == "Healthy", "healthy", "fibrotic")]
cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$barcode <- NULL
cell.data$condition <- NULL
blood <- AddMetaData(blood, cell.data, col.name = "condition")

# re-order condition factor levels 
# re-order condition factor #
blood$condition = factor(blood$condition)
blood$condition <- factor(blood$condition, levels(blood$condition)[c(2,1)])

### feature selection: HVGs ###
liver <- FindVariableFeatures(liver, selection.method = "vst", nfeatures = 4000)
lung <- FindVariableFeatures(lung, selection.method = "vst", nfeatures = 4000)
kidney <- FindVariableFeatures(kidney, selection.method = "vst", nfeatures = 4000)
blood <- FindVariableFeatures(blood, selection.method = "vst", nfeatures = 4000)
HVG.list <- c()
HVG.list <- append(HVG.list, VariableFeatures(liver))
HVG.list <- append(HVG.list, VariableFeatures(lung))
HVG.list <- append(HVG.list, VariableFeatures(kidney))
HVG.list <- append(HVG.list, VariableFeatures(blood))
# evaluate HVGs
overall.HVGs <- unique(HVG.list[which(table(HVG.list) > 1)]) # >1: 2962 HVGs, >2: 933 HVGs, >3: 162 HVGs

### merge data ###
MP <- merge(liver, c(lung, kidney, blood)) # 148072 cells

### set HVGs ###
# set HVGs
VariableFeatures(MP) <- overall.HVGs

### scale data ###
MP <- ScaleData(MP) # uncorrected

### dimensionality reduction: PCA ###
MP <- RunPCA(MP, features = VariableFeatures(object = MP))

# elbow plot
pca.elbow.plot <- ElbowPlot(MP, ndims = 50, reduction = "pca")
png(paste0(MP.path,"pca.elbow.plot.png"), width=1000,height=600,units="px")
print(pca.elbow.plot)
dev.off()

# visualize PCs
PC1_2.samples.plot <- DimPlot(object = MP, reduction = "pca", pt.size = .1, group.by = "patient.ID")
png(paste0(MP.path,"PC1_2.samples.png"), width=1000,height=1000,units="px")
print(PC1_2.samples.plot)
dev.off()

PC1_2.studies.plot <- DimPlot(object = MP, reduction = "pca", pt.size = .1, group.by = "study")
png(paste0(MP.path,"PC1_2.studies.png"), width=1000,height=1000,units="px")
print(PC1_2.studies.plot)
dev.off()

PC1_2.organs.plot <- DimPlot(object = MP, reduction = "pca", pt.size = .1, group.by = "organ")
png(paste0(MP.path,"PC1_2.organs.png"), width=1000,height=1000,units="px")
print(PC1_2.organs.plot)
dev.off()

PC1_2.condition.plot <- DimPlot(object = MP, reduction = "pca", pt.size = .1, group.by = "condition")
png(paste0(MP.path,"PC1_2.condition.png"), width=1000,height=1000,units="px")
print(PC1_2.condition.plot)
dev.off()

# save genes making up the PCs  
sink(paste0(MP.path, "PC_genes.txt"))
print(MP[["pca"]], dims = 1:50, nfeatures = 20)
sink()

### integration with harmony ###
# integrate samples
MP.harmony <- MP %>% RunHarmony(c("patient.ID", "organ"), theta = c(2,2), reduction.save = "harmony_theta2", plot_convergence = TRUE) # harmonize organs

# harmony elbow plot
harmony.elbow.plot <- ElbowPlot(MP.harmony, ndims = 50, reduction = "harmony_theta2")
png(paste0(MP.path,"harmony_theta2.elbow.plot.png"), width=1000,height=600,units="px")
print(harmony.elbow.plot)
dev.off()

# explore harmony coordinates
# visualize PCs 
harmony.PC1_2.samples.plot <- DimPlot(object = MP.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "patient.ID")
png(paste0(MP.path,"harmony_theta2.PC1_2.samples.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.samples.plot)
dev.off()

harmony.PC1_2.studies.plot <- DimPlot(object = MP.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "study")
png(paste0(MP.path,"harmony_theta2.PC1_2.studies.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.studies.plot)
dev.off()

harmony.PC1_2.organs.plot <- DimPlot(object = MP.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "organ")
png(paste0(MP.path,"harmony_theta2.PC1_2.organs.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.organs.plot)
dev.off()

harmony.PC1_2.condition.plot <- DimPlot(object = MP.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "condition")
png(paste0(MP.path,"harmony_theta2.PC1_2.condition.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.condition.plot)
dev.off()

# save genes making up the PCs
sink(paste0(MP.path, "harmony_PC_genes.txt"))
print(MP.harmony[["harmony_theta2"]], dims = 1:50, nfeatures = 20)
sink()

### save data ###
saveRDS(MP.harmony, paste0(MP.path, "MP_harmony.rds"))

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
#dims <- c(22,31,38,45,50) # harmony corrected PCA
dims <- c(12,23,31,38,42,48,50)
for(d in dims){
  
  # create folder
  dir.create(paste0(MP.path, "dim", d, "_annotation"))
  dim.path <- paste0(MP.path, "dim", d, "_annotation/")
  
  # run UMAP
  MP.harmony <- RunUMAP(MP.harmony, reduction = "harmony_theta2", dims = 1:d, seed.use=1)
  
  # no.1 overview
  study.batch.plot <- DimPlot(MP.harmony, reduction = "umap", group.by = "study", pt.size = 0.1)
  organ.batch.plot <- DimPlot(MP.harmony, reduction = "umap", group.by = "organ", pt.size = 0.1)
  condition.batch.plot <-  DimPlot(MP.harmony, reduction = "umap", group.by = "condition", pt.size = 0.1)
  batch1.plot <- study.batch.plot + organ.batch.plot + condition.batch.plot
  png(paste0(dim.path, "UMAP_dim", d, ".batch.png"), width=1800, height=600, units="px")
  print(batch1.plot)
  dev.off()
  
  ## no.2 - colored by patients
  #batch2.plot <- DimPlot(MP.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.01)
  #png(paste0(dim.path, "UMAP_dim", d, ".patient.batch1.png"), width=1600, height=1000, units="px")
  #print(batch2.plot)
  #dev.off()
  #
  ## no.3 - split by patients
  #batch3.plot <- DimPlot(MP.harmony, reduction = "umap", group.by = "study", split.by = "patient.ID", pt.size = 0.01, ncol = 6)
  #png(paste0(dim.path, "UMAP_dim", d, ".patient.batch2.png"), width=1800, height=1500, units="px")
  #print(batch3.plot)
  #dev.off()
  
  # no.4 - colored by organ
  batch4.plot <- DimPlot(MP.harmony, reduction = "umap", group.by = "organ", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".study.batch1.png"), width=1600, height=1000, units="px")
  print(batch4.plot)
  dev.off()
  
  # no.5 - split by organ
  batch5.plot <- DimPlot(MP.harmony, reduction = "umap", group.by = "study", split.by = "organ", pt.size = 0.01, ncol = 4)
  png(paste0(dim.path, "UMAP_dim", d, ".study.batch2.png"), width=1800, height=900, units="px")
  print(batch5.plot)
  dev.off()
  
  # no.6 - colored by lineage
  batch6.plot <- DimPlot(MP.harmony, reduction = "umap", group.by = "celltype_integrated", label=TRUE, repel=TRUE, pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".celltype.batch1.png"), width=1600, height=1000, units="px")
  print(batch6.plot)
  dev.off()
  
  # no.7 - colored by lineage, split by organ
  batch7.plot <- DimPlot(MP.harmony, reduction = "umap", split.by = "organ", group.by = "celltype_integrated", label=TRUE, repel=TRUE, pt.size = 0.01, ncol=2) + NoLegend()
  png(paste0(dim.path, "UMAP_dim", d, ".celltype.batch2.png"), width=1200, height=1200, units="px")
  print(batch7.plot)
  dev.off()
  
  # no.8 - lineage overview
  lineage.overview <- organ.batch.plot + batch6.plot
  png(paste0(dim.path, "UMAP_dim", d, ".lineage.overview.png"), width=2000, height=1000, units="px")
  print(lineage.overview)
  dev.off()
  
  ## QC metrics at cluster level ##
  cluster.qc.heatmap <- FeaturePlot(MP.harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "cluster.qc.heatmap.png"), width=1500,height=500,units="px")
  print(cluster.qc.heatmap)
  dev.off()
  
  ## gene expression heatmaps ##
  # feature plot with immune cell marker
  immunecell.markers <- FeaturePlot(MP.harmony, features = c("PTPRC"), pt.size = 0.5, ncol = 1) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "immunecell.markers.png"), width=1000,height=1000,units="px")
  print(immunecell.markers)
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
res <- seq(0.1,1.5,0.1)
for(d in dims){
  
  # set path variable
  dim.path <- paste0(MP.path, "dim", d, "_annotation/")
  
  # visualization with UMAP
  MP.harmony <- RunUMAP(MP.harmony, reduction = "harmony_theta2", dims=1:d, seed.use=1)
  
  # plot batch effect
  batch.plot <- DimPlot(MP.harmony, reduction = "umap", group.by = "organ", pt.size = 0.1)
  
  # create kNN graph
  MP.harmony <- FindNeighbors(MP.harmony, reduction = "harmony_theta2", dims = 1:d)
  
  for (r in res) {
    
    MP.harmony <- FindClusters(MP.harmony, reduction = "harmony_theta2", resolution = r)
    umap.plot <- DimPlot(MP.harmony, reduction = "umap", label = TRUE, pt.size = 0.1)
    
    # create eval plot
    eval.plot <- umap.plot + batch.plot
    png(paste0(dim.path, "UMAP_dim", d, "_res", r, ".png"), width=1800, height=800, units="px")
    print(eval.plot)
    dev.off()
    
    ## clustering at patient level
    #patient.clustering <- DimPlot(MP.harmony, group.by = "seurat_clusters", split.by = "patient.ID", pt.size = 0.1, ncol = 10)
    #png(paste0(dim.path,"UMAP_dim", d, "_res", r,".patient.clustering.png"), width=2000,height=400,units="px")
    #print(patient.clustering)
    #dev.off()
    
  }
}

3

# add metadata for visualisation no.1
cell.data <- data.table(barcode = colnames(MP.harmony),
                        celltype = MP.harmony$celltype_integrated)
cell.data[, celltype_visualisation := ifelse(celltype %in% c("Combined Liver Macrophage 1", "Combined Lung Macrophage 1", "Combined Kidney Macrophage"), paste0(celltype), "Other")]
cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$celltype <- NULL
cell.data$barcode <- NULL
MP.harmony <- AddMetaData(MP.harmony, cell.data, col.name = "celltype_visualisation_1")

conserved.plot1 <- DimPlot(MP.harmony, reduction = "umap", 
                           group.by = "celltype_visualisation_1",
                           pt.size = 0.1, cols = c(viridis(3),"#999999"),
                           label = TRUE, repel = TRUE) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 12))
png(paste0(MP.path, "dim50_annotation/conserved_celltype_1_umap.png"), width=1000,height=1000,units="px")
print(conserved.plot1)
dev.off()

# add metadata for visualisation no.2
cell.data <- data.table(barcode = colnames(MP.harmony),
                        celltype = MP.harmony$celltype_integrated)
cell.data[, celltype_visualisation := ifelse(celltype %in% c("Combined Liver KC1", "Combined Liver KC3", "Combined Lung Alveolar Macrophage 1", "Combined Lung Alveolar Macrophage 3"), paste0(celltype), "Other")]
cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$celltype <- NULL
cell.data$barcode <- NULL
MP.harmony <- AddMetaData(MP.harmony, cell.data, col.name = "celltype_visualisation_2")

conserved.plot2 <- DimPlot(MP.harmony, reduction = "umap", group.by = "celltype_visualisation_2", cols = c(viridis(4),"#999999"), label = TRUE, repel = TRUE, pt.size = 0.1) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 12))
png(paste0(MP.path, "dim50_annotation/conserved_celltype_2_umap.png"), width=1000,height=1000,units="px")
print(conserved.plot2)
dev.off()
