library(Seurat)
library(harmony)
library(future)
library(reticulate)
library(umap)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(scales)
library(patchwork)
library(data.table)

### path variables ###
healthy.lung.path <- "/home/s1987963/ds_group/Niklas/healthy_organs/healthy_lung_annotated.rds"
fibrotic.lung.path <- "/home/s1987963/ds_group/Niklas/fibrotic_organs/fibrotic_lung_annotated.rds"
lung.path <- "/home/s1987963/ds_group/Niklas/all_lung/"
harmony.samples.path <- "/home/s1987963/ds_group/Niklas/all_lung/harmonize_samples/HVGs_conditions/"

### load data ###
healthy.lung <- readRDS(healthy.lung.path) # 104884 cells
fibrotic.lung <- readRDS(fibrotic.lung.path) # 146046 cells

### merge data ###
lung <- merge(healthy.lung, fibrotic.lung) # 250930 cells

# add cell types to meta data
cell.data <- data.table(barcode = colnames(lung),
                        celltype = Idents(lung))
cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$barcode <- NULL
lung <- AddMetaData(lung, cell.data, col.name = "celltype")

### exclude 5134 doublets from downstream analyses ###
lung <- subset(lung, subset = lineage != "Doublet") # 245796 cells

### feature selection: HVGs ###
healthy.lung <- FindVariableFeatures(healthy.lung, selection.method = "vst", nfeatures = 4000)
fibrotic.lung <- FindVariableFeatures(fibrotic.lung, selection.method = "vst", nfeatures = 4000)
lung.HVGs <- intersect(VariableFeatures(fibrotic.lung), VariableFeatures(healthy.lung)) # 2042
VariableFeatures(lung) <- lung.HVGs

### scale data ###
lung <- ScaleData(lung) # uncorrected

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

PC1_2.condition.plot <- DimPlot(object = lung, reduction = "pca", pt.size = .1, group.by = "condition")
png(paste0(harmony.samples.path,"PC1_2.condition.png"), width=1000,height=1000,units="px")
print(PC1_2.condition.plot)
dev.off()

# save genes making up the PCs  
sink(paste0(harmony.samples.path, "PC_genes.txt"))
print(lung[["pca"]], dims = 1:50, nfeatures = 20)
sink()

### integration with harmony ###
# harmonize samples
lung.harmony <- lung %>% RunHarmony("patient.ID", theta = 2, reduction.save = "harmony_theta2", plot_convergence = TRUE) # harmonize all 60 samples independently

## harmony elbow plot ##
harmony.elbow.plot <- ElbowPlot(lung.harmony, ndims = 50, reduction = "harmony_theta2")
png(paste0(harmony.samples.path,"harmony_theta2.elbow.plot.png"), width=1000,height=600,units="px")
print(harmony.elbow.plot)
dev.off()

## explore harmony coordinates ##
## visualize PCs ##
harmony.PC1_2.samples.plot <- DimPlot(object = lung.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "patient.ID")
png(paste0(harmony.samples.path,"harmony_theta2.PC1_2.samples.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.samples.plot)
dev.off()

harmony.PC1_2.condition.plot <- DimPlot(object = lung.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "condition")
png(paste0(harmony.samples.path,"harmony_theta2.PC1_2.condition.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.condition.plot)
dev.off()

harmony.PC1_2.studies.plot <- DimPlot(object = lung.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "study")
png(paste0(harmony.samples.path,"harmony_theta2.PC1_2.studies.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.studies.plot)
dev.off()

### save genes making up the PCs ### 
sink(paste0(harmony.samples.path, "harmony_PC_genes.txt"))
print(lung.harmony[["harmony_theta2"]], dims = 1:50, nfeatures = 20)
sink()

### save data ###
saveRDS(lung.harmony, paste0(harmony.samples.path, "all_lung_harmony.rds"))

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

# visualize marker gene expression
dims <- c(50,40,31)
for(d in dims){
  
  # create folder
  dir.create(paste0(harmony.samples.path, "dim", d, "_annotation"))
  dim.path <- paste0(harmony.samples.path, "dim", d, "_annotation/")
  
  # run UMAP
  lung.harmony <- RunUMAP(lung.harmony, reduction = "harmony_theta2", dims = 1:d, seed.use=1)
  
  ## visualize batch effect 
  # no.1 overview
  sample.batch.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
  condition.batch.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = "condition", pt.size = 0.1)
  batch1.plot <- sample.batch.plot + condition.batch.plot
  png(paste0(dim.path, "UMAP_dim", d, ".batch.png"), width=1800, height=900, units="px")
  print(batch1.plot)
  dev.off()
  
  # no.2 - colored by patients
  batch2.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".patient.batch1.png"), width=1600, height=1000, units="px")
  print(batch2.plot)
  dev.off()
  
  # no.3 - split by patients
  batch3.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = "study", split.by = "patient.ID", pt.size = 0.01, ncol = 6)
  png(paste0(dim.path, "UMAP_dim", d, ".patient.batch2.png"), width=1800, height=600, units="px")
  print(batch3.plot)
  dev.off()
  
  # no.4 - colored by study
  batch4.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = "study", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".study.batch1.png"), width=1600, height=1000, units="px")
  print(batch4.plot)
  dev.off()
  
  # no.5 - split by study
  batch5.plot <- DimPlot(lung.harmony, reduction = "umap", split.by = "study", pt.size = 0.01, ncol = 3)
  png(paste0(dim.path, "UMAP_dim", d, ".study.batch2.png"), width=1800, height=600, units="px")
  print(batch5.plot)
  dev.off()
  
  # no.6 - colored by condition
  batch6.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = "condition", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".condition.batch1.png"), width=1600, height=1000, units="px")
  print(batch6.plot)
  dev.off()
  
  # no.7 - split by condition
  batch7.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = "study", split.by = "condition", pt.size = 0.01, ncol = 2)
  png(paste0(dim.path, "UMAP_dim", d, ".condition.batch2.png"), width=1800, height=900, units="px")
  print(batch7.plot)
  dev.off()
  
  # no.8 - colored by lineage
  batch8.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = "lineage", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".lineage.batch1.png"), width=1600, height=1000, units="px")
  print(batch8.plot)
  dev.off()
  
  # no.9 - colored by annotation
  batch9.plot <- DimPlot(lung.harmony, reduction = "umap", group.by = "ident", label=TRUE, pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".lineage.batch2.png"), width=1600, height=1000, units="px")
  print(batch9.plot)
  dev.off()
  
  # no.10 - lineage overview
  lineage.overview <- condition.batch.plot + batch8.plot
  png(paste0(dim.path, "UMAP_dim", d, ".lineage.overview.png"), width=2000, height=1000, units="px")
  print(lineage.overview)
  dev.off()
  
  ## compare healthy vs. fibrotic ##
  compare <- DimPlot(lung.harmony, reduction = "umap", split.by = "condition", group.by = "condition", 
                     pt.size = 0.01) + NoLegend() + scale_color_viridis(discrete=TRUE)
  png(paste0(dim.path, "UMAP_dim", d, "healthy_fibrotic_1.png"), width=1800, height=900, units="px")
  print(compare)
  dev.off()
  
  ## compare original annotation ##
  annotation <- DimPlot(lung.harmony, reduction = "umap", split.by = "condition", group.by = "celltype",
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
res <- seq(0.5,2,0.1)
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
    umap.plot <- DimPlot(lung.harmony, reduction = "umap", label = TRUE, pt.size = 0.1)
    
    # create eval plot
    eval.plot <- umap.plot + batch.plot
    png(paste0(dim.path, "UMAP_dim", d, "_res", r, ".png"), width=1500, height=600, units="px")
    print(eval.plot)
    dev.off()
    
    ## clustering at patient level
    #patient.clustering <- DimPlot(lung.harmony, group.by = "seurat_clusters", split.by = "patient.ID", pt.size = 0.1, ncol = 10)
    #png(paste0(dim.path,"UMAP_dim", d, "_res", r,".patient.clustering.png"), width=2000,height=400,units="px")
    #print(patient.clustering)
    #dev.off()
    
  }
}

### preliminary clustering ###
lung.harmony <- FindNeighbors(lung.harmony, reduction = "harmony_theta2", dims = 1:50)
lung.harmony <- FindClusters(lung.harmony, reduction = "harmony_theta2", resolution = 1.3)
# run UMAP
lung.harmony <- RunUMAP(lung.harmony, reduction = "harmony_theta2", dims=1:50, seed.use=1)

### save data ###
saveRDS(lung.harmony, paste0(harmony.samples.path, "all_lung_harmony.rds"))

# compare PTPRC expression across clusters
umap.plot <- DimPlot(lung.harmony, reduction = "umap", label = T, label.size = 6, pt.size = 0.1)
ptprc.plot <- VlnPlot(object = lung.harmony, features = c("PTPRC"), group.by = "seurat_clusters", pt.size = 0.1) + NoLegend()
immunecell.markers <- FeaturePlot(lung.harmony, features = c("PTPRC"), pt.size = 0.5, ncol = 1) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
immuneclusters.plot <- umap.plot + immunecell.markers - ptprc.plot + plot_layout(ncol=1, widths=c(2,1))
png(paste0(harmony.samples.path, "dim50_annotation/immuneclusters.png"), width=1800,height=1200,units="px")
print(immuneclusters.plot)
dev.off()

## compute cluster marker genes ###
lung.markers <- FindAllMarkers(lung.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
lung.top50.markers <- lung.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
write.csv(lung.markers, file = paste0(harmony.samples.path, "dim50_annotation/ALL_marker_genes.csv"))
write.csv(lung.top50.markers, file = paste0(harmony.samples.path, "dim50_annotation/top50_marker_genes.csv"))

### cell type annotation ###
cluster.annotation <- c("Combined Lung Alveolar Macrophage 1", "Combined Lung ATII 1", "Combined Lung Ciliated 1", 
                        "Combined Lung ATII 2", "Combined Lung Monocyte", "Combined Lung Alveolar Macrophage 2", 
                        "Combined Lung Venous/Peribronchial Endo", "Combined Lung Goblet 1", "Combined Lung Macrophage 2",
                        "Combined Lung CD4+ T cell", "Combined Lung Macrophage 1", "Combined Lung cDC2", 
                        "Combined Lung CD8+ T cell", "Combined Lung Mesenchyme 1", "Combined Lung Alveolar Macrophage 3",
                        "Combined Lung ATI", "Combined Lung NK cell", "Combined Lung Proliferating", 
                        "Combined Lung Basal cell 1", "Combined Lung Macrophage 3", "Combined Lung Mesenchyme 2",
                        "Combined Lung ArtEndo", "Combined Lung Alveolar Macrophage 4", "Combined Lung Ciliated 2",
                        "Combined Lung LymphEndo", "Combined Lung CapEndo 1", "Combined Lung B cell 1",
                        "Combined Lung Plasma cell 1", "Combined Lung Mast cell 1", "Combined Lung Ciliated 3", 
                        "Combined Lung doublet", "Combined Lung NK/MP doublet", "Combined Lung ATII 3", 
                        "Combined Lung Macrophage 4", "Combined Lung IFN-primed Macrophage", "Combined Lung Basal cell 2", 
                        "Combined T cell/plasma cell doublet", "Combined Lung ATII 4", "Combined Lung pDC", 
                        "Combined Lung Macrophage 5", "Combined Lung Plasma cell 2", "Combined Lung CapEndo 2",
                        "Combined Lung Plasma cell 3", "Combined Lung B cell 2", "Combined Lung Mesothelia",
                        "Combined Lung Mast cell 2", "Combined Lung Basal cell 3")
names(cluster.annotation) <- levels(lung.harmony)
lung.harmony <- RenameIdents(lung.harmony, cluster.annotation)

# add cell types to meta data
cell.data <- data.table(barcode = colnames(lung.harmony),
                        celltype = Idents(lung.harmony))
cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$barcode <- NULL
lung.harmony <- AddMetaData(lung.harmony, cell.data, col.name = "celltype_integrated")

# order cell type factors
lung.harmony$celltype_integrated = factor(lung.harmony$celltype_integrated)
lung.harmony$celltype_integrated <- factor(lung.harmony$celltype_integrated, levels(lung.harmony$celltype_integrated)[c(16, # AT I
                                                                                                               2,4,33,38, # AT II 1-4
                                                                                                               19,36,47, # Basal 1-3
                                                                                                               3,24,30, # Ciliated 1-3
                                                                                                               8, # Goblet
                                                                                                               22, # ArtEndo
                                                                                                               26, 42, # CapEndo 1-2
                                                                                                               7, # Venous/Peribronchial Endo
                                                                                                               25, # LymphEndo
                                                                                                               14, 21, # Mesenchyme 1-2 
                                                                                                               45, # Mesothelia
                                                                                                               27, 44, # B cell 1-2
                                                                                                               28,41,43, # Plasma cell 1-3
                                                                                                               10, # CD4+ T cell
                                                                                                               13, # CD8+ T cell
                                                                                                               17, # NK cell 1
                                                                                                               29, 46, # Mast cell 1-2
                                                                                                               5, # Monocyte
                                                                                                               11,9,20,34,40, # Macrophage 1-5
                                                                                                               35, # IFN-primed Macrophage
                                                                                                               1,6,15,23, # Alveolar Macrophage 1-4
                                                                                                               12, # cDC2
                                                                                                               39, # pDC
                                                                                                               18, # Proliferating
                                                                                                               31,32,37 # Doublets
                                                                                                               )])

### cell lineage annotation ###
cell.data <- data.table(barcode = colnames(lung.harmony),
                        celltype = Idents(lung.harmony))

lineage.annotation <- c("MP", "Epithelia", "Epithelia", "Epithelia", "MP", 
                        "MP", "Endothelia", "Epithelia", "MP", "T cell",
                        "MP", "MP", "T cell", "Mesenchyme", "MP",
                        "Epithelia", "NK cell", "Proliferating", "Epithelia", "MP",
                        "Mesenchyme", "Endothelia", "MP", "Epithelia", "Endothelia",
                        "Endothelia", "B cell", "Plasma cell", "Mast cell", "Epithelia",
                        "Doublet", "Doublet", "Epithelia", "MP", "MP", 
                        "Epithelia", "Doublet", "Epithelia", "MP", "MP", 
                        "Plasma cell", "Endothelia", "Plasma cell", "B cell", "Mesenchyme",
                        "Mast cell", "Epithelia")

lineage.data <- data.table(celltype = cluster.annotation, lineage = lineage.annotation)
meta.data <- merge(cell.data, lineage.data, by = "celltype")
meta.data <- data.frame(meta.data, row.names = meta.data$barcode)
meta.data$barcode <- NULL
meta.data$celltype <- NULL
lung.harmony <- AddMetaData(lung.harmony, meta.data, col.name = "lineage_integrated")
# order lineage factors
lung.harmony$lineage_integrated = factor(lung.harmony$lineage_integrated)
lung.harmony$lineage_integrated<- factor(lung.harmony$lineage_integrated, levels(lung.harmony$lineage_integrated)[c(4,3,6,1,9,11,8,5,7,10,2)])

# 

# check lineage annotation #
lineage.table1 <- table(lung.harmony$lineage, lung.harmony$lineage_integrated)
write.csv(lineage.table1, file = paste0(harmony.samples.path, "dim50_annotation/lineage_table_1.csv"))

lineage.table2 <- table(lung.harmony$condition, lung.harmony$lineage_integrated)
write.csv(lineage.table2, file = paste0(harmony.samples.path, "dim50_annotation/lineage_table_2.csv"))

lineage.table3 <- table(lung.harmony$patient.ID, lung.harmony$lineage_integrated)
write.csv(lineage.table3, file = paste0(harmony.samples.path, "dim50_annotation/lineage_table_3.csv"))

### save data ###
saveRDS(lung.harmony, "/home/s1987963/ds_group/Niklas/combined_organs/combined_lung_annotated.rds")

### subset MPs ###
lung.MPs <- subset(lung.harmony, subset = lineage_integrated == "MP") # 95936 cells

### save MPs ###
saveRDS(lung.harmony, "/home/s1987963/ds_group/Niklas/MP_lung/combined_lung_MPs_annotated.rds")
