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
healthy.liver.path <- "/home/s1987963/ds_group/Niklas/healthy_organs/healthy_liver_annotated.rds"
fibrotic.liver.path <- "/home/s1987963/ds_group/Niklas/fibrotic_organs/fibrotic_liver_annotated.rds"
liver.path <- "/home/s1987963/ds_group/Niklas/all_liver/"
harmony.samples.path <- "/home/s1987963/ds_group/Niklas/all_liver/harmonize_samples/HVGs_conditions/"

### load data ###
healthy.liver <- readRDS(healthy.liver.path) # 53337 cells
fibrotic.liver <- readRDS(fibrotic.liver.path) # 37071 cells

#### split healthy data ###
ramachandran.healthy <- subset(healthy.liver, subset = study == "ramachandran_liver")
ramachandran.healthy$dataset <- "ramachandran_healthy"
macparland.healthy <- subset(healthy.liver, subset = study == "macparland_liver")
macparland.healthy$dataset <- "macparland_healthy"
ramachandran.fibrotic <- subset(fibrotic.liver, subset = study == "ramachandran_liver")
ramachandran.fibrotic$dataset <- "ramachandran_fibrotic"

### merge data ###
healthy.liver <- merge(ramachandran.healthy, macparland.healthy)
fibrotic.liver <- ramachandran.fibrotic

### exclude 140 doublets from downstream analyses ###
healthy.liver <- subset(healthy.liver, subset = lineage != "Doublet")
fibrotic.liver <- subset(fibrotic.liver, subset = lineage != "Doublet")

#### normalization ###
#ramachandran.healthy <- NormalizeData(ramachandran.healthy, normalization.method = "LogNormalize", scale.factor = 10000)
#macparland.healthy <- NormalizeData(macparland.healthy, normalization.method = "LogNormalize", scale.factor = 10000)
#ramachandran.fibrotic <- NormalizeData(ramachandran.fibrotic, normalization.method = "LogNormalize", scale.factor = 10000)

### normalization ###
healthy.liver <- NormalizeData(healthy.liver, normalization.method = "LogNormalize", scale.factor = 10000)
fibrotic.liver <- NormalizeData(fibrotic.liver, normalization.method = "LogNormalize", scale.factor = 10000)

### merge data ###
liver <- merge(healthy.liver, fibrotic.liver) # 90268 cells

### feature selection: HVGs ###
healthy.liver <- FindVariableFeatures(healthy.liver, selection.method = "vst", nfeatures = 4000)
fibrotic.liver <- FindVariableFeatures(fibrotic.liver, selection.method = "vst", nfeatures = 4000)
liver.HVGs <- intersect(VariableFeatures(fibrotic.liver), VariableFeatures(healthy.liver)) # 2042
VariableFeatures(liver) <- liver.HVGs

### scale data ###
liver <- ScaleData(liver) # uncorrected

### dimensionality reduction: PCA ###
liver <- RunPCA(liver, features = VariableFeatures(object = liver))

# elbow plot
pca.elbow.plot <- ElbowPlot(liver, ndims = 50, reduction = "pca")
png(paste0(harmony.samples.path,"pca.elbow.plot.png"), width=1000,height=600,units="px")
print(pca.elbow.plot)
dev.off()

# visualize PCs
PC1_2.samples.plot <- DimPlot(object = liver, reduction = "pca", pt.size = .1, group.by = "patient.ID")
png(paste0(harmony.samples.path,"PC1_2.samples.png"), width=1000,height=1000,units="px")
print(PC1_2.samples.plot)
dev.off()

PC1_2.studies.plot <- DimPlot(object = liver, reduction = "pca", pt.size = .1, group.by = "dataset")
png(paste0(harmony.samples.path,"PC1_2.studies.png"), width=1000,height=1000,units="px")
print(PC1_2.studies.plot)
dev.off()

# save genes making up the PCs  
sink(paste0(harmony.samples.path, "PC_genes.txt"))
print(liver[["pca"]], dims = 1:50, nfeatures = 20)
sink()

### integration with harmony ###
# harmonize samples
liver.harmony <- liver %>% RunHarmony("patient.ID", theta = 2, reduction.save = "harmony_theta2", plot_convergence = TRUE) # harmonize all 12 samples independently

## harmony elbow plot ##
harmony.elbow.plot <- ElbowPlot(liver.harmony, ndims = 50, reduction = "harmony_theta2")
png(paste0(harmony.samples.path,"harmony_theta2.elbow.plot.png"), width=1000,height=600,units="px")
print(harmony.elbow.plot)
dev.off()

## explore harmony coordinates ##
## visualize PCs ##
harmony.PC1_2.samples.plot <- DimPlot(object = liver.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "patient.ID")
png(paste0(harmony.samples.path,"harmony_theta2.PC1_2.samples.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.samples.plot)
dev.off()

harmony.PC1_2.studies.plot <- DimPlot(object = liver.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "dataset")
png(paste0(harmony.samples.path,"harmony_theta2.PC1_2.studies.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.studies.plot)
dev.off()

### save genes making up the PCs ### 
sink(paste0(harmony.samples.path, "harmony_PC_genes.txt"))
print(liver.harmony[["harmony_theta2"]], dims = 1:50, nfeatures = 20)
sink()

### save data ###
saveRDS(liver.harmony, paste0(harmony.samples.path, "all_liver_harmony.rds"))

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
#dims <- c(32,40,50,19) # harmonize samples
dims <- c(30,38,44,50)
for(d in dims){
  
  # create folder
  dir.create(paste0(harmony.samples.path, "dim", d, "_annotation"))
  dim.path <- paste0(harmony.samples.path, "dim", d, "_annotation/")
  
  # run UMAP
  liver.harmony <- RunUMAP(liver.harmony, reduction = "harmony_theta2", dims = 1:d, seed.use=1)
  
  ## visualize batch effect 
  # no.1 overview
  sample.batch.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
  dataset.batch.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "dataset", pt.size = 0.1)
  condition.batch.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "condition", pt.size = 0.1)
  batch1.plot <- sample.batch.plot + dataset.batch.plot + condition.batch.plot
  png(paste0(dim.path, "UMAP_dim", d, ".batch.png"), width=1800, height=600, units="px")
  print(batch1.plot)
  dev.off()
  
  # no.2 - colored by patients
  batch2.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".patient.batch1.png"), width=1600, height=1000, units="px")
  print(batch2.plot)
  dev.off()
  
  # no.3 - split by patients
  batch3.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "study", split.by = "patient.ID", pt.size = 0.01, ncol = 6)
  png(paste0(dim.path, "UMAP_dim", d, ".patient.batch2.png"), width=1800, height=600, units="px")
  print(batch3.plot)
  dev.off()
  
  # no.4 - colored by dataset
  batch4.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "dataset", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".study.batch1.png"), width=1600, height=1000, units="px")
  print(batch4.plot)
  dev.off()
  
  # no.5 - split by dataset
  batch5.plot <- DimPlot(liver.harmony, reduction = "umap", split.by = "dataset", pt.size = 0.01, ncol = 3)
  png(paste0(dim.path, "UMAP_dim", d, ".study.batch2.png"), width=1800, height=600, units="px")
  print(batch5.plot)
  dev.off()
  
  # no.6 - colored by condition
  batch6.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "condition", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".condition.batch1.png"), width=1600, height=1000, units="px")
  print(batch6.plot)
  dev.off()
  
  # no.7 - split by condition
  batch7.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "study", split.by = "condition", pt.size = 0.01, ncol = 2)
  png(paste0(dim.path, "UMAP_dim", d, ".condition.batch2.png"), width=1800, height=900, units="px")
  print(batch7.plot)
  dev.off()
  
  # no.8 - colored by lineage
  batch8.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "lineage", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".lineage.batch1.png"), width=1600, height=1000, units="px")
  print(batch8.plot)
  dev.off()
  
  # no.9 - colored by annotation
  batch9.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "ident", label=TRUE, pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".lineage.batch2.png"), width=1600, height=1000, units="px")
  print(batch9.plot)
  dev.off()
  
  # no.10 - lineage overview
  lineage.overview <- condition.batch.plot + batch8.plot
  png(paste0(dim.path, "UMAP_dim", d, ".lineage.overview.png"), width=2000, height=1000, units="px")
  print(lineage.overview)
  dev.off()
  
  ## compare healthy vs. fibrotic ##
  compare <- DimPlot(liver.harmony, reduction = "umap", split.by = "condition", group.by = "condition", 
                     pt.size = 0.01) + NoLegend() + scale_color_viridis(discrete=TRUE)
  png(paste0(dim.path, "UMAP_dim", d, "healthy_fibrotic_1.png"), width=1800, height=900, units="px")
  print(compare)
  dev.off()
  
  ## compare original annotation ##
  annotation <- DimPlot(liver.harmony, reduction = "umap", split.by = "condition", group.by = "celltype",
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
  cluster.qc.heatmap <- FeaturePlot(liver.harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "cluster.qc.heatmap.png"), width=1500,height=500,units="px")
  print(cluster.qc.heatmap)
  dev.off()
  
  # feature plot with immune cell marker
  immunecell.markers <- FeaturePlot(liver.harmony, features = c("PTPRC"), pt.size = 0.5, ncol = 1) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "immunecell.markers.png"), width=1000,height=1000,units="px")
  print(immunecell.markers)
  dev.off()
  
  # feature plot with MNP markers
  MNP.markers <- FeaturePlot(liver.harmony, features = MNP.genes, pt.size = 0.5, ncol = 6) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "MNP.markers.png"), width=1800,height=900,units="px")
  print(MNP.markers)
  dev.off()
  
  # feature plot with cDC 1 markers
  cDC1.markers <- FeaturePlot(liver.harmony, features = cDC1.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "cDC1.markers.png"), width=1600,height=400,units="px")
  print(cDC1.markers)
  dev.off()
  
  # feature plot with cDC 2 markers
  cDC2.markers <- FeaturePlot(liver.harmony, features = cDC2.genes, pt.size = 0.5, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path, "UMAP_dim", d, "cDC2.markers.png"), width=1200,height=800,units="px")
  print(cDC2.markers)
  dev.off()
  
  # feature plot with pDC markers
  pDC.markers <- FeaturePlot(liver.harmony, features = pDC.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "pDC.markers.png"), width=1600,height=400,units="px")
  print(pDC.markers)
  dev.off()
  
  # feature plot with T cell markers
  Tcell.markers <- FeaturePlot(liver.harmony, features = Tcell.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "Tcell.markers.png"), width=1600,height=800,units="px")
  print(Tcell.markers)
  dev.off()
  
  # feature plot with NK cell markers
  NK.markers <- FeaturePlot(liver.harmony, features = NK.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "NK.markers.png"), width=1600,height=800,units="px")
  print(NK.markers)
  dev.off()
  
  # feature plot with B cell markers
  Bcell.markers <- FeaturePlot(liver.harmony, features = Bcell.genes, pt.size = 0.5, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "Bcell.markers.png"), width=1200,height=400,units="px")
  print(Bcell.markers)
  dev.off()
  
  # feature plot with plasma cell markers
  plasmacell.markers <- FeaturePlot(liver.harmony, features = plasmacell.genes, pt.size = 0.5, ncol = 5) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "plasmacell.markers.png"), width=1500,height=300,units="px")
  print(plasmacell.markers)
  dev.off()
  
  # feature plot with mast cell markers
  mastcell.markers <- FeaturePlot(liver.harmony, features = mastcell.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"mastcell.markers.png"), width=1600,height=400,units="px")
  print(mastcell.markers)
  dev.off()
  
  # feature plot with epithelial cell markers
  epithelial.markers <- FeaturePlot(liver.harmony, features = epithelial.genes, pt.size = 0.5, ncol = 5) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"epithelial.markers.png"), width=1500,height=400,units="px")
  print(epithelial.markers)
  dev.off()
  
  # feature plot with endothelial cell markers
  endothelial.markers <- FeaturePlot(liver.harmony, features = endothelial.genes, pt.size = 0.5, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"endothelial.markers.png"), width=1200,height=1200,units="px")
  print(endothelial.markers)
  dev.off()
  
  # feature plot with mesenchymal cell markers
  mesenchymal.markers <- FeaturePlot(liver.harmony, features = mesenchymal.genes, pt.size = 0.5, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"mesenchymal.markers.png"), width=1200,height=1200,units="px")
  print(mesenchymal.markers)
  dev.off()
  
  # feature plot with erythroid cell markers
  ery.markers <- FeaturePlot(liver.harmony, features = ery.genes, pt.size = 0.5, ncol = 2) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"ery.markers.png"), width=800,height=400,units="px")
  print(ery.markers)
  dev.off()
  
  # proliferating cells
  proliferation.markers <- FeaturePlot(liver.harmony, features = proliferation.genes, pt.size = 0.5, ncol = 2) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"proliferation.markers.png"), width=800,height=400,units="px")
  print(proliferation.markers)
  dev.off()
  
}

### clustering ###
res <- seq(0.5,2.0,0.1)
for(d in dims){
  
  # set path variable
  dim.path <- paste0(harmony.samples.path, "dim", d, "_annotation/")
  
  # calculate UMAP
  liver.harmony <- RunUMAP(liver.harmony, reduction = "harmony_theta2", dims=1:d, seed.use=1)
  
  # construct kNN graph
  liver.harmony <- FindNeighbors(liver.harmony, reduction = "harmony_theta2", dims = 1:d)
  
  # plot batch effect
  batch.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "condition", pt.size = 0.01)
  
  for (r in res) {
    # compute clusters
    liver.harmony <- FindClusters(liver.harmony, reduction = "harmony_theta2", resolution = r)
    
    # plot clusters
    cluster.plot <- DimPlot(liver.harmony, reduction = "umap", label = TRUE, pt.size = 0.1)
    
    # create eval plot
    eval.plot <- cluster.plot + batch.plot
    png(paste0(dim.path, "UMAP_dim", d, "_res", r, ".png"), width=1500, height=600, units="px")
    print(eval.plot)
    dev.off()
  }
}

### fibrotic vs healthy comparison ###
compare <- DimPlot(liver.harmony, reduction = "umap", split.by = "condition", group.by = "condition", 
                   pt.size = 0.01) + NoLegend() + scale_color_viridis(discrete=TRUE)
png(paste0(harmony.samples.path, "dim50_annotation/figure1.png"), width=1800, height=900, units="px")
print(compare)
dev.off()

annotation <- DimPlot(liver.harmony, reduction = "umap", split.by = "condition", group.by = "celltype",
                      label = TRUE, repel = TRUE, pt.size = 0.01) + NoLegend()
png(paste0(harmony.samples.path, "dim50_annotation/figure2.png"), width=1800, height=900, units="px")
print(annotation)
dev.off()

plot <- compare / annotation
png(paste0(harmony.samples.path, "dim50_annotation/figure.png"), width=2000, height=2000, units="px")
print(plot)
dev.off()

### preliminary clustering ###
liver.harmony <- FindNeighbors(liver.harmony, reduction = "harmony_theta2", dims = 1:50)
liver.harmony <- FindClusters(liver.harmony, reduction = "harmony_theta2", resolution = 2.0)
# run UMAP
liver.harmony <- RunUMAP(liver.harmony, reduction = "harmony_theta2", dims=1:50, seed.use=1)

### save data ###
saveRDS(liver.harmony, paste0(harmony.samples.path, "all_liver_harmony.rds"))

# compare PTPRC expression across clusters
umap.plot <- DimPlot(liver.harmony, reduction = "umap", label = T, label.size = 6, pt.size = 0.1)
ptprc.plot <- VlnPlot(object = liver.harmony, features = c("PTPRC"), group.by = "seurat_clusters", pt.size = 0.1) + NoLegend()
immuneclusters.plot <- umap.plot + immunecell.markers - ptprc.plot + plot_layout(ncol=1, widths=c(2,1))
png(paste0(harmony.samples.path, "dim50_annotation/immuneclusters.png"), width=1800,height=1200,units="px")
print(immuneclusters.plot)
dev.off()

## compute cluster marker genes ###
liver.markers <- FindAllMarkers(liver.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
liver.top50.markers <- liver.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
write.csv(liver.markers, file = paste0(harmony.samples.path, "dim50_annotation/ALL_marker_genes.csv"))
write.csv(liver.top50.markers, file = paste0(harmony.samples.path, "dim50_annotation/top50_marker_genes.csv"))

### cell type annotation ###
cluster.annotation <- c("Combined Liver CD8+ T cell 1", "Combined Liver NK1", "Combined Liver CD4+ T cell 1",
                        "Combined Liver CD4+ T cell 2", "Combined Liver KC1", "Combined Liver NK2",
                        "Combined Liver Hepatocyte 1", "Combined Liver LSEndo 1", "Combined Liver Macrophage 1",
                        "Combined Liver Cholangiocyte 1", "Combined Liver Endo 1", "Combined Liver Monocyte 1",
                        "Combined Liver PV Endo", "Combined Liver B cell", "Combined Liver Endo 2",
                        "Combined Liver Mesenchyme 1", "Combined Liver HA Endo", "Combined Liver CV Endo",
                        "Combined Liver cDC2", "Combined Liver Monocyte 2", "Combined Liver Proliferating 1",
                        "Combined Liver Plasma cell 1", "Combined Liver NK Cell 3", "Combined Liver LSEndo 2",
                        "Combined Liver Mesenchyme 2", "Combined Liver CTLA4+ T cell", "Combined Liver LymphEndo",
                        "Combined Liver KC2", "Combined Liver pDC", "Combined Liver cDC1",
                        "Combined Liver Hepatocyte 2", "Combined Liver KC3", "Combined Liver Cholangiocyte 2",
                        "Combined Liver Proliferating 3", "Combined Liver Basophil", "Combined Liver Hepatocyte 3",
                        "Combined Liver Mast Cell", "Combined Liver Proliferating 4", "Combined Liver LSEndo 3",
                        "Combined Liver Cholangiocyte 3", "Combined Liver IFN primed T cell", "Combined Liver Hepatocyte 4",
                        "Combined Liver Red blood cell", "Combined Liver Cholangiocyte 4", "Combined Liver NK Cell 4"
)
names(cluster.annotation) <- levels(liver.harmony)
liver.harmony <- RenameIdents(liver.harmony, cluster.annotation)
# add cell types to meta data
cell.data <- data.table(barcode = colnames(liver.harmony),
                        celltype = Idents(liver.harmony))
cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$barcode <- NULL
liver.harmony <- AddMetaData(liver.harmony, cell.data, col.name = "celltype_integrated")

# save annotated UMAP
annotated.umap.plot <- DimPlot(liver.harmony, reduction = "umap", label = T, repel = TRUE, label.size = 5, pt.size = 0.1)
png(paste0(harmony.samples.path, "dim50_annotation/UMAP_annotated.png"), width=1800,height=1200,units="px")
print(annotated.umap.plot)
dev.off()

### cell lineage annotation ###
cell.data <- data.table(barcode = colnames(liver.harmony),
                        celltype = Idents(liver.harmony))

lineage.annotation <- c("T cell","NK cell","T cell","T cell","MP",
                        "NK cell","Epithelia","Endothelia","MP","Epithelia",
                        "Endothelia","MP","Endothelia","B cell","Endothelia",
                        "Mesenchyme","Endothelia","Endothelia","MP","MP",
                        "Proliferating","Plasma cell", "NK cell","Endothelia","Mesenchyme",
                        "T cell","Endothelia","MP","MP","MP",
                        "Epithelia","MP", "Epithelia","Proliferating","Basophil",
                        "Epithelia","Mast cell","Proliferating","Endothelia","Epithelia",
                        "T cell","Epithelia","Red blood cell","Epithelia","NK cell"
)

lineage.data <- data.table(celltype = cluster.annotation, lineage = lineage.annotation)
meta.data <- merge(cell.data, lineage.data, by = "celltype")
meta.data <- data.frame(meta.data, row.names = meta.data$barcode)
meta.data$barcode <- NULL
meta.data$celltype <- NULL
liver.harmony <- AddMetaData(liver.harmony, meta.data, col.name = "lineage_integrated")

# save annotated UMAP
annotated.umap.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "lineage_integrated", label = T, repel = TRUE, label.size = 10, pt.size = 0.1)
png(paste0(harmony.samples.path, "dim50_annotation/UMAP_annotated.lineage.png"), width=1800,height=1200,units="px")
print(annotated.umap.plot)
dev.off()

# check lineage annotation #
lineage.table1 <- table(liver.harmony$lineage, liver.harmony$lineage_integrated)
write.csv(lineage.table1, file = paste0(harmony.samples.path, "dim50_annotation/lineage_table_1.csv"))

lineage.table2 <- table(liver.harmony$condition, liver.harmony$lineage_integrated)
write.csv(lineage.table2, file = paste0(harmony.samples.path, "dim50_annotation/lineage_table_2.csv"))

lineage.table3 <- table(liver.harmony$patient.ID, liver.harmony$lineage_integrated)
write.csv(lineage.table3, file = paste0(harmony.samples.path, "dim50_annotation/lineage_table_3.csv"))

### save R session ###
save.image(file = "/home/s1987963/ds_group/Niklas/all_liver/harmonize_samples/HVGs_conditions/combined_liver_harmony.RData")

### save data ###
saveRDS(liver.harmony, "/home/s1987963/ds_group/Niklas/combined_organs/combined_liver_annotated.rds")

