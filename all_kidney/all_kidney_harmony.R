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
healthy.kidney.path <- "/home/s1987963/ds_group/Niklas/healthy_organs/healthy_kidney_annotated.rds"
fibrotic.kidney.path <- "/home/s1987963/ds_group/Niklas/fibrotic_organs/fibrotic_kidney_annotated.rds"
kidney.path <- "/home/s1987963/ds_group/Niklas/all_kidney/"
harmony.samples.path <- "/home/s1987963/ds_group/Niklas/all_kidney/harmonize_samples/HVGs_conditions/"

### load data ###
healthy.kidney <- readRDS(healthy.kidney.path) # 61059 cells
fibrotic.kidney <- readRDS(fibrotic.kidney.path) # 9945 cells

### add missing meta data
fibrotic.kidney$condition <- "fibrotic"

### merge data ###
kidney <- merge(healthy.kidney, fibrotic.kidney) # 71004 cells

### feature selection: HVGs ###
healthy.kidney <- FindVariableFeatures(healthy.kidney, selection.method = "vst", nfeatures = 4000)
fibrotic.kidney <- FindVariableFeatures(fibrotic.kidney, selection.method = "vst", nfeatures = 4000)
kidney.HVGs <- intersect(VariableFeatures(fibrotic.kidney), VariableFeatures(healthy.kidney)) # 2042
VariableFeatures(kidney) <- kidney.HVGs

### scale data ###
kidney <- ScaleData(kidney) # uncorrected

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

PC1_2.condition.plot <- DimPlot(object = kidney, reduction = "pca", pt.size = .1, group.by = "condition")
png(paste0(harmony.samples.path,"PC1_2.condition.png"), width=1000,height=1000,units="px")
print(PC1_2.condition.plot)
dev.off()

# save genes making up the PCs  
sink(paste0(harmony.samples.path, "PC_genes.txt"))
print(kidney[["pca"]], dims = 1:50, nfeatures = 20)
sink()

### integration with harmony ###
# harmonize samples
kidney.harmony <- kidney %>% RunHarmony("patient.ID", theta = 2, reduction.save = "harmony_theta2", plot_convergence = TRUE) # harmonize all 60 samples independently

## harmony elbow plot ##
harmony.elbow.plot <- ElbowPlot(kidney.harmony, ndims = 50, reduction = "harmony_theta2")
png(paste0(harmony.samples.path,"harmony_theta2.elbow.plot.png"), width=1000,height=600,units="px")
print(harmony.elbow.plot)
dev.off()

## explore harmony coordinates ##
## visualize PCs ##
harmony.PC1_2.samples.plot <- DimPlot(object = kidney.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "patient.ID")
png(paste0(harmony.samples.path,"harmony_theta2.PC1_2.samples.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.samples.plot)
dev.off()

harmony.PC1_2.condition.plot <- DimPlot(object = kidney.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "condition")
png(paste0(harmony.samples.path,"harmony_theta2.PC1_2.condition.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.condition.plot)
dev.off()

harmony.PC1_2.studies.plot <- DimPlot(object = kidney.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "study")
png(paste0(harmony.samples.path,"harmony_theta2.PC1_2.studies.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.studies.plot)
dev.off()

### save genes making up the PCs ### 
sink(paste0(harmony.samples.path, "harmony_PC_genes.txt"))
print(kidney.harmony[["harmony_theta2"]], dims = 1:50, nfeatures = 20)
sink()

### save data ###
saveRDS(kidney.harmony, paste0(harmony.samples.path, "all_kidney_harmony.rds"))

### explore different numbers of harmony PCs ###
# MNP genes
MNP.genes <- c("CD14", "FCGR3A","CSF1R", "CD68", "LYZ", 
               "CCR2", "CX3CR1", "S100A8", "MARCO", "MNDA", 
               "TIMD4", "CD163", "C1QB", "F13A1",	"MCEMP1",	
               "TREM2",	"CD9",	"VCAN")

MNP.kidney.genes <- c("CD14", "FCGR3A", "CSF2RA", "FCGR2A",
                      "C1QA", "RGS1", "CD52", "SEPP1")

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
dims <- c(50,47,41,33)
for(d in dims){
  
  # create folder
  dir.create(paste0(harmony.samples.path, "dim", d, "_annotation"))
  dim.path <- paste0(harmony.samples.path, "dim", d, "_annotation/")
  
  # run UMAP
  kidney.harmony <- RunUMAP(kidney.harmony, reduction = "harmony_theta2", dims = 1:d, seed.use=1)
  
  ## visualize batch effect 
  # no.1 overview
  sample.batch.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
  condition.batch.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = "condition", pt.size = 0.1)
  batch1.plot <- sample.batch.plot + condition.batch.plot
  png(paste0(dim.path, "UMAP_dim", d, ".batch.png"), width=1800, height=900, units="px")
  print(batch1.plot)
  dev.off()
  
  # no.2 - colored by patients
  batch2.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".patient.batch1.png"), width=1600, height=1000, units="px")
  print(batch2.plot)
  dev.off()
  
  # no.3 - split by patients
  batch3.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = "study", split.by = "patient.ID", pt.size = 0.01, ncol = 6)
  png(paste0(dim.path, "UMAP_dim", d, ".patient.batch2.png"), width=1800, height=600, units="px")
  print(batch3.plot)
  dev.off()
  
  # no.4 - colored by study
  batch4.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = "study", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".study.batch1.png"), width=1600, height=1000, units="px")
  print(batch4.plot)
  dev.off()
  
  # no.5 - split by study
  batch5.plot <- DimPlot(kidney.harmony, reduction = "umap", split.by = "study", pt.size = 0.01, ncol = 3)
  png(paste0(dim.path, "UMAP_dim", d, ".study.batch2.png"), width=1800, height=600, units="px")
  print(batch5.plot)
  dev.off()
  
  # no.6 - colored by condition
  batch6.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = "condition", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".condition.batch1.png"), width=1600, height=1000, units="px")
  print(batch6.plot)
  dev.off()
  
  # no.7 - split by condition
  batch7.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = "study", split.by = "condition", pt.size = 0.01, ncol = 2)
  png(paste0(dim.path, "UMAP_dim", d, ".condition.batch2.png"), width=1800, height=900, units="px")
  print(batch7.plot)
  dev.off()
  
  # no.8 - colored by lineage
  batch8.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = "lineage", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".lineage.batch1.png"), width=1600, height=1000, units="px")
  print(batch8.plot)
  dev.off()
  
  # no.9 - colored by annotation
  batch9.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = "ident", label=TRUE, pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".lineage.batch2.png"), width=1600, height=1000, units="px")
  print(batch9.plot)
  dev.off()
  
  # no.10 - lineage overview
  lineage.overview <- condition.batch.plot + batch8.plot
  png(paste0(dim.path, "UMAP_dim", d, ".lineage.overview.png"), width=2000, height=1000, units="px")
  print(lineage.overview)
  dev.off()
  
  ## compare healthy vs. fibrotic ##
  compare <- DimPlot(kidney.harmony, reduction = "umap", split.by = "condition", group.by = "condition", 
                     pt.size = 0.01) + NoLegend() + scale_color_viridis(discrete=TRUE)
  png(paste0(dim.path, "UMAP_dim", d, "healthy_fibrotic_1.png"), width=1800, height=900, units="px")
  print(compare)
  dev.off()
  
  ## compare original annotation ##
  annotation <- DimPlot(kidney.harmony, reduction = "umap", split.by = "condition", group.by = "celltype",
                        label = TRUE, repel = TRUE, pt.size = 0.01) + NoLegend()
  png(paste0(dim.path, "UMAP_dim", d, "healthy_fibrotic_2.png"), width=1800, height=900, units="px")
  print(annotation)
  dev.off()
  
  ## compare original annotation in healthy vs. fibrotic data ##
  plot <- compare / annotation
  png(paste0(dim.path,"UMAP_dim", d, "healthy_fibrotic_3.png"), width=2000, height=2000, units="px")
  print(plot)
  dev.off()
  
  ## QC metrics at cluster level ##
  cluster.qc.heatmap <- FeaturePlot(kidney.harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "cluster.qc.heatmap.png"), width=1500,height=500,units="px")
  print(cluster.qc.heatmap)
  dev.off()
  
  ## gene expression heatmaps ##
  # overview feature plot
  lineage.overview <- FeaturePlot(kidney.harmony, features = lineage.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "overview.markers.png"), width=1600,height=500,units="px")
  print(lineage.overview )
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
  
  # feature plot with kidney-specific MNP markers
  MNP.kidney.markers <- FeaturePlot(kidney.harmony, features = MNP.kidney.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "MNP.kidney.markers.png"), width=1600,height=800,units="px")
  print(MNP.kidney.markers)
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
  mesenchymal.markers <- FeaturePlot(kidney.harmony, features = mesenchymal.genes, pt.size = 0.5, ncol = 5) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"mesenchymal.markers.png"), width=1500,height=600,units="px")
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
  
  rm(cluster.qc.heatmap, lineage.overview, mesenchymal.markers, endothelial.markers, epithelial.markers, immunecell.markers,
     MNP.markers, cDC1.markers, cDC2.markers, pDC.markers, mastcell.markers, 
     Tcell.markers, NK.markers, Bcell.markers, plasmacell.markers,
     ery.markers, proliferation.markers)
}

### clustering ###
# only for selected number of dims (for reasons of computational efficiency!)
res <- seq(2.5,4.0,0.1)
for(d in dims){
  
  # set path variable
  dim.path <- paste0(harmony.samples.path, "dim", d, "_annotation/")
  
  # visualization with UMAP
  kidney.harmony <- RunUMAP(kidney.harmony, reduction = "harmony_theta2", dims=1:d, seed.use=1)
  
  # plot batch effect
  batch.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
  
  # create kNN graph
  kidney.harmony <- FindNeighbors(kidney.harmony, reduction = "harmony_theta2", dims = 1:d)
  
  for (r in res) {
    
    kidney.harmony <- FindClusters(kidney.harmony, reduction = "harmony_theta2", resolution = r)
    umap.plot <- DimPlot(kidney.harmony, reduction = "umap", label = TRUE, pt.size = 0.1)
    
    # create eval plot
    eval.plot <- umap.plot + batch.plot
    png(paste0(dim.path, "UMAP_dim", d, "_res", r, ".png"), width=1500, height=600, units="px")
    print(eval.plot)
    dev.off()
    
    ## clustering at patient level
    #patient.clustering <- DimPlot(kidney.harmony, group.by = "seurat_clusters", split.by = "patient.ID", pt.size = 0.1, ncol = 10)
    #png(paste0(dim.path,"UMAP_dim", d, "_res", r,".patient.clustering.png"), width=2000,height=400,units="px")
    #print(patient.clustering)
    #dev.off()
    
  }
}

### preliminary clustering ###
kidney.harmony <- FindNeighbors(kidney.harmony, reduction = "harmony_theta2", dims = 1:50)
kidney.harmony <- FindClusters(kidney.harmony, reduction = "harmony_theta2", resolution = 3.4)
# run UMAP
kidney.harmony <- RunUMAP(kidney.harmony, reduction = "harmony_theta2", dims=1:50, seed.use=1)

### save data ###
saveRDS(kidney.harmony, paste0(harmony.samples.path, "all_kidney_harmony.rds"))

# compare PTPRC expression across clusters
umap.plot <- DimPlot(kidney.harmony, reduction = "umap", label = T, label.size = 6, pt.size = 0.1)
ptprc.plot <- VlnPlot(object = kidney.harmony, features = c("PTPRC"), group.by = "seurat_clusters", pt.size = 0.1) + NoLegend()
immunecell.markers <- FeaturePlot(kidney.harmony, features = c("PTPRC"), pt.size = 0.5, ncol = 1) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
immuneclusters.plot <- umap.plot + immunecell.markers - ptprc.plot + plot_layout(ncol=1, widths=c(2,1))
png(paste0(harmony.samples.path, "dim50_annotation/immuneclusters.png"), width=1800,height=1200,units="px")
print(immuneclusters.plot)
dev.off()

## compute cluster marker genes ###
kidney.markers <- FindAllMarkers(kidney.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
kidney.top50.markers <- kidney.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
write.csv(kidney.markers, file = paste0(harmony.samples.path, "dim50_annotation/ALL_marker_genes.csv"))
write.csv(kidney.top50.markers, file = paste0(harmony.samples.path, "dim50_annotation/top50_marker_genes.csv"))

### cell type annotation ###
cluster.annotation <- c("Combined Kidney PCT 1", "Combined Kidney PCT 2", "Combined Kidney PCT 3",
                        "Combined Kidney PCT 4", "Combined Kidney Connecting Tubule 1", "Combined Kidney PCT 5", 
                        "Combined Kidney PCT 6", "Combined Kidney PCT 7", "Combined Kidney PCT 8", 
                        "Combined Kidney PCT 9", "Combined Kidney Distinct PCT 1", "Combined Kidney Thick Ascending LOH 1",
                        "Combined Kidney NK cell", "Combined Kidney T cell 1", "Combined Kidney Intercalated cell Type A 1",
                        "Combined Kidney Thick Ascending LOH 2", "Combined Kidney Principal cell 1", "Combined Kidney Endo 1",
                        "Combined Kidney PCT 10", "Combined Kidney DCT 1", "Combined Kidney Monocyte", 
                        "Combined Kidney Intercalated cell Type A 2", "Combined Kidney Intercalated cell","Combined Kidney DCT 2", 
                        "Combined Kidney PCT 11", "Combined Kidney Distinct PCT 2", "Combined Kidney Glomerular Endo", 
                        "Combined Kidney Thick Ascending LOH 3", "Combined Kidney Connecting Tubule 2", "Combined Kidney Epithelial Progenitor",
                        "Combined Kidney Thick Ascending LOH 4", "Combined Kidney Principal cell 2", "Combined Kidney PCT 12",
                        "Combined Kidney Intercalated cell Type B", "Combined Kidney Connecting Tubule 3", "Combined Kidney Podocyte",
                        "Combined Kidney Mesenchyme 1", "Combined Kidney Connecting Tubule 4", "Combined Kidney Principal cell 3",
                        "Combined Kidney Thick Ascending LOH 5", "Combined Kidney Capillary Endo", "Combined Kidney Macrophage",
                        "Combined Kidney PCT 13", "Combined Kidney DCT 3", "Combined Kidney PCT/collecting tubule doublet",
                        "Combined Kidney B cell", "Combined Kidney Pelvic Epithelium 1", "Combined Kidney T cell 2",
                        "Combined Kidney Mesenchyme 2", "Combined Kidney Pelvic Epithelium 2", "Combined Kidney Epithelial Progenitor/Principal cell/Mesenchyme doublet",
                        "Combined Kidney Intercalated cell Type A 3", "Combined Kidney PCT 14", "Combined Kidney vSMC",
                        "Combined Kidney Mesenchyme 3", "Combined Kidney PCT 15", "Combined Kidney Endo/Epithelial doublet"
)
names(cluster.annotation) <- levels(kidney.harmony)
kidney.harmony <- RenameIdents(kidney.harmony, cluster.annotation)

# add cell types to meta data
cell.data <- data.table(barcode = colnames(kidney.harmony),
                        celltype = Idents(kidney.harmony))
cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$barcode <- NULL
kidney.harmony <- AddMetaData(kidney.harmony, cell.data, col.name = "celltype_integrated")

# order cell type factors
kidney.harmony$celltype_integrated = factor(kidney.harmony$celltype_integrated)
kidney.harmony$celltype_integrated <- factor(kidney.harmony$celltype_integrated, levels(kidney.harmony$celltype_integrated)[c(1,2,3,4,6,7,8,9,10,19,25,33,43,53,56,# PCT
                                                                                                                              11,26, # Distinct PCT
                                                                                                                              # descending LOH
                                                                                                                              12,16,28,31,40, # ascending LOH
                                                                                                                              20,24,44, #DCT
                                                                                                                              5,29,35,38, #CT
                                                                                                                              15,22,52, # IC A
                                                                                                                              34, # IC B
                                                                                                                              23, # IC
                                                                                                                              17,32,39, # principal cells
                                                                                                                              47,50,# pelvic epithelium
                                                                                                                              36, # podocyte
                                                                                                                              30, # epithelial progenitor cell
                                                                                                                              18,27,41, # endothelia
                                                                                                                              37,49,55,54, # mesenchyme 
                                                                                                                              46, # B cell
                                                                                                                              13,14,48, # NKT cell
                                                                                                                              21,42, # MP
                                                                                                                              45,51,57 # doublets
                                                                                                                              )])

### cell lineage annotation ###
cell.data <- data.table(barcode = colnames(kidney.harmony),
                        celltype = Idents(kidney.harmony))

lineage.annotation <- c("Epithelia", "Epithelia", "Epithelia", "Epithelia", "Epithelia", "Epithelia",
                        "Epithelia", "Epithelia", "Epithelia", "Epithelia", "Epithelia", "Epithelia",
                        "NK cell", "T cell", "Epithelia", "Epithelia", "Epithelia", "Endothelia",
                        "Epithelia", "Epithelia", "MP", "Epithelia", "Epithelia", "Epithelia",
                        "Epithelia", "Epithelia", "Endothelia", "Epithelia", "Epithelia", "Epithelia",
                        "Epithelia", "Epithelia", "Epithelia", "Epithelia", "Epithelia", "Epithelia",
                        "Mesenchyme", "Epithelia", "Epithelia", "Epithelia", "Endothelia", "MP",
                        "Epithelia", "Epithelia", "Doublet", "B cell", "Epithelia", "T cell",
                        "Mesenchyme", "Epithelia", "Doublet", "Epithelia", "Epithelia", "Mesenchyme",
                        "Mesenchyme", "Epithelia", "Doublet"
)

lineage.data <- data.table(celltype = cluster.annotation, lineage = lineage.annotation)
meta.data <- merge(cell.data, lineage.data, by = "celltype")
meta.data <- data.frame(meta.data, row.names = meta.data$barcode)
meta.data$barcode <- NULL
meta.data$celltype <- NULL
kidney.harmony <- AddMetaData(kidney.harmony, meta.data, col.name = "lineage_integrated")

# order lineage factors
kidney.harmony$lineage_integrated  <- factor(kidney.harmony$lineage_integrated)
kidney.harmony$lineage_integrated<- factor(kidney.harmony$lineage_integrated, levels(kidney.harmony$lineage_integrated)[c(4,3,5,1,8,7,6,2)])

# check lineage annotation #
lineage.table1 <- table(kidney.harmony$lineage, kidney.harmony$lineage_integrated)
write.csv(lineage.table1, file = paste0(harmony.samples.path, "dim50_annotation/lineage_table_1.csv"))

lineage.table2 <- table(kidney.harmony$condition, kidney.harmony$lineage_integrated)
write.csv(lineage.table2, file = paste0(harmony.samples.path, "dim50_annotation/lineage_table_2.csv"))

lineage.table3 <- table(kidney.harmony$patient.ID, kidney.harmony$lineage_integrated)
write.csv(lineage.table3, file = paste0(harmony.samples.path, "dim50_annotation/lineage_table_3.csv"))

### save data ###
saveRDS(kidney.harmony, "/home/s1987963/ds_group/Niklas/combined_organs/combined_kidney_annotated.rds")

### subset MPs ###
kidney.MPs <- subset(kidney.harmony, subset = lineage_integrated == "MP") #  cells

### save MPs ###
saveRDS(kidney.MPs, "/home/s1987963/ds_group/Niklas/MP_kidney/combined_kidney_MPs_annotated.rds")
