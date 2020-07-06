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
wilson.kidney.path <- "/home/s1987963/ds_group/Niklas/wilson_kidney/diabetes/wilson_kidney_diabetes_filtered.rds"
harmony.samples.path <- "/home/s1987963/ds_group/Niklas/fibrotic_kidney/harmonize_samples/"

### read files ###
kidney <- readRDS(wilson.kidney.path) # 9945 cells

### normalization ###
kidney <- NormalizeData(kidney, normalization.method = "LogNormalize", scale.factor = 10000)

### feature selection ###
## approach no.1: select top 2000 overall HVGs ##
kidney <- FindVariableFeatures(kidney, selection.method = "vst", nfeatures = 2000)

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

# save genes making up the PCs  
sink(paste0(harmony.samples.path, "PC_genes.txt"))
print(kidney[["pca"]], dims = 1:50, nfeatures = 20)
sink()

### integration with harmony ###
# harmonize samples
kidney.harmony <- kidney %>% RunHarmony("patient.ID", theta = 2, reduction.save = "harmony_theta2", plot_convergence = TRUE) # harmonize all 3 samples

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

### save genes making up the PCs ### 
sink(paste0(harmony.samples.path, "harmony_PC_genes.txt"))
print(kidney.harmony[["harmony_theta2"]], dims = 1:50, nfeatures = 20)
sink()

### save data ###
saveRDS(kidney, paste0(harmony.samples.path, "fibrotic_kidney_harmony.rds"))

### marker genes ###
# overview markers
lineage.genes <- c("PTPRC", "EPCAM", "PECAM1", "PDGFRA")

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
dims <- c(22,30,37,50)
for(d in dims){
  
  # create folder
  dir.create(paste0(harmony.samples.path, "dim", d, "_annotation"))
  dim.path <- paste0(harmony.samples.path, "dim", d, "_annotation/")
  
  # run UMAP
  kidney.harmony <- RunUMAP(kidney.harmony, reduction = "harmony_theta2", dims = 1:d, seed.use=1)
  
  ## batch plots ##
  # no.1 - colored by patients
  batch1.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".patient.batch1.png"), width=1600, height=1000, units="px")
  print(batch1.plot)
  dev.off()
  
  # no.2 - split by patients
  batch2.plot <- DimPlot(kidney.harmony, reduction = "umap", group.by = "study", split.by = "patient.ID", pt.size = 0.01, ncol = 4)
  png(paste0(dim.path, "UMAP_dim", d, ".patient.batch2.png"), width=1600, height=400, units="px")
  print(batch2.plot)
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
  
  # feature plot with kidney specific MNP markers
  kidneyMNP.markers <- FeaturePlot(kidney.harmony, features = MNP.kidney.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "kidneyMNP.markers.png"), width=1600,height=800,units="px")
  print(kidneyMNP.markers)
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
  png(paste0(dim.path,"UMAP_dim", d,"epithelial.markers.png"), width=1500,height=400,units="px")
  print(epithelial.markers)
  dev.off()
  
  # feature plot with endothelial cell markers
  endothelial.markers <- FeaturePlot(kidney.harmony, features = endothelial.genes, pt.size = 0.5, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"endothelial.markers.png"), width=1200,height=1200,units="px")
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
  
  rm(cluster.qc.heatmap, lineage.overview, mesenchymal.markers, endothelial.markers, epithelial.markers, immunecell.markers,
     MNP.markers, cDC1.markers, cDC2.markers, pDC.markers, mastcell.markers, 
     Tcell.markers, NK.markers, Bcell.markers, plasmacell.markers,
     ery.markers, proliferation.markers)
}

### clustering ###
# only for selected number of dims (for reasons of computational efficiency!)
dims <- c(30,37,50)
res <- seq(0.5,2.5,0.1)
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
    umap.plot <- DimPlot(kidney.harmony, reduction = "umap", label = F, pt.size = 0.1)
    
    # create eval plot
    eval.plot <- umap.plot + batch.plot
    png(paste0(dim.path, "UMAP_dim", d, "_res", r, ".png"), width=1500, height=600, units="px")
    print(eval.plot)
    dev.off()
    
    # clustering at patient level
    patient.clustering <- DimPlot(kidney.harmony, group.by = "seurat_clusters", split.by = "patient.ID", pt.size = 0.1, ncol = 5)
    png(paste0(dim.path,"UMAP_dim", d, "_res", r,".patient.clustering.png"), width=2000,height=400,units="px")
    print(patient.clustering)
    dev.off()
    
  }
}
