library(reticulate)
use_condaenv("r-reticulate", required = TRUE)
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

### path variables ###
lung.path <- "/Users/Niklas/Documents/Bioinformatics/MacrophageAtlas/reyfman_filtered/"

### prepare data for scrublet ###
lung.data <- Read10X_h5(paste0(lung.path, "fibrotic_05_filtered.h5"))
lung <- CreateSeuratObject(counts = lung.data, project = "reyfman_lung", min.cells = 3, min.features = 200)

### normalisation ###
lung <- NormalizeData(lung, normalization.method = "LogNormalize", scale.factor = 10000)

### feature selection ###
lung <- FindVariableFeatures(lung, selection.method = "vst", nfeatures = 2000)

### scale data ###
lung <- ScaleData(lung)

### dimensionality reduction: PCA ###
lung <- RunPCA(lung, features = VariableFeatures(object = lung))
ElbowPlot(lung, ndims = 50)

### clustering ###
dims <- c(10,14)
res <- seq(0.5,1.0,0.1)

for(d in dims){
  lung <- RunUMAP(lung, dims=1:d, seed.use=1)
  for (r in res) {
    lung <- FindNeighbors(lung, dims = 1:d)
    lung <- FindClusters(lung, resolution = r)
    umap.plot <- DimPlot(lung, reduction = "umap", label = F, pt.size = 0.1)
    #batch.plot <- DimPlot(lung, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
    #eval.plot <- umap.plot + batch.plot
    png(paste0(lung.path, "fibrotic5_UMAP_dim", d, "_res", r, ".png"), width=1000, height=800, units="px")
    print(umap.plot)
    dev.off()
  }
}

### definite clustering ###
lung <- RunUMAP(lung, dims=1:10, seed.use = 1 )
lung <- FindNeighbors(lung, dims = 1:10)
lung <- FindClusters(lung, resolution = 0.5)

### tSNE ###
lung <- RunTSNE(lung, dims = 1:10, seed.use = 1)

### save for scrublet ###
save(lung, file = paste0(lung.path, "fibrotic5_seurat.RData"))

### load packages ###
scr <- import("scrublet")
plt <- import("matplotlib.pyplot")

### load data ###
load(file= paste0(lung.path, "fibrotic5_seurat.RData"))

### remove doublets ###
scrubbed <- scr$Scrublet(t(as.matrix(lung@assays$RNA@data)), 
                         expected_doublet_rate=0.1,
                         sim_doublet_ratio=5)

tmp <- scrubbed$scrub_doublets()

scrubbed$plot_histogram()
plt$savefig(paste0(lung.path,"hist_doublet_threshold.png"))

lung$scrublet_auto <- tmp[[2]]
TSNEPlot(lung, group.by="scrublet_auto") 
#dev.print(png, paste0(lung.path,"/tsne_doublet_autothreshold.png"), width=5,height=5,res=300,units="in")
DimPlot(lung, group.by="scrublet_auto")
#dev.print(png, paste0(lung.path,"/umap_doublet_autothreshold.png"), width=5,height=5,res=300,units="in")

lung$scrublet_manual <- tmp[[1]] > 0.15
TSNEPlot(lung, group.by="scrublet_manual") + NoLegend()
#dev.print(png, paste0(lung.path,"/tsne_doublet_manualthreshold.png"), width=5,height=5,res=300,units="in")
DimPlot(lung, group.by="scrublet_manual") 
#dev.print(png, paste0(lung.path,"/umap_doublet_manualthreshold.png"), width=5,height=5,res=300,units="in")
