library(Seurat)
library(SingleCellExperiment)
library(celda)
library(umap)
library(reticulate)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

### load data ###
lung.data.dir <- "smb://cmvm.datastore.ed.ac.uk/cmvm/scs/groups/pramacha-GROUP/Niklas/decontX/GSM3926545_Hum1"
lung.data <- Read10X(lung.data.dir)

### create non-sparse matrix for decontX ###
lung.counts <- as.matrix(lung.data)
storage.mode(lung.counts) <- "integer"

### intialize seurat object ###
lung <- CreateSeuratObject(lung.data, project = "lung", min.cells=3, min.features=200)

### process lung data ###
### normalisation ###
lung <- NormalizeData(lung, normalization.method = "LogNormalize", scale.factor = 10000)
### feature selection ###
lung <- FindVariableFeatures(lung, selection.method = "vst", nfeatures = 2000)
### scale data ###
lung <- ScaleData(lung, vars.to.regress = "nFeature_RNA")
### dimensionality reduction: PCA ###
lung <- RunPCA(lung, features = VariableFeatures(object = lung))
### UMAP ###
lung <- RunUMAP(lung, dims=1:10, seed.use=1)
### clustering ###
lung <- FindNeighbors(lung, dims = 1:10)
lung <- FindClusters(lung, resolution = 0.9)

### decontX ###
lung.counts.contamination <- decontX(counts = lung.counts, z = NULL)