library(Seurat)
library(SoupX)
library(celda)
library(umap)
library(reticulate)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
#library(Matrix)
#rowSums = Matrix::rowSums
#colSums = Matrix::colSums

### load data ###
lung.data.dir <- "/home/s1987963/processed_data/raredon_lung/healthy/GSM3926545_Hum1"
lung.data <- Read10X(lung.data.dir)

d <- as.matrix(lung.data)
lung.toc <- lung.tod

lung <- CreateSeuratObject(lung.tod, project = "lung", min.cells=3, min.features=200)
lung.counts <- GetAssayData(object = lung, assay = "RNA", slot = "data")
lung.counts <- as.matrix(lung.counts)

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

### SoupX ###
lung.sc <- SoupChannel(lung.tod, lung.toc, calcSoupProfile=TRUE, keepDroplets = TRUE)

# add clusters #
lung.sc <- setClusters(lung.sc, lung$seurat_clusters)

# estimating contamination fraction
lung.sc <- autoEstCont(lung.sc)

# manual
head(lung.sc$soupProfile[order(lung.sc$soupProfile$est, decreasing = TRUE), ], n = 50)







