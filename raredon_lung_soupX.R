library(Seurat)
library(SoupX)
library(rhdf5)
library(umap)
library(reticulate)
library(dplyr)
library(ggplot2)
library(patchwork)


### load data ###
lung.data.dir <- "/home/s1987963/processed_data/raredon_lung/healthy/GSM3926545_Hum1"
lung.data <- Read10X(lung.data.dir)

lung <- CreateSeuratObject(lung.data, project = "lung", min.cells=3, min.features=200)

### SoupX ###
lung.sc <- SoupChannel(lung.data, lung[["RNA"]])
