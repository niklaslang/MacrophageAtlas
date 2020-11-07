library(Seurat)
library(scater)
library(SingleCellExperiment)
library(data.table)
library(edgeR)
library(reticulate)
library(umap)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)
library(pheatmap)
library(scales)
library(RColorBrewer)
library(viridis)
library(patchwork)
library(EnhancedVolcano)

### path variables ###
lung.path <- "/home/s1987963/ds_group/Niklas//MP_lung/marker_genes/lung_MP_wilcoxon_marker.csv"
liver.path <- "/home/s1987963/ds_group/Niklas//MP_liver/marker_genes/liver_MP_wilcoxon_marker.csv"
kidney.path <- "/home/s1987963/ds_group/Niklas//MP_kidney/marker_genes/kidney_MP_wilcoxon_marker.csv"
blood.path <- "/home/s1987963/ds_group/Niklas//MP_blood/marker_genes/blood_MP_wilcoxon_marker.csv"

### read marker genes ###
lung.marker <- fread(lung.path, stringsAsFactors = TRUE)
liver.marker <- fread(liver.path, stringsAsFactors = TRUE)
kidney.marker <- fread(kidney.path, stringsAsFactors = TRUE)
blood.marker <- fread(blood.path, stringsAsFactors = TRUE)

### macrophage1 subset marker ###
lung.mac1 <- as.vector(lung.marker[lung.marker[, cluster == "Combined Lung Macrophage 1"]]$gene)
liver.mac1 <- as.vector(liver.marker[liver.marker[, cluster == "Combined Liver Macrophage 1"]]$gene)
kidney.mac1 <- as.vector(kidney.marker[kidney.marker[, cluster == "Combined Kidney Macrophage"]]$gene)

### identify conserved marker of macrophage 1 ###
cons.M1.markers <- intersect(lung.mac1, liver.mac1)
cons.markers2 <- intersect(cons.M1.markers, kidney.mac1)

### macrophage1 subset marker ###
lung.AM1 <- as.vector(lung.marker[lung.marker[, cluster == "Combined Lung Alveolar Macrophage 1"]]$gene)
liver.KC1 <- as.vector(liver.marker[liver.marker[, cluster == "Combined Liver KC1"]]$gene)

### identify conserved marker of macrophage 1 ###
cons.TM1.markers <- intersect(lung.AM1, liver.KC1)

### macrophage3 subset marker ###
lung.AM3 <- as.vector(lung.marker[lung.marker[, cluster == "Combined Lung Alveolar Macrophage 3"]]$gene)
liver.KC3 <- as.vector(liver.marker[liver.marker[, cluster == "Combined Liver KC3"]]$gene)

### identify conserved marker of macrophage 1 ###
cons.TM3.markers <- intersect(lung.AM3, liver.KC3)

### unique TM 1 and TM 3 markers ###
TM1.diff.markers <- setdiff(cons.TM1.markers, cons.TM3.markers)
TM3.diff.markers <- setdiff(cons.TM3.markers, cons.TM1.markers)
unique.TM1 <- cons.TM1.markers[which(cons.TM1.markers %in% TM1.diff.markers)]
unique.TM3 <- cons.TM3.markers[which(cons.TM3.markers %in% TM3.diff.markers)]
