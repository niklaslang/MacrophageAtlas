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
MP.path <- "/home/s1987963/ds_group/Niklas/combined_MPs/combined_lung_MPs_annotated.rds"
marker.path <- "/home/s1987963/ds_group/Niklas/MP_lung/marker_genes/"

### load data ###
MP <- readRDS(MP.path) # 95936 cells

### find marker genes using wilcoxon rank sum test ###
MP.wilcoxon.marker <- FindAllMarkers(MP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
write.csv(MP.wilcoxon.marker, file = paste0(marker.path, "lung_MP_wilcoxon_marker.csv"))

### find marker genes using MAST statistical framework ###
MP.MAST.marker <- FindAllMarkers(MP, only.pos = TRUE, test.use = "MAST", min.pct = 0.25, logfc.threshold = 0.25) 
write.csv(MP.MAST.marker, file = paste0(marker.path, "lung_MP_MAST_marker.csv"))