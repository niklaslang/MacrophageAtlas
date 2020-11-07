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
fig.path <- "/home/s1987963/ds_group/Niklas/mesenchyme_liver/harmonize_samples/figures/"
dir.create(fig.path)

### read liver mesenchyme data ###
liver.mesenchyme <- readRDS("/home/s1987963/ds_group/Niklas/combined_mesenchyme/combined_liver_mesenchyme_subclustered_annotated.rds")

### remove doublets ###
liver.mesenchyme <- subset(liver.mesenchyme, subset = lineage_integrated != "Doublet")

### figure 1: annotation ###
fig1 <- DimPlot(liver.mesenchyme, reduction = "umap", label = TRUE, label.size = 4, repel = TRUE, pt.size = 0.5) + NoLegend() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
png(paste0(fig.path, "figure_1.png"), width=1000, height=1000, units="px")
print(fig1)
dev.off()