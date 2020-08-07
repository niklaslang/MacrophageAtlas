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
fig.path <- "/home/s1987963/ds_group/Niklas/all_lung/harmonize_samples/figures/"

### read lung data ###
lung <- readRDS("/home/s1987963/ds_group/Niklas/combined_organs/combined_lung_annotated.rds")

### remove doublets ###
lung <- subset(lung, subset = lineage_integrated != "Doublet")

### reorder condition factor ###
lung$study <- factor(lung$study)
lung$study <- factor(lung$study, levels(lung$study)[c(2,3,1)])

### figure 1: annotation ###
fig1 <- DimPlot(lung, reduction = "umap", label = TRUE, label.size = 4, repel = TRUE, pt.size = 0.01) + NoLegend() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
png(paste0(fig.path, "figure_1.png"), width=1000, height=1000, units="px")
print(fig1)
dev.off()

### figure 2: healthy/fibrotic batch ###
fig2 <- DimPlot(lung, reduction = "umap", split.by = "condition", group.by = "condition", ncol = 1,
                cols = c(viridis(11)[2],viridis(11)[10]),
                pt.size = 0.01) + NoLegend() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
png(paste0(fig.path, "figure_2.png"), width=500, height=1000, units="px")
print(fig2)
dev.off()

### figure 3: study batch ###
fig3 <- DimPlot(lung, reduction = "umap", group.by = "study", split.by = "study", ncol = 1,
                pt.size = 0.01) + scale_color_viridis(option="plasma",discrete=TRUE, name = "Study") + NoLegend() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
png(paste0(fig.path, "figure_3.png"), width=500, height=1000, units="px")
print(fig3)
dev.off()

### figure 5: overall annotation figure ###
panel1 <- fig1
panel2 <- fig2
panel3 <- fig3
fig5 <- panel1 + panel2 + panel3 + plot_layout(widths = c(6,3,2))
png(paste0(fig.path, "figure_5.png"), width=1300,height=800,units="px")
print(fig5)
dev.off()

### figure 6: healthy / fibrotic side by side comparison ###
fig6 <- DimPlot(lung, reduction = "umap", split.by = "condition", group.by = "condition", 
                pt.size = 0.01, cols = c(viridis(11)[2],viridis(11)[10])) + NoLegend() + theme(text = element_text(size = 20),
                                                                                               plot.title = element_text(size = 20),
                                                                                               axis.text.x = element_text(size = 16),
                                                                                               axis.text.y = element_text(size = 16), 
                                                                                               axis.title.x = element_text(size = 20),
                                                                                               axis.title.y = element_text(size = 20),
                                                                                               legend.title = element_text(size = 20),
                                                                                               legend.text = element_text(size = 20))
png(paste0(fig.path, "figure_6.png"), width=1600,height=800,units="px")
print(fig6)
dev.off()

### figure 7: healthy / fibrotic annotation side by side ###
fig7 <- DimPlot(lung, reduction = "umap", split.by = "condition", group.by = "celltype", 
                repel = TRUE, label = TRUE, label.size = 3.5, pt.size = 0.01) + NoLegend() + 
  theme(text = element_text(size = 20),
        plot.title = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
png(paste0(fig.path, "figure_7.png"), width=1600,height=800,units="px")
print(fig7)
dev.off()

### figure 8: overall comparison of healthy vs fibrotic annotation ###
fig8 <- fig6 / fig7
png(paste0(fig.path, "figure_8.png"), width=1200,height=1200,units="px")
print(fig8)
dev.off()