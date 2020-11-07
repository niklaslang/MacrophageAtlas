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
fig.path <- "/home/s1987963/ds_group/Niklas/all_blood/harmonize_samples/figures/"

### read blood data ###
blood <- readRDS("/home/s1987963/ds_group/Niklas/combined_organs/combined_blood_annotated.rds")

### remove doublets ###
blood <- subset(blood, subset = lineage_integrated != "Doublet")

### remove erythrocytes ###
blood <- subset(blood, subset = celltype_integrated != "Combined Blood Red Blood cell") # 72619 cells

### figure 1: annotation ###
fig1 <- DimPlot(blood, reduction = "umap", label = TRUE, label.size = 4, repel = TRUE, pt.size = 0.01) + NoLegend() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
png(paste0(fig.path, "figure_1.png"), width=1000, height=1000, units="px")
print(fig1)
dev.off()

### figure 2: healthy/fibrotic batch ###
fig2 <- DimPlot(blood, reduction = "umap", group.by = "condition", ncol = 1,
                pt.size = 0.001) + scale_color_viridis(option="plasma",discrete=TRUE, name = "Condition") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
png(paste0(fig.path, "figure_2.png"), width=1000, height=1000, units="px")
print(fig2)
dev.off()

### figure 3: study batch ###
fig3 <- DimPlot(blood, reduction = "umap", group.by = "study", ncol = 1,cols = c(plasma(11)[9], plasma(11)[3]),
                pt.size = 0.01) + labs(color = "Study") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
png(paste0(fig.path, "figure_3.png"), width=1000, height=1000, units="px")
print(fig3)
dev.off()

### figure 5: overall annotation figure ###
panel1 <- fig1
panel2 <- fig2 / fig3
fig5 <- panel1 + panel2 + plot_layout(widths = c(2,1))
png(paste0(fig.path, "figure_5.png"), width=1300,height=800,units="px")
print(fig5)
dev.off()

### figure 6: healthy / fibrotic side by side comparison ###
fig6 <- DimPlot(blood, reduction = "umap", split.by = "condition", group.by = "condition", 
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
fig7 <- DimPlot(blood, reduction = "umap", split.by = "condition", group.by = "celltype", 
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

### HEATMAP ###
MP <- subset(blood, subset = lineage_integrated == "MP")
MP <- ScaleData(MP)
MP$celltype_integrated <- factor(MP$celltype_integrated, levels(MP$celltype_integrated)[c(11,12,13,14,15,16,17,18)])
MP.features <- c("CD14", "LYZ", "VCAN", "MNDA", "FCGR3A", "CSF1R", "CD68", "CLEC9A", "IDO1", "CD1C", "CLEC10A", "IRF7", "IRF8")
MP.fig1 <- DoHeatmap(subset(MP, downsample = 100), group.by = "celltype_integrated", features = MP.features, size = 5) + 
  theme(text = element_text(size = 20),
        plot.title = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
png(paste0(fig.path, "MP_figure_1.png"), width=1000,height=1000,units="px")
print(MP.fig1)
dev.off()

### FEATURE PLOT ###
MP.fig2 <- FeaturePlot(blood, features = MP.features, pt.size = 0.5, ncol = 4)  + 
  theme(text = element_text(size = 20),
        plot.title = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(fig.path, "MP_figure_2.png"), width=1200,height=900,units="px")
print(MP.fig2)
dev.off()

### MP specific markers ###
MP.specific.features <- c("TREM2", "CD9",
                          "LGALS3", "ACP5", "ANXA2", "ATP6V1F", "C15orf48", "CSTB", "FABP5", "GPNMB", "GSN","HCST", "ITGB2", "SDS")

### dot plot ###
MP.fig5 <- DotPlot(MP, features = rev(MP.specific.features), cols = "RdYlBu", group.by = "celltype_integrated") + RotatedAxis() +
  theme(text = element_text(size = 20),
        plot.title = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
png(paste0(fig.path, "MP_figure_5.png"), width=1200,height=600,units="px")
print(MP.fig5)
dev.off()
