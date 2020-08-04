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
fig.path <- "/home/s1987963/ds_group/Niklas/all_kidney/harmonize_samples/figures/"

### read kidney data ###
kidney <- readRDS("/home/s1987963/ds_group/Niklas/combined_organs/combined_kidney_annotated.rds")

### remove doublets ###
kidney <- subset(kidney, subset = lineage_integrated != "Doublet")

### add technology meta data: scRNA-seq vs snRNA-seq ###
meta.data <- data.table(barcode = colnames(kidney),
                        study = kidney$study)
meta.data[, tech := ifelse(study == "wilson_kidney", "snRNA-seq", "scRNA-seq")]
meta.data <- data.frame(meta.data, row.names = meta.data$barcode)
meta.data$study <- NULL
meta.data$barcode <- NULL
kidney <- AddMetaData(kidney, meta.data, col.name = "tech")

### reorder condition factor ###
kidney$condition <- factor(kidney$condition)
kidney$condition <- factor(kidney$condition, levels(kidney$condition)[c(2,1)])

### figure 1: annotation ###
fig1 <- DimPlot(kidney, reduction = "umap", label = TRUE, label.size = 5, repel = TRUE, pt.size = 0.01) + NoLegend() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
png(paste0(fig.path, "figure_1.png"), width=1000, height=1000, units="px")
print(fig1)
dev.off()

### figure 2: healthy/fibrotic batch ###
fig2 <- DimPlot(kidney, reduction = "umap", group.by = "condition", 
                pt.size = 0.01) + scale_color_viridis(discrete=TRUE, name = "Condition") +
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
fig3 <- DimPlot(kidney, reduction = "umap", group.by = "study", 
                pt.size = 0.01) + scale_color_viridis(option="plasma",discrete=TRUE, name = "Study")+
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
png(paste0(fig.path, "figure_3.png"), width=1000, height=1000, units="px")
print(fig3)
dev.off()

### figure 4: scRNA-seq vs snRNA-seq ###
fig4 <- DimPlot(kidney, reduction = "umap", group.by = "tech", cols = c(plasma(11)[3],plasma(11)[9]),
                pt.size = 0.01) + labs(color = "Technology") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
png(paste0(fig.path, "figure_4.png"), width=1000, height=1000, units="px")
print(fig4)
dev.off()

### figure 5: overall annotation figure ###
panel1 <- fig1
panel2 <- fig2 / fig3 / fig4
fig5 <- panel1 + panel2 + plot_layout(widths = c(3,1))
png(paste0(fig.path, "figure_5.png"), width=1300,height=900,units="px")
print(fig5)
dev.off()

### figure 6: healthy / fibrotic side by side comparison ###
fig6 <- DimPlot(kidney, reduction = "umap", split.by = "condition", group.by = "condition", 
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
fig7 <- DimPlot(kidney, reduction = "umap", split.by = "condition", group.by = "celltype", 
                repel = TRUE, label = TRUE, label.size = 4, pt.size = 0.01) + NoLegend() + 
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