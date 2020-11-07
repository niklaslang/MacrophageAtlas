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
fig.path <- "/home/s1987963/ds_group/Niklas/all_liver/harmonize_samples/figures/"

### read liver data ###
liver <- readRDS("/home/s1987963/ds_group/Niklas/combined_organs/combined_liver_annotated.rds")

### reorder condition factor ###
liver$condition <- factor(liver$condition)
liver$condition <- factor(liver$condition, levels(liver$condition)[c(2,1)])

### reorder lineage factor ###
liver$lineage_integrated <- factor(liver$lineage_integrated)
liver$lineage_integrated <- factor(liver$lineage_integrated, levels(liver$lineage_integrated)[c(4,3,6,1,9,12,8,2,5,7,10,11)])

### correct cell-type meta data
cell.data <- data.table(barcode = colnames(liver),
                        celltype = Idents(liver),
                        lineage = liver$lineage_integrated)
cell.data[, celltype_integrated := ifelse(celltype == "Combined Liver NK Cell 4", "Combined Liver NK cell 4", paste0(celltype))]
cell.data[, celltype_integrated := ifelse(celltype_integrated == "Combined Liver NK Cell 3", "Combined Liver NK cell 3", paste0(celltype_integrated))]
cell.data[, celltype_integrated := ifelse(celltype_integrated == "Combined Liver NK1", "Combined Liver NK cell 1", paste0(celltype_integrated))]
cell.data[, celltype_integrated := ifelse(celltype_integrated == "Combined Liver NK2", "Combined Liver NK cell 2", paste0(celltype_integrated))]
cell.data[, celltype_integrated := ifelse(celltype_integrated == "Combined Liver CD8+ T cell 1", "Combined Liver CD8+ T cell", paste0(celltype_integrated))]
cell.data[, celltype_integrated := ifelse(celltype_integrated == "Combined Liver Plasma cell 1", "Combined Liver Plasma cell", paste0(celltype_integrated))]
cell.data[, celltype_integrated := ifelse(celltype_integrated == "Combined Liver Proliferating 3", "Combined Liver Proliferating 2", paste0(celltype_integrated))]
cell.data[, celltype_integrated := ifelse(celltype_integrated == "Combined Liver Proliferating 4", "Combined Liver Proliferating 3", paste0(celltype_integrated))]
cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$celltype <- NULL
cell.data$lineage <- NULL
cell.data$barcode <- NULL
liver <- AddMetaData(liver, cell.data, col.name = "celltype_integrated")

### reorder celltype factor ###
liver$celltype_integrated <- factor(liver$celltype_integrated)
liver$celltype_integrated <- factor(liver$celltype_integrated, levels(liver$celltype_integrated)[c(17,18,19,20, # Hepatocyte 1-4
                                                                                                   8,9,10,11, # Cholangiocytes 1-4
                                                                                                   14,15, # Endo 1-2
                                                                                                   25,26,27, # LSEndo 1-3
                                                                                                   28, # LymphEndo
                                                                                                   16, #HA Endo
                                                                                                   44, # PV Endo
                                                                                                   13, # CV Endo
                                                                                                   31,32, # Mesenchyme 1-2
                                                                                                   1, # B cell
                                                                                                   40, # Plasma cell
                                                                                                   3,4,# CD4+ T cell 1-2
                                                                                                   5, # CD8+ T cell
                                                                                                   12, # CTLA4+ T cell
                                                                                                   21, # IFN primed T cell
                                                                                                   35,36,37,38, # NK cell 1-4
                                                                                                   2, # Basophil
                                                                                                   30, # Mast cell
                                                                                                   33,34, # Monocyte 1-2
                                                                                                   29, # Macrophage 1 
                                                                                                   22,23,24, # KC1-3
                                                                                                   6, # cDC1
                                                                                                   7, # cDC2
                                                                                                   39, # pDC
                                                                                                   41,42,43, # Proliferating 1-4
                                                                                                   45 # RBC
                                                                                                   )])

### remove erythrocytes ###
liver <- subset(liver, subset = celltype_integrated != "Combined Liver Red blood cell")

### figure 1: annotation ###
fig1 <- DimPlot(liver, reduction = "umap", label = TRUE, label.size = 5, repel = TRUE, pt.size = 0.01) + NoLegend() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
png(paste0(fig.path, "figure_1.png"), width=1000, height=1000, units="px")
print(fig1)
dev.off()

### figure 2: healthy/fibrotic batch ###
fig2 <- DimPlot(liver, reduction = "umap", split.by = "condition", group.by = "study", ncol=1, 
                cols = c(viridis(21)[2],viridis(21)[20]),
                pt.size = 0.001) + labs(color = "Study") + #scale_color_viridis(discrete=TRUE, name = "Condition") +
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
fig3 <- DimPlot(liver, reduction = "umap", group.by = "study", cols = c(plasma(11)[3],plasma(11)[9]),
                pt.size = 0.001) + labs(color = "Study") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
png(paste0(fig.path, "figure_3.png"), width=1000, height=1000, units="px")
print(fig3)
dev.off()

### figure 4: overall annotation figure ###
panel1 <- fig1
panel2 <- fig2
fig4 <- panel1 + panel2 + plot_layout(widths = c(2,1))
png(paste0(fig.path, "figure_4.png"), width=1300,height=800,units="px")
print(fig4)
dev.off()

### figure 6: healthy / fibrotic side by side comparison ###
fig6 <- DimPlot(liver, reduction = "umap", split.by = "condition", group.by = "condition", 
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
fig7 <- DimPlot(liver, reduction = "umap", split.by = "condition", group.by = "celltype", 
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

### prepare data for marker genes ###
MP <- subset(liver, subset = lineage_integrated == "MP")
MP$celltype_integrated <- factor(MP$celltype_integrated, levels(MP$celltype_integrated)[c(33,34,35,36,37,38,39,40,41)])

MP.features <- c("CD68", 
                 "S100A8","S100A9", "VCAN",
                 "MARCO", "CD5L", "CD163", "TIMD4",
                 "LYZ", "MNDA",
                 "TREM2", "CD9",
                 "XCR1", "CLEC9A", "IDO1",
                 "CD1C","CLEC10A",
                 "LILRA4","IRF7", "IRF8")

### HEATMAP ###
MP.fig1 <- DoHeatmap(subset(MP, downsample = 1000), group.by = "celltype_integrated", features = MP.features, size = 5) + 
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

### violin plot ###
## helper function 1 ##
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, group.by = "celltype_integrated", ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 90), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

panel1 <- StackedVlnPlot(obj = MP, features = MP.features[1:10])
panel2 <- StackedVlnPlot(obj = MP, features = MP.features[11:20])
MP.fig2 <- panel1 | panel2
png(paste0(fig.path, "MP_figure_2.png"), width=1200,height=1200,units="px")
print(MP.fig2)
dev.off()

### FEATURE UMAP HEATMAP ###
MP.fig3 <- FeaturePlot(subset(liver, downsample = 1000), features = MP.features, ncol = 3) & scale_colour_viridis()
  theme(text = element_text(size = 20),
        plot.title = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) + plot_layout(guides = 'collect')
png(paste0(fig.path, "MP_figure_3.png"), width=1200,height=1200,units="px")
print(MP.fig3)
dev.off()

### MP specific markers ###
MP.specific.features <- c("CD163", "MARCO", "CD5L", "TIMD4",
                          "C1QA", "C1QB", "C1QC", "LIPA", "APOE", "MRC1",
                          "TREM2", "CD9")
MP.fig4 <- FeaturePlot(liver, features = MP.specific.features, pt.size = 0.5, ncol = 5)  + 
  theme(text = element_text(size = 20),
        plot.title = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
png(paste0(fig.path, "MP_figure_4.png"), width=1200,height=600,units="px")
print(MP.fig4)
dev.off()

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

#MP.panel1 <- StackedVlnPlot(obj = MP, features = MP.specific.features[1:5])
#MP.panel2 <- StackedVlnPlot(obj = MP, features = MP.specific.features[6:10])
#MP.fig4 <- MP.panel1 | MP.panel2