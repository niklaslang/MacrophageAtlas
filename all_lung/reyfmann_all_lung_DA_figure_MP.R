library(Seurat)
library(scater)
library(SingleCellExperiment)
library(edgeR)
library(reticulate)
library(umap)
library(dplyr)
library(magrittr)
library(Matrix)
library(data.table)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(scales)
library(RColorBrewer)
library(viridis)
library(ggsci)
library(patchwork)

### path variables ###
lung.path <- "/home/s1987963/ds_group/Niklas/combined_organs/combined_lung_annotated.rds"
DA.path <- "/home/s1987963/ds_group/Niklas/all_lung/DA/MP_figures/"

### read files ###
lung <- readRDS(lung.path)

######################################
### DATA CLEANING AND MANIPULATING ###
######################################
### subset data ###
# only use reyfman data for downstream analyses
lung <- subset(lung, subset = study == "reyfman_lung")
lung <- subset(lung, subset = lineage_integrated != "Red blood cell" & lineage_integrated != "Doublet") # 70581 cells
# MP subset
MP <- subset(lung, subset = lineage_integrated == "MP") # 32073 cells
MP <- subset(MP, subset = celltype_integrated != "Combined Lung Macrophage 5") # 31834 cells

# add metadata for visualisation
cell.data <- data.table(barcode = colnames(lung),
                        celltype = lung$celltype_integrated,
                        lineage = lung$lineage_integrated)
cell.data[, celltype_ordered := ifelse(lineage == "MP" & celltype != "Combined Lung Macrophage 5", paste0(celltype), "Other")]
cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$celltype <- NULL
cell.data$lineage <- NULL
cell.data$barcode <- NULL
lung <- AddMetaData(lung, cell.data, col.name = "celltype_ordered")

# reorder lineage factors
lung$celltype_ordered <- factor(lung$celltype_ordered)
lung$celltype_ordered <- factor(lung$celltype_ordered, levels(lung$celltype_ordered)[c(11,7,8,9,10,6,1,2,3,4,5,12,13)])

### convert seurat objects to sce ###
MP.sce <- as.SingleCellExperiment(MP)

### create meta data data frame ###
colData(MP.sce) %>% 
  as.data.frame %>% 
  transmute(
    Condition = condition, 
    Patient.ID = patient.ID,
    Sample.ID = patient.ID,
    Cluster.ID = seurat_clusters,
    Celltype = celltype_integrated,
    Lineage = lineage_integrated) %>%
  mutate_all(as.factor) %>% 
  set_rownames(colnames(MP.sce)) %>% 
  DataFrame -> colData(MP.sce)

# view data frame
head(colData(MP.sce))

# store number of cluster and number of samples
n.clusters <- length(kids <- set_names(levels(MP.sce$Celltype)[c(31,32,33,34,35,37,38,39,40,41,42,43)]))
n.samples <- length(sids <- set_names(levels(MP.sce$Patient.ID)))

# summary of experimental design
m <- match(sids, MP.sce$Sample.ID)
n_cells <- as.numeric(table(MP.sce$Sample.ID))

(ei <- data.frame(colData(MP.sce)[m, ], 
                  n_cells, row.names = NULL) %>% 
    select(-c("Celltype", "Lineage", "Cluster.ID")))

# calculate cell type abundance per sample
n_cells <- table(MP.sce$Celltype, MP.sce$Patient.ID)[c(31,32,33,34,35,37,38,39,40,41,42,43),]

# calculate cell/lineage proportions across samples
freqs_cells <- prop.table(n_cells, margin = 1)

# prepare data.frames for visualisation
df.celltypes <- data.frame(
  Frequency = as.numeric(freqs_cells), 
  Celltype = rep(kids, n.samples),
  Sample.ID = rep(sids, each = n.clusters))
m <- match(df.celltypes$Sample.ID, ei$Sample.ID)
df.celltypes$Condition <- ei$Condition[m]

### re-order celltype factor ###
df.celltypes$Celltype <- factor(df.celltypes$Celltype, levels(df.celltypes$Celltype)[c(11,7,8,9,10,6,1,2,3,4,5,12)])

# create celltype composition data table: healthy vs fibrotic
df.MP.fibrotic <- df.celltypes[df.celltypes$Condition == "fibrotic",]
df.MP.fibrotic$Pair.ID <- rep(1:8, each = 12)
df.MP.fibrotic$Pair.ID <- paste0(df.MP.fibrotic$Pair.ID, "_", rep(1:12, 8))
df.MP.fibrotic$fibrotic_frequency <- df.MP.fibrotic$Frequency
MP.fibrotic.dt <- data.table(df.MP.fibrotic[c("Celltype", "Pair.ID", "fibrotic_frequency")])
df.MP.healthy <- df.celltypes[df.celltypes$Condition == "healthy",]
df.MP.healthy$Pair.ID <- rep(1:8, each = 12)
df.MP.healthy$Pair.ID <- paste0(df.MP.healthy$Pair.ID, "_", rep(1:12, 8))
df.MP.healthy$healthy_frequency <- df.MP.healthy$Frequency
MP.healthy.dt <- data.table(df.MP.healthy[c("Celltype", "Pair.ID", "healthy_frequency")])

# create data table with composition pairs
abundance.pairs.dt <- merge(MP.fibrotic.dt, MP.healthy.dt, by = c("Pair.ID", "Celltype"))

# quantify changes in cell type composition in each cell lineage
abundance.pairs.dt[, ratio := fibrotic_frequency/healthy_frequency, by = Pair.ID]

# log transform changes in cell type composition in each cell lineage
abundance.pairs.dt[, log.ratio := log2(ratio), by = Pair.ID]

# compute CI of log2 changes
abundance.pairs.dt[, sd := sd(log.ratio), by= Celltype]

# convert to data.frame for plotting 
abundance.pairs <- data.frame(abundance.pairs.dt)

# reorder celltype factors 
abundance.pairs$Celltype <- factor(abundance.pairs$Celltype, levels(abundance.pairs$Celltype)[c(12,11,10,9,8,7,6,5,4,3,2,1)])

#################
### DA FIGURE ###
#################

# seurat colours
seurat.colours <- hue_pal()(12)

# plot UMAP
MP.celltype.umap <- DimPlot(lung, reduction = "umap", group.by = "celltype_ordered", cols = c(seurat.colours,"#999999"), label = TRUE, repel = TRUE, pt.size = 0.1) +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20, vjust = -1),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 16))
png(paste0(DA.path, "reyfman_celltype_umap.png"), width=1000,height=1000,units="px")
print(MP.celltype.umap)
dev.off()

MP.celltype.barplot <- ggplot(df.celltypes, aes(x = Condition, y = Frequency, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic() + 
  labs(x = "Condition", y = "Proportion of each MP population", fill = "Cell type") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20, vjust = -1),
        axis.title.y = element_text(size = 20),
        legend.position="none")
png(paste0(DA.path, "reyfman_celltype_barplot.png"), width=500,height=1000,units="px")
print(MP.celltype.barplot)
dev.off()

MP.celltype.boxplot <- ggplot(abundance.pairs, aes(x = log.ratio, y = Celltype, fill = Celltype)) +
  geom_boxplot(outlier.colour = NA) +
  scale_fill_manual(values = seurat.colours[c(7,8,9,10,11,6,2,3,4,5,1,12)]) +
  stat_boxplot(geom ='errorbar') +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_segment(aes(x=-0.5, y=12.5, xend=-3, yend=12.5), lineend = "round", linejoin = "mitre", arrow=arrow(), size=0.1) + 
  geom_segment(aes(x=0.5, y=12.5, xend=3, yend=12.5), lineend = "round", linejoin = "mitre", arrow=arrow(), size=0.1) +
  scale_x_continuous(limits = c(-6, 6)) +
  theme_classic() +
  labs(x = "log2(relative proportion of cell type\n fibrotic/healthy)", y = "Cell type", title = "Decreasing | Increasing \n cell type proportion in fibrosis") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20, vjust = -1),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") +
  geom_text(x = 4.5, y = 10.90, label = "***", size = 10) +
  geom_text(x = 4, y = 8.90, label = "*", size = 10) +
  geom_text(x = 5.5, y = 3.90, label = "*", size = 10) +
  geom_text(x = 3, y = 1.90, label = "*", size = 10) +
  png(paste0(DA.path, "reyfman_celltype_boxplot.png"), width=500,height=1000,units="px")
print(MP.celltype.boxplot)
dev.off()

MP.figure <- MP.celltype.umap + MP.celltype.barplot + MP.celltype.boxplot + plot_layout(widths = c(4,1,2), guides = 'collect')
png(paste0(DA.path, "reyfman_figure.png"), width=1700,height=800,units="px")
print(MP.figure)
dev.off()
