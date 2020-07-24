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
liver.path <- "/home/s1987963/ds_group/Niklas/combined_organs/combined_liver_annotated.rds"
DA.path <- "/home/s1987963/ds_group/Niklas/all_liver/DA/NKT_figures/"

### read files ###
liver <- readRDS(liver.path)

######################################
### DATA CLEANING AND MANIPULATING ###
######################################
### subset data ###
# only use ramachandran data for downstream analyses
liver <- subset(liver, subset = study == "ramachandran_liver") # 82132 cells
liver <- subset(liver, subset = patient.ID != "Healthy_6" & patient.ID != "Healthy_7" & lineage_integrated != "Red blood cell") # 73091 cells
# add metadata for visualisation
cell.data <- data.table(barcode = colnames(liver),
                        celltype = Idents(liver),
                        lineage = liver$lineage_integrated)
cell.data[, celltype_ordered := ifelse(lineage == "T cell" | lineage == "NK cell", paste0(celltype), "Other")]
cell.data[, celltype_ordered := ifelse(celltype_ordered == "Combined Liver NK Cell 4", "Other", paste0(celltype_ordered))]
cell.data[, celltype_ordered := ifelse(celltype_ordered == "Combined Liver NK Cell 3", "Combined Liver NK cell 3", paste0(celltype_ordered))]
cell.data[, celltype_ordered := ifelse(celltype_ordered == "Combined Liver NK1", "Combined Liver NK cell 1", paste0(celltype_ordered))]
cell.data[, celltype_ordered := ifelse(celltype_ordered == "Combined Liver NK2", "Combined Liver NK cell 2", paste0(celltype_ordered))]
cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$celltype <- NULL
cell.data$lineage <- NULL
cell.data$barcode <- NULL
liver <- AddMetaData(liver, cell.data, col.name = "celltype_ordered")
# reorder lineage factors
liver$celltype_ordered = factor(liver$celltype_ordered)

# NKT subset
NKT <- subset(liver, subset = celltype_ordered != "Other")
#NKT <- subset(liver, subset = lineage_integrated == "T cell" | lineage_integrated == "NK cell" & celltype_ordered != )

### convert seurat objects to sce ###
NKT.sce <- as.SingleCellExperiment(NKT)

### create meta data data frame ###
colData(NKT.sce) %>% 
  as.data.frame %>% 
  transmute(
    Condition = condition, 
    Patient.ID = patient.ID,
    Sample.ID = patient.ID,
    Cluster.ID = seurat_clusters,
    Celltype = celltype_ordered,
    Lineage = lineage_integrated) %>%
  mutate_all(as.factor) %>% 
  set_rownames(colnames(NKT.sce)) %>% 
  DataFrame -> colData(NKT.sce)

# view data frame
head(colData(NKT.sce))

# store number of cluster and number of samples
n.clusters <- length(kids <- set_names(levels(NKT.sce$Celltype)[c(1,2,3,4,5,6,7,8)]))
n.samples <- length(sids <- set_names(levels(NKT.sce$Patient.ID)))

# summary of experimental design
m <- match(sids, NKT.sce$Sample.ID)
n_cells <- as.numeric(table(NKT.sce$Sample.ID))

(ei <- data.frame(colData(NKT.sce)[m, ], 
                  n_cells, row.names = NULL) %>% 
    select(-c("Celltype", "Lineage", "Cluster.ID")))

# calculate cell type abundance per sample
n_cells <- table(NKT.sce$Celltype, NKT.sce$Patient.ID)[c(1,2,3,4,5,6,7,8),]

# calculate cell/lineage proportions across samples
freqs_cells <- prop.table(n_cells, margin = 1)

# prepare data.frames for visualisation
df.celltypes <- data.frame(
  Frequency = as.numeric(freqs_cells), 
  Celltype = rep(kids, n.samples),
  Sample.ID = rep(sids, each = n.clusters))
m <- match(df.celltypes$Sample.ID, ei$Sample.ID)
df.celltypes$Condition <- ei$Condition[m]

# reorder lineage factors
df.celltypes$Celltype <- factor(df.celltypes$Celltype, levels(df.celltypes$Celltype))
# reorder condition factors
df.celltypes$Condition <- factor(df.celltypes$Condition, levels(df.celltypes$Condition)[c(2,1)])

# create celltype composition data table: healthy vs fibrotic
df.NKT.fibrotic <- df.celltypes[df.celltypes$Condition == "fibrotic",]
df.NKT.fibrotic$Pair.ID <- rep(1:5, each = 8)
df.NKT.fibrotic$Pair.ID <- paste0(df.NKT.fibrotic$Pair.ID, "_", rep(1:8, 5))
df.NKT.fibrotic$fibrotic_frequency <- df.NKT.fibrotic$Frequency
NKT.fibrotic.dt <- data.table(df.NKT.fibrotic[c("Celltype", "Pair.ID", "fibrotic_frequency")])
df.NKT.healthy <- df.celltypes[df.celltypes$Condition == "healthy",]
df.NKT.healthy$Pair.ID <- rep(1:5, each = 8)
df.NKT.healthy$Pair.ID <- paste0(df.NKT.healthy$Pair.ID, "_", rep(1:8, 5))
df.NKT.healthy$healthy_frequency <- df.NKT.healthy$Frequency
NKT.healthy.dt <- data.table(df.NKT.healthy[c("Celltype", "Pair.ID", "healthy_frequency")])

# create data table with composition pairs
abundance.pairs.dt <- merge(NKT.fibrotic.dt, NKT.healthy.dt, by = c("Pair.ID", "Celltype"))

# quantify changes in cell type composition in each cell lineage
abundance.pairs.dt[, ratio := fibrotic_frequency/healthy_frequency, by = Pair.ID]

# log transform changes in cell type composition in each cell lineage
abundance.pairs.dt[, log.ratio := log2(ratio), by = Pair.ID]

# compute CI of log2 changes
abundance.pairs.dt[, sd := sd(log.ratio), by= Celltype]

# convert to data.frame for plotting 
abundance.pairs <- data.frame(abundance.pairs.dt)

# reorder celltype factors 
abundance.pairs$Celltype <- factor(abundance.pairs$Celltype, levels(abundance.pairs$Celltype)[c(8,7,6,5,4,3,2,1)])

#################
### DA FIGURE ###
#################

# seurat colours
seurat.colours <- hue_pal()(8)

# plot UMAP
NKT.celltype.umap <- DimPlot(liver, reduction = "umap", group.by = "celltype_ordered", cols = c(seurat.colours,"#999999"), label = TRUE, repel = TRUE, pt.size = 0.1) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 16, vjust = -1),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 12))
png(paste0(DA.path, "celltype_umap.png"), width=1000,height=1000,units="px")
print(NKT.celltype.umap)
dev.off()

NKT.celltype.barplot <- ggplot(df.celltypes, aes(x = Condition, y = Frequency, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic() + 
  labs(x = "Condition", y = "Proportion of each NKT cell population", fill = "Cell type") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 16, vjust = -1),
        axis.title.y = element_text(size = 16),
        legend.position="none")
png(paste0(DA.path, "celltype_barplot.png"), width=500,height=1000,units="px")
print(NKT.celltype.barplot)
dev.off()

NKT.celltype.boxplot <- ggplot(abundance.pairs, aes(x = log.ratio, y = Celltype, fill = Celltype)) +
  geom_boxplot(outlier.colour = NA) +
  scale_fill_manual(values = rev(seurat.colours)) +
  stat_boxplot(geom ='errorbar') +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_x_continuous(limits = c(-4, 4)) +
  geom_segment(aes(x=-0.5, y=8.5, xend=-3.5, yend=8.5), lineend = "round", linejoin = "mitre", arrow=arrow(), size=0.2) + 
  geom_segment(aes(x=0.5, y=8.5, xend=3.5, yend=8.5), lineend = "round", linejoin = "mitre", arrow=arrow(), size=0.2) + 
  theme_classic() +
  labs(x = "log2(relative proportion of cell type\n fibrotic/healthy)", y = "Cell type", title = " Increasing | Decreasing\n cell type proportion in fibrosis") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 16, vjust = -1),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none")
png(paste0(DA.path, "celltype_boxplot.png"), width=500,height=1000,units="px")
print(NKT.celltype.boxplot)
dev.off()

NKT.figure <- NKT.celltype.umap + NKT.celltype.barplot + NKT.celltype.boxplot + plot_layout(widths = c(4,1,2), guides = 'collect')
png(paste0(DA.path, "figure.png"), width=1600,height=800,units="px")
print(NKT.figure)
dev.off()