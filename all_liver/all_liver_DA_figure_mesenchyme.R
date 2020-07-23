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
DA.path <- "/home/s1987963/ds_group/Niklas/all_liver/DA/mesenchyme_figures/"

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
cell.data[, celltype_ordered := ifelse(lineage == "Mesenchyme", paste0(celltype), "Other")]
cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$celltype <- NULL
cell.data$lineage <- NULL
cell.data$barcode <- NULL
liver <- AddMetaData(liver, cell.data, col.name = "celltype_ordered")

# mesenchyme subset
mesenchyme <- subset(liver, subset = lineage_integrated== "Mesenchyme") # 2578 cells

# reorder lineage factors
liver$celltype_ordered = factor(liver$celltype_ordered)

### convert seurat objects to sce ###
mesenchyme.sce <- as.SingleCellExperiment(mesenchyme)

### create meta data data frame ###
colData(mesenchyme.sce) %>% 
  as.data.frame %>% 
  transmute(
    Condition = condition, 
    Patient.ID = patient.ID,
    Sample.ID = patient.ID,
    Cluster.ID = seurat_clusters,
    Celltype = celltype_ordered,
    Lineage = lineage_integrated) %>%
  mutate_all(as.factor) %>% 
  set_rownames(colnames(mesenchyme.sce)) %>% 
  DataFrame -> colData(mesenchyme.sce)

# view data frame
head(colData(mesenchyme.sce))

# store number of cluster and number of samples
n.clusters <- length(kids <- set_names(levels(mesenchyme.sce$Celltype)))
n.samples <- length(sids <- set_names(levels(mesenchyme.sce$Patient.ID)))

# summary of experimental design
m <- match(sids, mesenchyme.sce$Sample.ID)
n_cells <- as.numeric(table(mesenchyme.sce$Sample.ID))

(ei <- data.frame(colData(mesenchyme.sce)[m, ], 
                  n_cells, row.names = NULL) %>% 
    select(-c("Celltype", "Lineage", "Cluster.ID")))

# calculate cell type abundance per sample
n_cells <- table(mesenchyme.sce$Celltype, mesenchyme.sce$Patient.ID)

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
df.mesenchyme.fibrotic <- df.celltypes[df.celltypes$Condition == "fibrotic",]
df.mesenchyme.fibrotic$Pair.ID <- rep(1:5, each = 2)
df.mesenchyme.fibrotic$Pair.ID <- paste0(df.mesenchyme.fibrotic$Pair.ID, "_", rep(1:2, 5))
df.mesenchyme.fibrotic$fibrotic_frequency <- df.mesenchyme.fibrotic$Frequency
mesenchyme.fibrotic.dt <- data.table(df.mesenchyme.fibrotic[c("Celltype", "Pair.ID", "fibrotic_frequency")])
df.mesenchyme.healthy <- df.celltypes[df.celltypes$Condition == "healthy",]
df.mesenchyme.healthy$Pair.ID <- rep(1:5, each = 2)
df.mesenchyme.healthy$Pair.ID <- paste0(df.mesenchyme.healthy$Pair.ID, "_", rep(1:2, 5))
df.mesenchyme.healthy$healthy_frequency <- df.mesenchyme.healthy$Frequency
mesenchyme.healthy.dt <- data.table(df.mesenchyme.healthy[c("Celltype", "Pair.ID", "healthy_frequency")])

# create data table with composition pairs
abundance.pairs.dt <- merge(mesenchyme.fibrotic.dt, mesenchyme.healthy.dt, by = c("Pair.ID", "Celltype"))

# quantify changes in cell type composition in each cell lineage
abundance.pairs.dt[, ratio := fibrotic_frequency/healthy_frequency, by = Pair.ID]

# log transform changes in cell type composition in each cell lineage
abundance.pairs.dt[, log.ratio := log2(ratio), by = Pair.ID]

# compute CI of log2 changes
abundance.pairs.dt[, sd := sd(log.ratio), by= Celltype]

# convert to data.frame for plotting 
abundance.pairs <- data.frame(abundance.pairs.dt)

# reorder celltype factors 
abundance.pairs$Celltype <- factor(abundance.pairs$Celltype, levels(abundance.pairs$Celltype)[c(2,1)])

#################
### DA FIGURE ###
#################

# seurat colours
seurat.colours <- hue_pal()(2)

# plot UMAP
mesenchyme.celltype.umap <- DimPlot(liver, reduction = "umap", group.by = "celltype_ordered", cols = c(seurat.colours,"#999999"), label = TRUE, repel = TRUE, pt.size = 0.1) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 16, vjust = -1),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 12))
png(paste0(DA.path, "celltype_umap.png"), width=1000,height=1000,units="px")
print(mesenchyme.celltype.umap)
dev.off()

mesenchyme.celltype.barplot <- ggplot(df.celltypes, aes(x = Condition, y = Frequency, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic() + 
  labs(x = "Condition", y = "Proportion of each mesenchyme population", fill = "Cell type") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 16, vjust = -1),
        axis.title.y = element_text(size = 16),
        legend.position="none")
png(paste0(DA.path, "celltype_barplot.png"), width=500,height=1000,units="px")
print(mesenchyme.celltype.barplot)
dev.off()

mesenchyme.celltype.boxplot <- ggplot(abundance.pairs, aes(x = log.ratio, y = Celltype, fill = Celltype)) +
  geom_boxplot(outlier.colour = NA) +
  scale_fill_manual(values = rev(seurat.colours)) +
  stat_boxplot(geom ='errorbar') +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_x_continuous(limits = c(-3, 3)) +
  theme_classic() +
  labs(x = "log2(relative proportion of cell type\n fibrotic/healthy)", y = "Cell type") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 16, vjust = -1),
        axis.title.y = element_text(size = 16),
        legend.position="none")
png(paste0(DA.path, "celltype_boxplot.png"), width=500,height=1000,units="px")
print(mesenchyme.celltype.boxplot)
dev.off()

patch1 <- mesenchyme.celltype.umap
patch2 <- mesenchyme.celltype.barplot / mesenchyme.celltype.boxplot
mesenchyme.figure <- patch1 + patch2 + plot_layout(widths = c(4,1), guides = 'collect')
png(paste0(DA.path, "figure.png"), width=1400,height=800,units="px")
print(mesenchyme.figure)
dev.off()