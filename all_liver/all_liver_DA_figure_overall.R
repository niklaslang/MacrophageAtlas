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
DA.path <- "/home/s1987963/ds_group/Niklas/all_liver/DA/overall_figures/"

### read files ###
liver <- readRDS(liver.path)

######################################
### DATA CLEANING AND MANIPULATING ###
######################################
### subset data ###
# only use ramachandran data for downstream analyses
liver <- subset(liver, subset = study == "ramachandran_liver") # 82132 cells
liver <- subset(liver, subset = patient.ID != "Healthy_6" & patient.ID != "Healthy_7" & lineage_integrated != "Red blood cell") # 73091 cells

# reorder lineage factors
liver$lineage_ordered = factor(liver$lineage_integrated)
liver$lineage_ordered <- factor(liver$lineage_ordered, levels(liver$lineage_ordered)[c(4,3,6,1,9,11,8,2,5,7,10)])

### convert seurat objects to sce ###
liver.sce <- as.SingleCellExperiment(liver)

### create meta data data frame ###
colData(liver.sce) %>% 
  as.data.frame %>% 
  transmute(
    Condition = condition, 
    Patient.ID = patient.ID,
    Sample.ID = patient.ID,
    Cluster.ID = seurat_clusters,
    Celltype = celltype_integrated,
    Lineage = lineage_integrated) %>%
  mutate_all(as.factor) %>% 
  set_rownames(colnames(liver.sce)) %>% 
  DataFrame -> colData(liver.sce)

# view data frame
head(colData(liver.sce))

# store number of cluster and number of samples
n.clusters <- length(kids <- set_names(levels(liver.sce$Celltype)))
n.lineages <- length(lids <- set_names(levels(liver.sce$Lineage)))
n.samples <- length(sids <- set_names(levels(liver.sce$Patient.ID)))

# summary of experimental design
m <- match(sids, liver.sce$Sample.ID)
n_cells <- as.numeric(table(liver.sce$Sample.ID))
(ei <- data.frame(colData(liver.sce)[m, ], 
                  n_cells, row.names = NULL) %>% 
    select(-c("Celltype", "Lineage", "Cluster.ID")))

### DA visualisation ###
# calculate lineage abundance per sample
n_lineages <- table(liver.sce$Lineage, liver.sce$Patient.ID)

# calculate lineage proportions across samples
freqs_lineages <- prop.table(n_lineages, margin = 1)

# prepare data.frames for visualisation
df.lineages <- data.frame(
  Frequency = as.numeric(freqs_lineages), 
  Lineage = rep(lids, n.samples),
  Sample.ID = rep(sids, each = n.lineages))
m <- match(df.lineages$Sample.ID, ei$Sample.ID)
df.lineages$Condition <- ei$Condition[m]

# reorder lineage factors
df.lineages$Lineage <- factor(df.lineages$Lineage, levels(df.lineages$Lineage)[c(4,3,6,1,9,11,8,2,5,7,10)])
# reorder condition factors
df.lineages$Condition <- factor(df.lineages$Condition, levels(df.lineages$Condition)[c(2,1)])

# create abundance table
abundance.table <- table(liver.sce$Lineage, liver.sce$Sample.ID) 
abundances.table <- unclass(abundance.table)

# create celltype composition data table: healthy vs fibrotic
df.fibrotic.lineages <- df.lineages[df.lineages$Condition == "fibrotic",]
df.fibrotic.lineages$Pair.ID <- rep(1:5, each = 11)
df.fibrotic.lineages$Pair.ID <- paste0(df.fibrotic.lineages$Pair.ID, "_", rep(1:11, 5))
df.fibrotic.lineages$fibrotic_frequency <- df.fibrotic.lineages$Frequency
fibrotic.lineages.dt <- data.table(df.fibrotic.lineages[c("Lineage", "Pair.ID", "fibrotic_frequency")])
df.healthy.lineages <- df.lineages[df.lineages$Condition == "healthy",]
df.healthy.lineages$Pair.ID <- rep(1:5, each = 11)
df.healthy.lineages$Pair.ID <- paste0(df.healthy.lineages$Pair.ID, "_", rep(1:11, 5))
df.healthy.lineages$healthy_frequency <- df.healthy.lineages$Frequency
healthy.lineages.dt <- data.table(df.healthy.lineages[c("Lineage", "Pair.ID", "healthy_frequency")])

# create data table with composition pairs
abundance.pairs.dt <- merge(fibrotic.lineages.dt, healthy.lineages.dt, by = c("Pair.ID", "Lineage"))

# quantify changes in cell type composition in each cell lineage
abundance.pairs.dt[, ratio := fibrotic_frequency/healthy_frequency, by = Pair.ID]

# log transform changes in cell type composition in each cell lineage
abundance.pairs.dt[, log.ratio := log2(ratio), by = Pair.ID]

# compute CI of log2 changes
abundance.pairs.dt[, sd := sd(log.ratio), by= Lineage]

# convert to data.frame for plotting 
abundance.pairs <- data.frame(abundance.pairs.dt)
# reorder lineage factors 
abundance.pairs$Lineage <- factor(abundance.pairs$Lineage, levels(abundance.pairs$Lineage)[c(11,10,9,8,7,6,5,4,3,2,1)])

#################
### DA FIGURE ###
#################

# seurat colours
seurat.colours <- hue_pal()(11)

# plot UMAP - coloured by lineages
overall.lineage.umap <- DimPlot(liver, reduction = "umap", group.by = "lineage_ordered", label = TRUE, repel = TRUE, pt.size = 0.1) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 16, vjust = -1),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 12))
png(paste0(DA.path, "lineage_umap.png"), width=1000,height=1000,units="px")
print(overall.lineage.umap)
dev.off()

# plot cell lineage composition in health vs disease
overall.lineage.barplot <- ggplot(df.lineages, aes(x = Condition, y = Frequency, fill = Lineage)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic() + 
  labs(x = "Condition", y = "Proportion of each cell lineage", fill = "Cell lineage") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 16, vjust = -1),
        axis.title.y = element_text(size = 16),
        legend.position="none")
png(paste0(DA.path, "lineage_barplot.png"), width=500,height=1000,units="px")
print(overall.lineage.barplot)
dev.off()

# plot changes in cell type composition in health vs disease
overall.lineage.boxplot <- ggplot(abundance.pairs, aes(x = log.ratio, y = Lineage, fill = Lineage)) +
  geom_boxplot(outlier.colour = NA) +
  scale_fill_manual(values = rev(seurat.colours)) +
  stat_boxplot(geom ='errorbar') +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_x_continuous(limits = c(-3, 3)) +
  theme_classic() +
  labs(x = "log2(relative proportion of cell lineage\n fibrotic/healthy)", y = "Lineage", fill = "Cell lineage") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 16, vjust = -1),
        axis.title.y = element_text(size = 16),
        legend.position="none")
png(paste0(DA.path, "lineage_boxplot.png"), width=500,height=1000,units="px")
print(overall.lineage.boxplot)
dev.off()

overall.figure <- overall.lineage.umap + overall.lineage.barplot + overall.lineage.boxplot + plot_layout(widths = c(4,1,2), guides = 'collect')
png(paste0(DA.path, "figure.png"), width=1600,height=800,units="px")
print(overall.figure)
dev.off()