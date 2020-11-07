library(Seurat)
library(scater)
library(SingleCellExperiment)
library(edgeR)
library(reticulate)
library(umap)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)
library(pheatmap)
library(scales)
library(RColorBrewer)
library(viridis)
library(ggsci)
library(patchwork)

### path variables ###
lung.path <- "/home/s1987963/ds_group/Niklas/combined_organs/combined_lung_annotated.rds"
DA.path <- "/home/s1987963/ds_group/Niklas/all_lung/DA/all_studies/"

### read files ###
lung <- readRDS(lung.path)

### subset data ###
# only use reyfman data for downstream analyses
lung <- subset(lung, subset = study == "reyfman_lung")
lung <- subset(lung, subset = lineage_integrated != "Red blood cell" & lineage_integrated != "Doublet") # 70581 cells

# endothelial cells
endothelia <- subset(lung, subset = lineage_integrated == "Endothelia") # 1238 cells

# MPs
MP <- subset(lung, subset = lineage_integrated == "MP") # 32073 cells

# NK/T cells
NKTcell <- subset(lung, subset = lineage_integrated == "T cell" | lineage_integrated == "NK cell") # 1195 cells
Tcell <- subset(lung, subset = lineage_integrated == "T cell") # 1027 cells
NKcell <- subset(lung, subset = lineage_integrated == "NK cell") # 168 cells

# mesenchymal cells
mesenchyme <- subset(lung, subset = lineage_integrated == "Mesenchyme") # 558 cells

### convert seurat objects to sce ###
lung.sce <- as.SingleCellExperiment(lung)
endothelia.sce <- as.SingleCellExperiment(endothelia)
MP.sce <- as.SingleCellExperiment(MP)
NKTcell.sce <- as.SingleCellExperiment(NKTcell)
NKcell.sce <- as.SingleCellExperiment(NKcell)
Tcell.sce <- as.SingleCellExperiment(Tcell)
mesenchyme.sce <- as.SingleCellExperiment(mesenchyme)

#################
### FUNCTIONS ###
#################

composition_barplot <- function(data, level, split, x_label, legend_title){
  plot <- ggplot(data, aes(x = level, y = Frequency, fill = split)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_classic() + 
    labs(x = x_label, fill = legend_title) +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12), 
          axis.title = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12)) +
    scale_fill_viridis(discrete=TRUE)
  return(plot)
}

abundance_barplot <- function(data, level, split, x_label, legend_title){
  plot <- ggplot(data, aes(x = level, y = Frequency, fill = split)) +
    geom_bar(stat = "identity", position = "fill") +
    facet_wrap(~ Condition, scales = "free_x") +
    labs(x = x_label, fill = legend_title) +
    theme_classic() + 
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12), 
          axis.title = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          strip.text.x = element_text(size=16)) +
    scale_fill_viridis(discrete=TRUE)
  return(plot)
}

###############
### OVERALL ###
###############

### create meta data data frame ###
colData(lung.sce) %>% 
  as.data.frame %>% 
  transmute(
    Condition = condition, 
    Patient.ID = patient.ID,
    Sample.ID = patient.ID,
    Cluster.ID = seurat_clusters,
    Celltype = celltype_integrated,
    Lineage = lineage_integrated) %>%
  mutate_all(as.factor) %>% 
  set_rownames(colnames(lung.sce)) %>% 
  DataFrame -> colData(lung.sce)

# view data frame
head(colData(lung.sce))

# re-order cell lineage factor
lung.sce$Lineage <- factor(lung.sce$Lineage, levels(lung.sce$Lineage)[c(1,2,3,4,5,6,7,8,9,10)])

# store number of cluster and number of samples
n.clusters <- length(kids <- set_names(levels(lung.sce$Celltype)))
n.lineages <- length(lids <- set_names(levels(lung.sce$Lineage)))
n.samples <- length(sids <- set_names(levels(lung.sce$Patient.ID)))

# summary of experimental design
m <- match(sids, lung.sce$Sample.ID)
n_cells <- as.numeric(table(lung.sce$Sample.ID))
(ei <- data.frame(colData(lung.sce)[m, ], 
                  n_cells, row.names = NULL) %>% 
    select(-c("Celltype", "Lineage", "Cluster.ID")))

### DA visualisation ###
# calculate cell type abundance per sample
n_cells <- table(lung.sce$Celltype, lung.sce$Patient.ID)
# calculate lineage abundance per sample
n_lineages <- table(lung.sce$Lineage, lung.sce$Patient.ID)

# calculate cell/lineage proportions across samples
freqs_cells <- prop.table(n_cells, margin = 1)
freqs_lineages <- prop.table(n_lineages, margin = 1)

# prepare data.frames for visualisation
df.celltypes <- data.frame(
  Frequency = as.numeric(freqs_cells), 
  Celltype = rep(kids, n.samples),
  Sample.ID = rep(sids, each = n.clusters))
m <- match(df.celltypes$Sample.ID, ei$Sample.ID)
df.celltypes$Condition <- ei$Condition[m]

df.lineages <- data.frame(
  Frequency = as.numeric(freqs_lineages), 
  Lineage = rep(lids, n.samples),
  Sample.ID = rep(sids, each = n.lineages))
m <- match(df.lineages$Sample.ID, ei$Sample.ID)
df.lineages$Condition <- ei$Condition[m]

### visualisation ###
DA.path <- "/home/s1987963/ds_group/Niklas/all_lung/DA/overall_figures/"
# visualise lineage composition at patient level
overall_lineage_patient.composition <- composition_barplot(df.lineages, df.lineages$Lineage, df.lineages$Sample.ID,
                                                           "Cell lineage", "Patient ID")
png(paste0(DA.path, "reyfman_overall_lineage_patient.composition.barplot.png"), width=1200,height=600,units="px") # width = ncols X 100 + 100 for legend
print(overall_lineage_patient.composition)
dev.off()

# visualise lineage composition at condition level
overall_lineage_condition.composition <- composition_barplot(df.lineages, df.lineages$Lineage, df.lineages$Condition,
                                                             "Cell lineage", "Condition")
png(paste0(DA.path, "reyfman_overall_lineage_condition.composition.barplot.png"), width=1200,height=600,units="px")
print(overall_lineage_condition.composition)
dev.off()

# compare composition plots
overall_lineage_composition <- overall_lineage_condition.composition + overall_lineage_patient.composition
png(paste0(DA.path, "reyfman_overall_lineage_composition.barplot.png"), width=1200,height=600,units="px")
print(overall_lineage_composition)
dev.off()

# barplot of relative celltype abundances per sample
overall_celltype_abundance <- abundance_barplot(df.celltypes, df.celltypes$Sample.ID, df.celltypes$Celltype, "Patient ID", "Cell type")
png(paste0(DA.path, "reyfman_overall_celltype_abundance.barplot.png"), width=1200,height=600,units="px")
print(overall_celltype_abundance)
dev.off()

# barplot of relative lineage abundances per sample
overall_lineage_abundance <- abundance_barplot(df.lineages, df.lineages$Sample.ID, df.lineages$Lineage, "Patient ID", "Cell lineage")
png(paste0(DA.path, "reyfman_overall_lineage_abundance.barplot.png"), width=1100,height=600,units="px")
print(overall_lineage_abundance)
dev.off()

### DA analysis ###
# abundance table
abundances <- table(lung.sce$Lineage, lung.sce$Sample.ID) 
abundances <- unclass(abundances) 
head(abundances)

# attaching some column metadata
extra.info <- colData(lung.sce)[match(colnames(abundances), lung.sce$Sample.ID),]
y.ab <- DGEList(abundances, samples=extra.info)
y.ab

# filter out low-abundance labels
keep <- filterByExpr(y.ab, group=y.ab$samples$Condition)
y.ab <- y.ab[keep,]
summary(keep)

# create design matrix
design <- model.matrix(~ factor(Condition), y.ab$samples)

# estimate the NB dipersion for each cell type
y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)
plotBCV(y.ab, cex=1)

# estimate QL dispersion for each cell type
fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)
summary(fit.ab$df.prior)
plotQLDisp(fit.ab, cex=1)

# test for differences in abundance between healthy and fibrotic tissue
results <- glmQLFTest(fit.ab, coef=ncol(design))
results <- topTags(results)
write.csv(results$table, file = paste0(DA.path, "reyfman_overall_lineage_abundance.csv"))

###########
### MPs ###
###########

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
n.clusters <- length(kids <- set_names(levels(MP.sce$Celltype)[c(31,32,33,34,35,36,37,38,39,40,41,42,43)]))
n.samples <- length(sids <- set_names(levels(MP.sce$Patient.ID)))

# summary of experimental design
m <- match(sids, MP.sce$Sample.ID)
n_cells <- as.numeric(table(MP.sce$Sample.ID))

(ei <- data.frame(colData(MP.sce)[m, ], 
                  n_cells, row.names = NULL) %>% 
    select(-c("Celltype", "Lineage", "Cluster.ID")))

# calculate cell type abundance per sample
n_cells <- table(MP.sce$Celltype, MP.sce$Patient.ID)[c(31,32,33,34,35,36,37,38,39,40,41,42,43),]

# calculate cell/lineage proportions across samples
freqs_cells <- prop.table(n_cells, margin = 1)

# prepare data.frames for visualisation
df.celltypes <- data.frame(
  Frequency = as.numeric(freqs_cells), 
  Celltype = rep(kids, n.samples),
  Sample.ID = rep(sids, each = n.clusters))
m <- match(df.celltypes$Sample.ID, ei$Sample.ID)
df.celltypes$Condition <- ei$Condition[m]

### visualisation ###
DA.path <- "/home/s1987963/ds_group/Niklas/all_lung/DA/MP_figures/"
# visualise celltype composition at patient level
MP_celltype_patient.composition <- composition_barplot(df.celltypes, df.celltypes$Celltype, df.celltypes$Sample.ID, "Cell type", "Patient ID")
png(paste0(DA.path, "reyfman_MP_celltype_patient.composition.barplot.png"), width=1000,height=600,units="px")
print(MP_celltype_patient.composition)
dev.off()

# visualise celltype composition at condition level
MP_celltype_condition.composition <- composition_barplot(df.celltypes, df.celltypes$Celltype, df.celltypes$Condition, "Cell type", "Condition")
png(paste0(DA.path, "reyfman_MP_celltype_condition.composition.barplot.png"), width=1000,height=600,units="px")
print(MP_celltype_condition.composition)
dev.off()

# compare composition plots
MP_celltype_composition <- MP_celltype_condition.composition + MP_celltype_patient.composition
png(paste0(DA.path, "reyfman_MP_celltype_composition.barplot.png"), width=1000,height=600,units="px")
print(MP_celltype_composition)
dev.off()

# visualise cell type abundance
# barplot of relative  celltype abundances
MP_celltype_abundance <- abundance_barplot(df.celltypes, df.celltypes$Sample.ID, df.celltypes$Celltype, "Patient ID", "Cell type")
png(paste0(DA.path, "reyfman_MP_celltype_abundance.barplot.png"), width=1100,height=600,units="px")
print(MP_celltype_abundance)
dev.off()

### DA analysis ###
# abundance table
abundances <- table(MP.sce$Celltype, MP.sce$Sample.ID) 
abundances <- unclass(abundances) 
head(abundances)

# attaching some column metadata.
extra.info <- colData(lung.sce)[match(colnames(abundances), lung.sce$Sample.ID),]
y.ab <- DGEList(abundances, samples=extra.info)
y.ab

# filter out low-abundance labels
keep <- filterByExpr(y.ab, group=y.ab$samples$Condition)
y.ab <- y.ab[keep,]
summary(keep)

# create design matrix
design <- model.matrix(~ factor(Condition), y.ab$samples)

# estimate the NB dipersion for each cell type
y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)
plotBCV(y.ab, cex=1)

# estimate QL dispersion for each cell type
fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)
summary(fit.ab$df.prior)
plotQLDisp(fit.ab, cex=1)

# test for differences in abundance between healthy and fibrotic tissue
res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))
results <- topTags(res)
write.csv(results$table, file = paste0(DA.path, "reyfman_MP_celltype_abundance.csv"))

##################
### ENDOTHELIA ###
##################

### create meta data data frame ###
colData(endothelia.sce) %>% 
  as.data.frame %>% 
  transmute(
    Condition = condition, 
    Patient.ID = patient.ID,
    Sample.ID = patient.ID,
    Cluster.ID = seurat_clusters,
    Celltype = celltype_integrated,
    Lineage = lineage_integrated) %>%
  mutate_all(as.factor) %>% 
  set_rownames(colnames(endothelia.sce)) %>% 
  DataFrame -> colData(endothelia.sce)

# view data frame
head(colData(endothelia.sce))

# store number of cluster and number of samples
n.clusters <- length(kids <- set_names(levels(endothelia.sce$Celltype)[c(13,14,15,16,17)]))
n.samples <- length(sids <- set_names(levels(endothelia.sce$Patient.ID)))

# summary of experimental design
m <- match(sids, endothelia.sce$Sample.ID)
n_cells <- as.numeric(table(endothelia.sce$Sample.ID))

(ei <- data.frame(colData(endothelia.sce)[m, ], 
                  n_cells, row.names = NULL) %>% 
    select(-c("Celltype", "Lineage", "Cluster.ID")))

### DA visualisation ###
# calculate cell type abundance per sample
n_cells <- table(endothelia.sce$Celltype, endothelia.sce$Patient.ID)[c(13,14,15,16,17),]

# calculate cell/lineage proportions across samples
freqs_cells <- prop.table(n_cells, margin = 1)

# prepare data.frames for visualisation
df.celltypes <- data.frame(
  Frequency = as.numeric(freqs_cells), 
  Celltype = rep(kids, n.samples),
  Sample.ID = rep(sids, each = n.clusters))
m <- match(df.celltypes$Sample.ID, ei$Sample.ID)
df.celltypes$Condition <- ei$Condition[m]

### visualisation ###
DA.path <- "/home/s1987963/ds_group/Niklas/all_lung/DA/endothelia_figures/"
# visualise celltype composition at patient level
endothelia_celltype_patient.composition <- composition_barplot(df.celltypes, df.celltypes$Celltype, df.celltypes$Sample.ID, "Cell type", "Patient ID")
png(paste0(DA.path, "reyfman_endothelia_celltype_patient.composition.barplot.png"), width=1000,height=600,units="px")
print(endothelia_celltype_patient.composition)
dev.off()

# visualise celltype composition at condition level
endothelia_celltype_condition.composition <- composition_barplot(df.celltypes, df.celltypes$Celltype, df.celltypes$Condition, "Cell type", "Condition")
png(paste0(DA.path, "reyfman_endothelia_celltype_condition.composition.barplot.png"), width=1000,height=600,units="px")
print(endothelia_celltype_condition.composition)
dev.off()

# compare composition plots
endothelia_celltype_composition <- endothelia_celltype_condition.composition + endothelia_celltype_patient.composition
png(paste0(DA.path, "reyfman_endothelia_celltype_composition.barplot.png"), width=1000,height=600,units="px")
print(endothelia_celltype_composition)
dev.off()

# barplot of relative  celltype abundances
endothelia_celltype_abundance <- abundance_barplot(df.celltypes, df.celltypes$Sample.ID, df.celltypes$Celltype, "Patient ID", "Cell type")
png(paste0(DA.path, "reyfman_endothelia_celltype_abundance.barplot.png"), width=1100,height=600,units="px")
print(endothelia_celltype_abundance)
dev.off()

### DA analysis ###
# abundance table
abundances <- table(endothelia.sce$Celltype, endothelia.sce$Sample.ID) 
abundances <- unclass(abundances) 
head(abundances)

# attaching some column metadata.
extra.info <- colData(endothelia.sce)[match(colnames(abundances), endothelia.sce$Sample.ID),]
y.ab <- DGEList(abundances, samples=extra.info)
y.ab

# filter out low-abundance labels
keep <- filterByExpr(y.ab, group=y.ab$samples$Condition)
y.ab <- y.ab[keep,]
summary(keep)

# create design matrix
design <- model.matrix(~ factor(Condition), y.ab$samples)

# estimate the NB dipersion for each cell type
y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)
plotBCV(y.ab, cex=1)

# estimate QL dispersion for each cell type
fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)
summary(fit.ab$df.prior)
plotQLDisp(fit.ab, cex=1)

# test for differences in abundance between healthy and fibrotic tissue
res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))
results <- topTags(res)
write.csv(results$table, file = paste0(DA.path, "reyfman_endothelia_celltype_abundance.csv"))

##################
### MESENCHYME ###
##################

### create meta data data frame ###
colData(mesenchyme.sce) %>% 
  as.data.frame %>% 
  transmute(
    Condition = condition, 
    Patient.ID = patient.ID,
    Sample.ID = patient.ID,
    Cluster.ID = seurat_clusters,
    Celltype = celltype_integrated,
    Lineage = lineage_integrated) %>%
  mutate_all(as.factor) %>% 
  set_rownames(colnames(mesenchyme.sce)) %>% 
  DataFrame -> colData(mesenchyme.sce)

# view data frame
head(colData(mesenchyme.sce))

# store number of cluster and number of samples
n.clusters <- length(kids <- set_names(levels(mesenchyme.sce$Celltype)[c(18,19,20)]))
n.samples <- length(sids <- set_names(levels(mesenchyme.sce$Patient.ID)))

# summary of experimental design
m <- match(sids, mesenchyme.sce$Sample.ID)
n_cells <- as.numeric(table(mesenchyme.sce$Sample.ID))

(ei <- data.frame(colData(mesenchyme.sce)[m, ], 
                  n_cells, row.names = NULL) %>% 
    select(-c("Celltype", "Lineage", "Cluster.ID")))

# calculate cell type abundance per sample
n_cells <- table(mesenchyme.sce$Celltype, mesenchyme.sce$Patient.ID)[c(18,19,20),]

# calculate cell/lineage proportions across samples
freqs_cells <- prop.table(n_cells, margin = 1)

# prepare data.frames for visualisation
df.celltypes <- data.frame(
  Frequency = as.numeric(freqs_cells), 
  Celltype = rep(kids, n.samples),
  Sample.ID = rep(sids, each = n.clusters))
m <- match(df.celltypes$Sample.ID, ei$Sample.ID)
df.celltypes$Condition <- ei$Condition[m]

### visualisation ###
DA.path <- "/home/s1987963/ds_group/Niklas/all_lung/DA/mesenchyme_figure/"
# visualise celltype composition at patient level
mesenchyme_celltype_patient.composition <- composition_barplot(df.celltypes, df.celltypes$Celltype, df.celltypes$Sample.ID, "Cell type", "Patient ID")
png(paste0(DA.path, "reyfman_mesenchyme_celltype_patient.composition.barplot.png"), width=300,height=600,units="px")
print(mesenchyme_celltype_patient.composition)
dev.off()

# visualise celltype composition at condition level
mesenchyme_celltype_condition.composition <- composition_barplot(df.celltypes, df.celltypes$Celltype, df.celltypes$Condition, "Cell type", "Condition")
png(paste0(DA.path, "reyfman_mesenchyme_celltype_condition.composition.barplot.png"), width=300,height=600,units="px")
print(mesenchyme_celltype_condition.composition)
dev.off()

# compare composition plots
mesenchyme_celltype_composition <- mesenchyme_celltype_condition.composition + mesenchyme_celltype_patient.composition
png(paste0(DA.path, "reyfman_mesenchyme_celltype_composition.barplot.png"), width=600,height=600,units="px")
print(mesenchyme_celltype_composition)
dev.off()

# visualise cell type abundance
# barplot of relative  celltype abundances
mesenchyme_celltype_abundance <- abundance_barplot(df.celltypes, df.celltypes$Sample.ID, df.celltypes$Celltype, "Patient ID", "Cell type")
png(paste0(DA.path, "reyfman_mesenchyme_celltype_abundance.barplot.png"), width=1100,height=600,units="px")
print(mesenchyme_celltype_abundance)
dev.off()

### DA analysis ###
# doesn't work for 2 cell types only :(