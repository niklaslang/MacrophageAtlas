library(Seurat)
library(scater)
library(SingleCellExperiment)
library(data.table)
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
library(patchwork)

### path variables ###
liver.path <- "/home/s1987963/ds_group/Niklas/combined_organs/combined_liver_annotated.rds"
DA.path <- "/home/s1987963/ds_group/Niklas/all_liver/DA/"

### read files ###
liver <- readRDS(liver.path)

### subset data ###
# only use ramachandran data for downstream analyses
liver <- subset(liver, subset = study == "ramachandran_liver") # 82132 cells
liver <- subset(liver, subset = patient.ID != "Healthy_6" & patient.ID != "Healthy_7" & lineage_integrated != "Red blood cell") # 73091 cells

# immune cells
immune <- subset(liver, subset = lineage_integrated %in% c("MP", "Basophil", "T cell", "B cell", "NK cell", "Mast cell", "Plasma cell")) # 52616 cells

# endothelial cells
endothelial <- subset(liver, subset = lineage_integrated == "Endothelia") # 12803 cells

# MPs
MP <- subset(liver, subset = lineage_integrated == "MP") # 11999 cells

# T cells
NKTcell <- subset(liver, subset = lineage_integrated == "T cell" | lineage_integrated == "NK cell") # 37253 cells

### convert seurat objects to sce ###
liver.sce <- as.SingleCellExperiment(liver)
immune.sce <- as.SingleCellExperiment(immune)
endothelia.sce <- as.SingleCellExperiment(endothelial)
MP.sce <- as.SingleCellExperiment(MP)
NKTcell.sce <- as.SingleCellExperiment(NKTcell)

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
# calculate cell type abundance per sample
n_cells <- table(liver.sce$Celltype, liver.sce$Patient.ID)
# calculate lineage abundance per sample
n_lineages <- table(liver.sce$Lineage, liver.sce$Patient.ID)

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
# visualise lineage composition at patient level
overall_lineage_patient.composition <- composition_barplot(df.lineages, df.lineages$Lineage, df.lineages$Sample.ID,
                                                           "Cell lineage", "Patient ID")
png(paste0(DA.path, "overall_lineage_patient.composition.barplot.png"), width=1200,height=600,units="px") # width = ncols X 100 + 100 for legend
print(overall_lineage_patient.composition)
dev.off()

# visualise lineage composition at condition level
overall_lineage_condition.composition <- composition_barplot(df.lineages, df.lineages$Lineage, df.lineages$Condition,
                                                             "Cell lineage", "Condition")
png(paste0(DA.path, "overall_lineage_condition.composition.barplot.png"), width=1200,height=600,units="px")
print(overall_lineage_condition.composition)
dev.off()

# compare composition plots
overall_lineage_composition <- overall_lineage_condition.composition + overall_lineage_patient.composition
png(paste0(DA.path, "overall_lineage_composition.barplot.png"), width=1200,height=600,units="px")
print(overall_lineage_composition)
dev.off()

# barplot of relative celltype abundances per sample
overall_celltype_abundance <- abundance_barplot(df.celltypes, df.celltypes$Sample.ID, df.celltypes$Celltype, "Patient ID", "Cell type")
png(paste0(DA.path, "overall_celltype_abundance.barplot.png"), width=1200,height=600,units="px")
print(overall_celltype_abundance)
dev.off()

# barplot of relative lineage abundances per sample
overall_lineage_abundance <- abundance_barplot(df.lineages, df.lineages$Sample.ID, df.lineages$Lineage, "Patient ID", "Cell lineage")
png(paste0(DA.path, "overall_lineage_abundance.barplot.png"), width=1100,height=600,units="px")
print(overall_lineage_abundance)
dev.off()

### DA analysis ###
# abundance table
abundances <- table(liver.sce$Lineage, liver.sce$Sample.ID) 
abundances <- unclass(abundances) 
head(abundances)

# attaching some column metadata
extra.info <- colData(liver.sce)[match(colnames(abundances), liver.sce$Sample.ID),]
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
topTags(res)
write.csv(res, file = paste0(DA.path, "overall_lineage_abundance.csv"))

##############
### IMMUNE ###
##############

### create meta data data frame ###
colData(immune.sce) %>% 
  as.data.frame %>% 
  transmute(
    Condition = condition, 
    Patient.ID = patient.ID,
    Sample.ID = patient.ID,
    Cluster.ID = seurat_clusters,
    Celltype = celltype_integrated,
    Lineage = lineage_integrated) %>%
  mutate_all(as.factor) %>% 
  set_rownames(colnames(immune.sce)) %>% 
  DataFrame -> colData(immune.sce)

# view data frame
head(colData(immune.sce))

# store number of cluster and number of samples
n.clusters <- length(kids <- set_names(levels(immune.sce$Celltype)[-c(7,8,10,11,13,15,16,17,18,21,24,25,27,31,33,34,36,38,39,40,42,43,44)]))
n.lineages <- length(lids <- set_names(levels(immune.sce$Lineage)))
n.samples <- length(sids <- set_names(levels(immune.sce$Patient.ID)))

# summary of experimental design
m <- match(sids, immune.sce$Sample.ID)
n_cells <- as.numeric(table(immune.sce$Sample.ID))

(ei <- data.frame(colData(immune.sce)[m, ], 
                  n_cells, row.names = NULL) %>% 
    select(-c("Celltype", "Lineage", "Cluster.ID")))

### DA visualisation ###
# calculate cell type abundance per sample
n_cells <- table(immune.sce$Celltype, immune.sce$Patient.ID)[-c(7,8,10,11,13,15,16,17,18,21,24,25,27,31,33,34,36,38,39,40,42,43,44),]
# calculate lineage abundance per sample
n_lineages <- table(immune.sce$Lineage, immune.sce$Patient.ID)

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
# visualise lineage composition at condition level
immune_lineage_patient.composition <- composition_barplot(df.lineages, df.lineages$Lineage, df.lineages$Sample.ID, "Cell lineage", "Patient ID")
png(paste0(DA.path, "immune_lineage_patient.composition.barplot.png"), width=800,height=600,units="px")
print(immune_lineage_patient.composition)
dev.off()

# visualise lineage composition at condition level
immune_lineage_condition.composition <- composition_barplot(df.lineages, df.lineages$Lineage, df.lineages$Condition, "Cell lineage", "Condition")
png(paste0(DA.path, "immune_lineage_condition.composition.barplot.png"), width=800,height=600,units="px")
print(immune_lineage_condition.composition)
dev.off()

# compare composition plots
immune_lineage_composition <- immune_lineage_condition.composition + immune_lineage_patient.composition
png(paste0(DA.path, "immune_lineage_composition.barplot.png"), width=800,height=600,units="px")
print(immune_lineage_composition)
dev.off()

# visualise celltype composition at condition level
immune_celltype_patient.composition <- composition_barplot(df.celltypes, df.celltypes$Celltype, df.celltypes$Sample.ID, "Cell types", "Patient ID")
png(paste0(DA.path, "immune_celltype_patient.composition.barplot.png"), width=800,height=600,units="px")
print(immune_celltype_patient.composition)
dev.off()

# visualise celltype composition at condition level
immune_celltype_condition.composition <- composition_barplot(df.celltypes, df.celltypes$Celltype, df.celltypes$Condition, "Cell type", "Condition")
png(paste0(DA.path, "immune_celltype_condition.composition.barplot.png"), width=800,height=600,units="px")
print(immune_celltype_condition.composition)
dev.off()

# compare composition plots
immune_celltype_composition <- immune_celltype_condition.composition + immune_celltype_patient.composition
png(paste0(DA.path, "immune_celltype_composition.barplot.png"), width=1200,height=600,units="px")
print(immune_celltype_composition)
dev.off()

# barplot of relative celltype abundances per sample
immune_celltype_abundance <- abundance_barplot(df.celltypes, df.celltypes$Sample.ID, df.celltypes$Celltype, "Patient ID", "Cell type")
png(paste0(DA.path, "immune_celltype_abundance.barplot.png"), width=1200,height=600,units="px")
print(immune_celltype_abundance)
dev.off()

# barplot of relative lineage abundances per sample
immune_lineage_abundance <- abundance_barplot(df.lineages, df.lineages$Sample.ID, df.lineages$Lineage, "Patient ID", "Cell lineage")
png(paste0(DA.path, "immune_lineage_abundance.barplot.png"), width=1100,height=600,units="px")
print(immune_lineage_abundance)
dev.off()

### DA analysis: cell level ###
# abundance table
abundances <- table(immune.sce$Celltype, immune.sce$Sample.ID) 
abundances <- unclass(abundances) 
head(abundances)

# attaching some column metadata.
extra.info <- colData(liver.sce)[match(colnames(abundances), liver.sce$Sample.ID),]
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
topTags(res)
write.csv(res, file = paste0(DA.path, "immmune_celltypes_abundance.csv"))

### DA analysis: lineage level ###
# abundance table
abundances <- table(immune.sce$Lineage, immune.sce$Sample.ID) 
abundances <- unclass(abundances) 
head(abundances)

# attaching some column metadata.
extra.info <- colData(liver.sce)[match(colnames(abundances), liver.sce$Sample.ID),]
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
topTags(res)
write.csv(res, file = paste0(DA.path, "immmune_lineage_abundance.csv"))

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
n.clusters <- length(kids <- set_names(levels(MP.sce$Celltype)[c(5,9,12,19,20,28,29,30,32)]))
n.samples <- length(sids <- set_names(levels(MP.sce$Patient.ID)))

# summary of experimental design
m <- match(sids, MP.sce$Sample.ID)
n_cells <- as.numeric(table(MP.sce$Sample.ID))

(ei <- data.frame(colData(MP.sce)[m, ], 
                  n_cells, row.names = NULL) %>% 
    select(-c("Celltype", "Lineage", "Cluster.ID")))

# calculate cell type abundance per sample
n_cells <- table(MP.sce$Celltype, MP.sce$Patient.ID)[c(5,9,12,19,20,28,29,30,32),]

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
# visualise celltype composition at patient level
MP_celltype_patient.composition <- composition_barplot(df.celltypes, df.celltypes$Celltype, df.celltypes$Sample.ID, "Cell type", "Patient ID")
png(paste0(DA.path, "MP_celltype_patient.composition.barplot.png"), width=1000,height=600,units="px")
print(MP_celltype_patient.composition)
dev.off()

# visualise celltype composition at condition level
MP_celltype_condition.composition <- composition_barplot(df.celltypes, df.celltypes$Celltype, df.celltypes$Condition, "Cell type", "Condition")
png(paste0(DA.path, "MP_celltype_condition.composition.barplot.png"), width=1000,height=600,units="px")
print(MP_celltype_condition.composition)
dev.off()

# compare composition plots
MP_celltype_composition <- MP_celltype_condition.composition + MP_celltype_patient.composition
png(paste0(DA.path, "MP_celltype_composition.barplot.png"), width=1000,height=600,units="px")
print(MP_celltype_composition)
dev.off()

# visualise cell type abundance
# barplot of relative  celltype abundances
MP_celltype_abundance <- abundance_barplot(df.celltypes, df.celltypes$Sample.ID, df.celltypes$Celltype, "Patient ID", "Cell type")
png(paste0(DA.path, "MP_celltype_abundance.barplot.png"), width=1100,height=600,units="px")
print(MP_celltype_abundance)
dev.off()

### DA analysis ###
# abundance table
abundances <- table(MP.sce$Celltype, MP.sce$Sample.ID) 
abundances <- unclass(abundances) 
head(abundances)

# attaching some column metadata.
extra.info <- colData(liver.sce)[match(colnames(abundances), liver.sce$Sample.ID),]
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
topTags(res)
write.csv(res, file = paste0(DA.path, "MP_celltype_abundance.csv"))

###############
### NKT CELLS ###
###############

### create meta data data frame ###
colData(NKTcell.sce) %>% 
  as.data.frame %>% 
  transmute(
    Condition = condition, 
    Patient.ID = patient.ID,
    Sample.ID = patient.ID,
    Cluster.ID = seurat_clusters,
    Celltype = celltype_integrated,
    Lineage = lineage_integrated) %>%
  mutate_all(as.factor) %>% 
  set_rownames(colnames(NKTcell.sce)) %>% 
  DataFrame -> colData(NKTcell.sce)

# view data frame
head(colData(NKTcell.sce))

# store number of cluster and number of samples
n.clusters <- length(kids <- set_names(levels(NKTcell.sce$Celltype)[c(1,2,3,4,6,23,26,41,45)]))
n.samples <- length(sids <- set_names(levels(NKTcell.sce$Patient.ID)))

# summary of experimental design
m <- match(sids, NKTcell.sce$Sample.ID)
n_cells <- as.numeric(table(NKTcell.sce$Sample.ID))

(ei <- data.frame(colData(NKTcell.sce)[m, ], 
                  n_cells, row.names = NULL) %>% 
    select(-c("Celltype", "Lineage", "Cluster.ID")))

### DA visualisation ###
# calculate cell type abundance per sample
n_cells <- table(NKTcell.sce$Celltype, NKTcell.sce$Patient.ID)[c(1,2,3,4,6,23,26,41,45),]

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
# visualise celltype composition at patient level
NKTcell_celltype_patient.composition <- composition_barplot(df.celltypes, df.celltypes$Celltype, df.celltypes$Sample.ID, "Cell type", "Patient ID")
png(paste0(DA.path, "NKTcell_celltype_patient.composition.barplot.png"), width=1000,height=600,units="px")
print(NKTcell_celltype_patient.composition)
dev.off()

# visualise celltype composition at condition level
NKTcell_celltype_condition.composition <- composition_barplot(df.celltypes, df.celltypes$Celltype, df.celltypes$Condition, "Cell type", "Condition")
png(paste0(DA.path, "NKTcell_celltype_condition.composition.barplot.png"), width=1000,height=600,units="px")
print(NKTcell_celltype_condition.composition)
dev.off()

# compare composition plots
NKTcell_celltype_composition <- NKTcell_celltype_condition.composition + NKTcell_celltype_patient.composition
png(paste0(DA.path, "NKTcell_celltype_composition.barplot.png"), width=1000,height=600,units="px")
print(NKTcell_celltype_composition)
dev.off()

# barplot of relative  celltype abundances
NKTcell_celltype_abundance <- abundance_barplot(df.celltypes, df.celltypes$Sample.ID, df.celltypes$Celltype, "Patient ID", "Cell type")
png(paste0(DA.path, "NKTcell_celltype_abundance.barplot.png"), width=1100,height=600,units="px")
print(NKTcell_celltype_abundance)
dev.off()

### DA analysis ###
# abundance table
abundances <- table(NKTcell.sce$Celltype, NKTcell.sce$Sample.ID) 
abundances <- unclass(abundances) 
head(abundances)

# attaching some column metadata.
extra.info <- colData(liver.sce)[match(colnames(abundances), liver.sce$Sample.ID),]
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
topTags(res)
write.csv(res, file = paste0(DA.path, "NKTcell_celltype_abundance.csv"))

###################
### ENDOTHELIAL ###
###################

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
n.clusters <- length(kids <- set_names(levels(endothelia.sce$Celltype)[c(8,11,13,15,17,18,24,27)]))
n.samples <- length(sids <- set_names(levels(endothelia.sce$Patient.ID)))

# summary of experimental design
m <- match(sids, endothelia.sce$Sample.ID)
n_cells <- as.numeric(table(endothelia.sce$Sample.ID))

(ei <- data.frame(colData(endothelia.sce)[m, ], 
                  n_cells, row.names = NULL) %>% 
    select(-c("Celltype", "Lineage", "Cluster.ID")))

### DA visualisation ###
# calculate cell type abundance per sample
n_cells <- table(endothelia.sce$Celltype, endothelia.sce$Patient.ID)[c(8,11,13,15,17,18,24,27),]

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
# visualise celltype composition at patient level
endothelia_celltype_patient.composition <- composition_barplot(df.celltypes, df.celltypes$Celltype, df.celltypes$Sample.ID, "Cell type", "Patient ID")
png(paste0(DA.path, "endothelia_celltype_patient.composition.barplot.png"), width=1000,height=600,units="px")
print(endothelia_celltype_patient.composition)
dev.off()

# visualise celltype composition at condition level
endothelia_celltype_condition.composition <- composition_barplot(df.celltypes, df.celltypes$Celltype, df.celltypes$Condition, "Cell type", "Condition")
png(paste0(DA.path, "endothelia_celltype_condition.composition.barplot.png"), width=1000,height=600,units="px")
print(endothelia_celltype_condition.composition)
dev.off()

# compare composition plots
endothelia_celltype_composition <- endothelia_celltype_condition.composition + endothelia_celltype_patient.composition
png(paste0(DA.path, "endothelia_celltype_composition.barplot.png"), width=1000,height=600,units="px")
print(endothelia_celltype_composition)
dev.off()

# barplot of relative  celltype abundances
endothelia_celltype_abundance <- abundance_barplot(df.celltypes, df.celltypes$Sample.ID, df.celltypes$Celltype, "Patient ID", "Cell type")
png(paste0(DA.path, "endothelia_celltype_abundance.barplot.png"), width=1100,height=600,units="px")
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
topTags(res)
write.csv(res, file = paste0(DA.path, "endothelia_celltype_abundance.csv"))

##################
### MESENCHYME ###
##################

#################
### EPITHELIA ###
#################