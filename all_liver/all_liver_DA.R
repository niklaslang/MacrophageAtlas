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

### convert seurat object to sce ###
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
m1 <- match(sids, liver.sce$Sample.ID)
m2 <- match(lids, liver.sce$Lineage)
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

# visualise cell type abundance
# prep data.frame for plotting
df <- data.frame(
  Frequency = as.numeric(freqs_cells), 
  Celltype = rep(kids, n.samples),
  Sample.ID = rep(sids, each = n.clusters))
m <- match(df$Sample.ID, ei$Sample.ID)
df$Condition <- ei$Condition[m]

# barplot of relative  celltype abundances
cell_abundance <- ggplot(df, aes(x = Sample.ID, y = Frequency, fill = Celltype)) +
  geom_bar(stat = "identity", color="black", position = "fill") +
  facet_wrap(~ Condition, scales = "free_x") +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12), 
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size=16))
png(paste0(DA.path, "cell.barplot.png"), width=1800,height=600,units="px")
print(cell_abundance)
dev.off()

# boxplot of relative cluster-abundances
cell_abundance2 <- ggplot(df, aes(x = Condition, y = Frequency, color = Condition)) +
  geom_boxplot(outlier.colour = NA) +  geom_jitter() +
  facet_wrap(~ Celltype, scales = "free_y", ncol = 4) +
  theme_classic()
png(paste0(DA.path, "cell.boxplot.png"), width=1800,height=1200,units="px")
print(cell_abundance2)
dev.off()

# visualise lineage abundance
# prep data.frame for plotting
df <- data.frame(
  Frequency = as.numeric(freqs_lineages), 
  Lineage = rep(lids, n.samples),
  Sample.ID = rep(sids, each = n.lineages))
m <- match(df$Sample.ID, ei$Sample.ID)
df$Condition <- ei$Condition[m]

# barplot of relative cluster-abundances
lineage_abundance <- ggplot(df, aes(x = Sample.ID, y = Frequency, fill = Lineage)) +
  geom_bar(stat = "identity", color="black", position = "fill") + 
  facet_wrap(~ Condition, scales = "free_x") +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12), 
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size=16))
png(paste0(DA.path, "lineage.histogram_persample.png"), width=1200,height=600,units="px")
print(lineage_abundance)
dev.off()

# boxplot of relative cluster-abundances
lineage_abundance2 <- ggplot(df, aes(x = Condition, y = Frequency, color = Condition)) +
  geom_boxplot(outlier.colour = NA) +  geom_jitter() +
  facet_wrap(~ Lineage, scales = "free_y", ncol = 4) +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12), 
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size=16))
png(paste0(DA.path, "lineage.boxplot.png"), width=1800,height=1200,units="px")
print(lineage_abundance2)
dev.off()

### DA analysis ###
# abundance table
abundances <- table(liver.sce$Celltype, liver.sce$Sample.ID) 
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

# assuming most labels do not change
y.ab2 <- calcNormFactors(y.ab)
y.ab2$samples$norm.factors
y.ab2 <- estimateDisp(y.ab2, design, trend="none")
fit.ab2 <- glmQLFit(y.ab2, design, robust=TRUE, abundance.trend=FALSE)
res2 <- glmQLFTest(fit.ab2, coef=ncol(design))
summary(decideTests(res2))
topTags(res2, n=10)

# removing epithelial, endothelial, mesenchymal and proliferating cells
offenders <- c("Combined Liver Hepatocyte 1", "Combined Liver Cholangiocyte 1",
               "Combined Liver Proliferating 1",  "Combined Liver Hepatocyte 2",
               "Combined Liver Cholangiocyte 2", "Combined Liver Cholangiocyte 3",
               "Combined Liver LSEndo 1", "Combined Liver Endo 1", "Combined Liver PV Endo",
               "Combined Liver Endo 2", "Combined Liver HA Endo", "Combined Liver CV Endo",
               "Combined Liver LSEndo 2", "Combined Liver LymphEndo", "Combined Liver Mesenchyme 1",
               "Combined Liver Mesenchyme 2", "Combined Liver Mesenchyme 3", "Combined Liver Proliferating 1",
               "Combined Liver Proliferating 2", "Combined Liver Proliferating 3")
y.ab3 <- y.ab[setdiff(rownames(y.ab), offenders), keep.lib.sizes=FALSE]
y.ab3$samples
y.ab3 <- estimateDisp(y.ab3, design, trend="none")
fit.ab3 <- glmQLFit(y.ab3, design, robust=TRUE, abundance.trend=FALSE)
res3 <- glmQLFTest(fit.ab3, coef=ncol(design))
summary(decideTests(res3))
topTags(res3, n=10)

# testing against a log-fold change threshold
res.lfc <- glmTreat(fit.ab, coef=ncol(design), lfc=1)
summary(decideTests(res.lfc))
topTags(res.lfc)
