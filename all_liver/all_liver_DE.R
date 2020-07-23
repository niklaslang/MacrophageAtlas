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
DE.path <- "/home/s1987963/ds_group/Niklas/all_liver/DE/"

### read files ###
liver <- readRDS(liver.path)

### subset data ###
# only use ramachandran data for downstream analyses
liver <- subset(liver, subset = study == "ramachandran_liver") # 82132 cells

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

### overview ###
table(liver.sce$celltype_integrated, liver.sce$condition)

# creating pseudo-bulk samples
liver.pseudobulk <- aggregateAcrossCells(liver.sce, id=colData(liver.sce)[,c("celltype_integrated", "patient.ID")])
liver.pseudobulk

# picking one cell type
label <- "Combined Liver cDC2"
pseudobulk.cell <- liver.pseudobulk[,label==liver.pseudobulk$celltype_integrated]

# Creating up a DGEList object for use in edgeR:
y <- DGEList(counts(pseudobulk.cell), samples=colData(pseudobulk.cell))
y

# filtering genes 
keep <- filterByExpr(y, group=pseudobulk.cell$condition)
y <- y[keep,]
summary(keep)

# compute normalization factors to account for composition bias
y <- calcNormFactors(y)
y$samples

# statistical modelling
# create design matrix
design <- model.matrix(~factor(condition), y$samples)
design

# estimate NP dispersion
y <- estimateDisp(y, design)
summary(y$trended.dispersion)
plotBCV(y)

# estimate the quasi-likelihood dispersion
fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$var.prior)
summary(fit$df.prior)
plotQLDisp(fit)

# test for differences 
res <- glmQLFTest(fit, coef=ncol(design))
res.dt <- data.table(res$table)
summary(decideTests(res))
top.genes <- topTags(res, n=100)
write.csv(top.genes$table, file = paste0(DE.path, "5_7/DE_combined_liver_cDC2.csv"))

### ALLTOGETHER ###
# Pulling out a sample-level 'targets' data.frame:
targets <- colData(liver.sce)[!duplicated(liver.sce$patient.ID),]

# Constructing the design matrix:
design <-  model.matrix(~ factor(condition), data=targets)
rownames(design) <- targets$patient.ID

# run scran's pseudoBulkDGE()
de.results <- pseudoBulkDGE(liver.pseudobulk, 
                            sample=liver.pseudobulk$patient.ID,
                            label=liver.pseudobulk$celltype_integrated,
                            design=design,
                            coef=ncol(design),
                            condition=targets$condition)
