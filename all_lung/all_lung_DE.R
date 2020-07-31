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
lung.path <- "/home/s1987963/ds_group/Niklas/combined_organs/combined_lung_annotated.rds"
DE.path <- "/home/s1987963/ds_group/Niklas/all_lung/DE/"

### read files ###
lung <- readRDS(lung.path)

### convert to sce ###
lung.sce <- as.SingleCellExperiment(lung)

### overview ###
table(lung.sce$celltype_integrated, lung.sce$condition)

# creating pseudo-bulk samples
lung.pseudobulk <- aggregateAcrossCells(lung.sce, id=colData(lung.sce)[,c("celltype_integrated", "patient.ID")])

# picking one cell type
label1 <- "Combined Lung Monocyte"
label2 <- "Combined Lung Macrophage 1"
label3 <- "Combined Lung Macrophage 2"
label4 <- "Combined Lung Macrophage 3"
label5 <- "Combined Lung Macrophage 4"
label6 <- "Combined Lung IFN-primed Macrophage"
label7 <- "Combined Lung Alveolar Macrophage 1"
label8 <- "Combined Lung Alveolar Macrophage 2"
label9 <- "Combined Lung Alveolar Macrophage 3"
label10 <- "Combined Lung Alveolar Macrophage 4"
label11 <- "Combined Lung cDC2"
label12 <- "Combined Lung pDC"
label.list <- c(label1, label2, label3, label4, label5, label6, label7, label8, label9, label10, label11, label12)

#######################
### no.1: monocytes ###
#######################

## prepare data ##
# subset data
pseudobulk.cell <- lung.pseudobulk[,label1==lung.pseudobulk$celltype_integrated]

# Creating up a DGEList object for use in edgeR:
y <- DGEList(counts(pseudobulk.cell), samples=colData(pseudobulk.cell))

# filtering genes 
keep <- filterByExpr(y, group=pseudobulk.cell$condition)
y <- y[keep,]

# compute normalization factors to account for composition bias
y <- calcNormFactors(y)

# statistical modelling ##
# design matrix
design <- model.matrix(~factor(condition), y$samples)

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
summary(decideTests(res)) # 991 down and 551 up
top.genes <- topTags(res, n=3000)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_lung_monocyte.csv"))

##########################
### no.2: macrophage 1 ###
##########################

## prepare data ##
# subset data
pseudobulk.cell <- lung.pseudobulk[,label2==lung.pseudobulk$celltype_integrated]

# Creating up a DGEList object for use in edgeR:
y <- DGEList(counts(pseudobulk.cell), samples=colData(pseudobulk.cell))

# filtering genes 
keep <- filterByExpr(y, group=pseudobulk.cell$condition)
y <- y[keep,]

# compute normalization factors to account for composition bias
y <- calcNormFactors(y)

# statistical modelling ##
# design matrix
design <- model.matrix(~factor(condition), y$samples)

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
summary(decideTests(res)) # 645 down and 492 up
top.genes <- topTags(res, n=3000)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_lung_macrophage_1.csv"))

###########################
### no.3: macrophages 2 ###
###########################

## prepare data ##
# subset data
pseudobulk.cell <- lung.pseudobulk[,label3==lung.pseudobulk$celltype_integrated]

# Creating up a DGEList object for use in edgeR:
y <- DGEList(counts(pseudobulk.cell), samples=colData(pseudobulk.cell))

# filtering genes 
keep <- filterByExpr(y, group=pseudobulk.cell$condition)
y <- y[keep,]

# compute normalization factors to account for composition bias
y <- calcNormFactors(y)

# statistical modelling ##
# design matrix
design <- model.matrix(~factor(condition), y$samples)

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
summary(decideTests(res)) # 1233 down and 870 up
top.genes <- topTags(res, n=3000)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_lung_macrophage_2.csv"))

###########################
### no.4: macrophages 3 ###
###########################

## prepare data ##
# subset data
pseudobulk.cell <- lung.pseudobulk[,label4==lung.pseudobulk$celltype_integrated]

# Creating up a DGEList object for use in edgeR:
y <- DGEList(counts(pseudobulk.cell), samples=colData(pseudobulk.cell))

# filtering genes 
keep <- filterByExpr(y, group=pseudobulk.cell$condition)
y <- y[keep,]

# compute normalization factors to account for composition bias
y <- calcNormFactors(y)

# statistical modelling ##
# design matrix
design <- model.matrix(~factor(condition), y$samples)

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
summary(decideTests(res)) # 256 down and 261 up
top.genes <- topTags(res, n=3000)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_lung_macrophage_3.csv"))

###########################
### no.5: macrophages 4 ###
###########################

## prepare data ##
# subset data
pseudobulk.cell <- lung.pseudobulk[,label5==lung.pseudobulk$celltype_integrated]

# Creating up a DGEList object for use in edgeR:
y <- DGEList(counts(pseudobulk.cell), samples=colData(pseudobulk.cell))

# filtering genes 
keep <- filterByExpr(y, group=pseudobulk.cell$condition)
y <- y[keep,]

# compute normalization factors to account for composition bias
y <- calcNormFactors(y)

# statistical modelling ##
# design matrix
design <- model.matrix(~factor(condition), y$samples)

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
summary(decideTests(res)) # 61 down and 98 up
top.genes <- topTags(res, n=3000)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_lung_macrophage_4.csv"))

###################################
### no.6: IFN-primed Macrophage ###
###################################

## prepare data ##
# subset data
pseudobulk.cell <- lung.pseudobulk[,label6==lung.pseudobulk$celltype_integrated]

# Creating up a DGEList object for use in edgeR:
y <- DGEList(counts(pseudobulk.cell), samples=colData(pseudobulk.cell))

# filtering genes 
keep <- filterByExpr(y, group=pseudobulk.cell$condition)
y <- y[keep,]

# compute normalization factors to account for composition bias
y <- calcNormFactors(y)

# statistical modelling ##
# design matrix
design <- model.matrix(~factor(condition), y$samples)

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
summary(decideTests(res)) # 9 down and 12 up
top.genes <- topTags(res, n=700)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_lung_IFN_primed_macrophage.csv"))

###################################
### no.7: alveolar macrophage 1 ###
###################################

## prepare data ##
# subset data
pseudobulk.cell <- lung.pseudobulk[,label7==lung.pseudobulk$celltype_integrated]

# Creating up a DGEList object for use in edgeR:
y <- DGEList(counts(pseudobulk.cell), samples=colData(pseudobulk.cell))

# filtering genes 
keep <- filterByExpr(y, group=pseudobulk.cell$condition)
y <- y[keep,]

# compute normalization factors to account for composition bias
y <- calcNormFactors(y)

# statistical modelling ##
# design matrix
design <- model.matrix(~factor(condition), y$samples)

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
summary(decideTests(res)) # 1484 down and 702 up
top.genes <- topTags(res, n=3000)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_lung_alveolar_macrophage_1.csv"))

###################################
### no.8: alveolar macrophage 2 ###
###################################

## prepare data ##
# subset data
pseudobulk.cell <- lung.pseudobulk[,label8==lung.pseudobulk$celltype_integrated]

# Creating up a DGEList object for use in edgeR:
y <- DGEList(counts(pseudobulk.cell), samples=colData(pseudobulk.cell))

# filtering genes 
keep <- filterByExpr(y, group=pseudobulk.cell$condition)
y <- y[keep,]

# compute normalization factors to account for composition bias
y <- calcNormFactors(y)

# statistical modelling ##
# design matrix
design <- model.matrix(~factor(condition), y$samples)

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
summary(decideTests(res)) # 813 down and 625 up
top.genes <- topTags(res, n=1500)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_lung_alveolar_macrophage_2.csv"))

###################################
### no.9: alveolar macrophage 3 ###
###################################

## prepare data ##
# subset data
pseudobulk.cell <- lung.pseudobulk[,label9==lung.pseudobulk$celltype_integrated]

# Creating up a DGEList object for use in edgeR:
y <- DGEList(counts(pseudobulk.cell), samples=colData(pseudobulk.cell))

# filtering genes 
keep <- filterByExpr(y, group=pseudobulk.cell$condition)
y <- y[keep,]

# compute normalization factors to account for composition bias
y <- calcNormFactors(y)

# statistical modelling ##
# design matrix
design <- model.matrix(~factor(condition), y$samples)

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
summary(decideTests(res)) # 369 down and 302 up
top.genes <- topTags(res, n=1000)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_lung_alveolar_macrophage_3.csv"))

###################################
### no.10: alveolar macrophage 4 ###
###################################

## prepare data ##
# subset data
pseudobulk.cell <- lung.pseudobulk[,label10==lung.pseudobulk$celltype_integrated]

# Creating up a DGEList object for use in edgeR:
y <- DGEList(counts(pseudobulk.cell), samples=colData(pseudobulk.cell))

# filtering genes 
keep <- filterByExpr(y, group=pseudobulk.cell$condition)
y <- y[keep,]

# compute normalization factors to account for composition bias
y <- calcNormFactors(y)

# statistical modelling ##
# design matrix
design <- model.matrix(~factor(condition), y$samples)

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
summary(decideTests(res)) # 162 down and  159 up
top.genes <- topTags(res, n=1000)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_lung_alveolar_macrophage_4.csv"))

###################
### no.11: cDC2 ###
###################

## prepare data ##
# subset data
pseudobulk.cell <- lung.pseudobulk[,label11==lung.pseudobulk$celltype_integrated]

# Creating up a DGEList object for use in edgeR:
y <- DGEList(counts(pseudobulk.cell), samples=colData(pseudobulk.cell))

# filtering genes 
keep <- filterByExpr(y, group=pseudobulk.cell$condition)
y <- y[keep,]

# compute normalization factors to account for composition bias
y <- calcNormFactors(y)

# statistical modelling ##
# design matrix
design <- model.matrix(~factor(condition), y$samples)

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
summary(decideTests(res)) # 1155 down and 829 up
top.genes <- topTags(res, n=2500)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_lung_cDC2.csv"))

##################
### no.12: pDC ###
##################

## prepare data ##
# subset data
pseudobulk.cell <- lung.pseudobulk[,label12==lung.pseudobulk$celltype_integrated]

# Creating up a DGEList object for use in edgeR:
y <- DGEList(counts(pseudobulk.cell), samples=colData(pseudobulk.cell))

# filtering genes 
keep <- filterByExpr(y, group=pseudobulk.cell$condition)
y <- y[keep,]

# compute normalization factors to account for composition bias
y <- calcNormFactors(y)

# statistical modelling ##
# design matrix
design <- model.matrix(~factor(condition), y$samples)

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
summary(decideTests(res)) # 11 down and 11 up
top.genes <- topTags(res, n=1000)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_lung_pDC.csv"))
