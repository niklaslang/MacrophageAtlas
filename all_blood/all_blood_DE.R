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
blood.path <- "/home/s1987963/ds_group/Niklas/combined_organs/combined_blood_annotated.rds"
DE.path <- "/home/s1987963/ds_group/Niklas/all_blood/DE/"

### read files ###
blood <- readRDS(blood.path)

### convert to sce ###
blood.sce <- as.SingleCellExperiment(blood)

### overview ###
table(blood.sce$celltype_integrated, blood.sce$condition)

# creating pseudo-bulk samples
blood.pseudobulk <- aggregateAcrossCells(blood.sce, id=colData(blood.sce)[,c("celltype_integrated", "patient.ID")])

# picking one cell type
label1 <- "Combined Blood CD14+ Monocyte 1"
label2 <- "Combined Blood CD14+ Monocyte 2"
label3 <- "Combined Blood CD14+ Monocyte 3"
label4 <- "Combined Blood CD16+ Monocyte"
label5 <- "Combined Blood cDC1"
label6 <- "Combined Blood cDC2"
label7 <- "Combined Blood pDC 1"
label8 <- "Combined Blood pDC 2"

##############################
### no.1: CD14+ Monocyte 1 ###
##############################

## prepare data ##
# subset data
pseudobulk.cell <- blood.pseudobulk[,label1==blood.pseudobulk$celltype_integrated]

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
summary(decideTests(res)) #  1071 down and 606 up
top.genes <- topTags(res, n=2000)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_blood_cd14+_monocyte_1.csv"))

##############################
### no.2: CD14+ Monocyte 2 ###
##############################

## prepare data ##
# subset data
pseudobulk.cell <- blood.pseudobulk[,label2==blood.pseudobulk$celltype_integrated]

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
summary(decideTests(res)) # 866 down and 570 up
top.genes <- topTags(res, n=1500)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_blood_cd14+_monocyte_2.csv"))

##############################
### no.3: CD14+ Monocyte 3 ###
##############################

## prepare data ##
# subset data
pseudobulk.cell <- blood.pseudobulk[,label3==blood.pseudobulk$celltype_integrated]

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
summary(decideTests(res)) # 67 down and 51 up
top.genes <- topTags(res, n=1000)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_blood_cd14+_monocyte_3.csv"))

############################
### no.4: CD16+ Monocyte ###
############################

## prepare data ##
# subset data
pseudobulk.cell <- blood.pseudobulk[,label4==blood.pseudobulk$celltype_integrated]

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
summary(decideTests(res)) # 605 down and 509 up
top.genes <- topTags(res, n=1500)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_blood_cd16+_monocyte.csv"))

###################
### no.5: cDC 1 ###
###################

## prepare data ##
# subset data
pseudobulk.cell <- blood.pseudobulk[,label5==blood.pseudobulk$celltype_integrated]

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
summary(decideTests(res)) # 96 down and 193 up
top.genes <- topTags(res, n=1000)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_blood_cDC1.csv"))

##################
### no.6: cDC2 ###
##################

## prepare data ##
# subset data
pseudobulk.cell <- blood.pseudobulk[,label6==blood.pseudobulk$celltype_integrated]

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
summary(decideTests(res)) # 538 down and 557 up
top.genes <- topTags(res, n=1500)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_blood_cDC2.csv"))

###################
### no.7: pDC 1 ###
###################

## prepare data ##
# subset data
pseudobulk.cell <- blood.pseudobulk[,label7==blood.pseudobulk$celltype_integrated]

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
summary(decideTests(res)) # 468 down and 354 up
top.genes <- topTags(res, n=1500)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_blood_pDC_1.csv"))

###################
### no.8: pDC 2 ###
###################

## prepare data ##
# subset data
pseudobulk.cell <- blood.pseudobulk[,label8==blood.pseudobulk$celltype_integrated]

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
summary(decideTests(res)) # 45 down and 151 up
top.genes <- topTags(res, n=1000)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_blood_pDC_2.csv"))
