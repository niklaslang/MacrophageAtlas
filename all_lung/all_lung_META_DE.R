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
library(EnhancedVolcano)

### path variables ###
lung.path <- "/home/s1987963/ds_group/Niklas/combined_organs/combined_lung_annotated.rds"
DE.path <- "/home/s1987963/ds_group/Niklas/all_lung/DE/"

### read files ###
lung <- readRDS(lung.path)

### create meta cell types ###
cell.data <- data.table(barcode = colnames(lung),
                        celltype = lung$celltype_integrated)
cell.data[, meta.celltype := ifelse(celltype %in% c("Combined Lung Macrophage 1",
                                                    "Combined Lung Macrophage 2",
                                                    "Combined Lung Macrophage 3",
                                                    "Combined Lung Macrophage 4",
                                                    "Combined Lung Macrophage 5",
                                                    "Combined Lung IFN-primed Macrophage"),
                                    "Combined Lung monocyte-derived Macrophage", paste0(celltype))]
cell.data[, meta.celltype := ifelse(celltype %in% c("Combined Lung Alveolar Macrophage 1",
                                                    "Combined Lung Alveolar Macrophage 2", 
                                                    "Combined Lung Alveolar Macrophage 3",
                                                    "Combined Lung Alveolar Macrophage 4"),
                                    "Combined Lung tissue resident Macrophage", paste0(meta.celltype))]
cell.data[, meta.celltype := ifelse(celltype %in% c("Combined Lung cDC2"),
                                    "Combined Lung cDC", paste0(meta.celltype))] 
cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$barcode <- NULL
cell.data$celltype <- NULL
lung <- AddMetaData(lung, cell.data, col.name = "meta.celltype")

### convert to sce ###
lung.sce <- as.SingleCellExperiment(lung)

### overview ###
table(lung.sce$meta.celltype, lung.sce$condition)

# creating pseudo-bulk samples
lung.pseudobulk <- aggregateAcrossCells(lung.sce, id=colData(lung.sce)[,c("meta.celltype", "patient.ID")])

# label meta cell types
label1 <- "Combined Lung monocyte-derived Macrophage"
label2 <- "Combined Lung tissue resident Macrophage"
label3 <- "Combined Lung Monocyte"
label4 <- "Combined Lung cDC"

##########################################
### no.1: monocyte-derived Macrophage  ###
##########################################

## prepare data ##
# subset data
pseudobulk.cell <- lung.pseudobulk[,label1==lung.pseudobulk$meta.celltype]

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
summary(decideTests(res)) #  2085 down and 1698 up
top.genes <- topTags(res, n=9265)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_lung_META_monocyte_derived_macrophages.csv"))

# save results in data frame 
monocyte_derived_macrophage.genes <- data.frame(gene = rownames(top.genes$table),
                                     p.adj = top.genes$table$FDR,
                                     logFC = top.genes$table$logFC)

####################
### volcano plot ###
####################
# prepare data #
p.adj.cutoff <- 0.05
logFC.cutoff <- 1.0

monocyte_derived_macrophage.genes$significance <- "NS"
monocyte_derived_macrophage.genes$significance[(abs(monocyte_derived_macrophage.genes$logFC) > logFC.cutoff)] <- "FC"
monocyte_derived_macrophage.genes$significance[(monocyte_derived_macrophage.genes$p.adj < p.adj.cutoff)] <- "FDR"
monocyte_derived_macrophage.genes$significance[(monocyte_derived_macrophage.genes$p.adj < p.adj.cutoff) & (abs(monocyte_derived_macrophage.genes$logFC) > logFC.cutoff)] <- "FC_FDR"
table(monocyte_derived_macrophage.genes$significance)
monocyte_derived_macrophage.genes$significance <- factor(monocyte_derived_macrophage.genes$significance, levels=c("NS", "FC", "FDR", "FC_FDR"))
# sort genes by ordered padj
monocyte_derived_macrophage.genes_ordered <- monocyte_derived_macrophage.genes[order(monocyte_derived_macrophage.genes$p.adj), ]
# create a column to indicate which genes to label
monocyte_derived_macrophage.genes_ordered$gene.labels <- ""
monocyte_derived_macrophage.genes_ordered$gene.labels <- monocyte_derived_macrophage.genes_ordered$gene %in% monocyte_derived_macrophage.genes_ordered$gene[1:100] 

# volcano plot
monocyte_derived_macrophage.plot1 <- ggplot(monocyte_derived_macrophage.genes_ordered, aes(x=logFC, y=-log10(p.adj))) +
  geom_point(aes(color=factor(significance)), alpha=1/2, size=2.0) +
  geom_text_repel(aes(x=logFC, y=-log10(p.adj), label = ifelse(gene.labels == T, paste0(gene),""))) +
  scale_color_manual( values=c(NS="grey30", FC="forestgreen", FDR="royalblue", FC_FDR="red2"), 
                      labels=c(NS="NS", 
                               FC="logFC", 
                               FDR="p-value (adj)", 
                               FC_FDR="logFC & p-value (adj)")) +
  scale_x_continuous(limits = c(-7.5,7.5)) +
  scale_y_continuous(limits = c(0,12.5)) +
  ggtitle("Lung monocyte-derived macrophage DE genes\nhealthy vs. fibrotic") +
  xlab(bquote(~Log[2]~ "fold change")) +
  ylab(bquote(~-Log[10]~adjusted~italic(P))) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.position="top",            #Moves the legend to the top of the plot
        legend.key=element_blank(),       #removes the border
        legend.key.size=unit(0.5, "cm"),  #Sets overall area/size of the legend
        legend.text=element_text(size=12), #Text size
        legend.title=element_blank()) +
  geom_vline(xintercept=c(-logFC.cutoff, logFC.cutoff), linetype="longdash", colour="black", size=0.4) +
  geom_hline(yintercept=-log10(p.adj.cutoff), linetype="longdash", colour="black", size=0.4)
png(paste0(DE.path, "lung_META_monocyte_derived_macrophage_volcano_1.png"), width=1200,height=800,units="px")
print(monocyte_derived_macrophage.plot1)
dev.off()

#############################
### enhanced volcano plot ###
#############################
monocyte_derived_macrophage.plot2 <- EnhancedVolcano(top.genes$table,
                                          lab = rownames(top.genes$table),
                                          x = 'logFC',
                                          y = 'FDR',
                                          pCutoff = 0.05,
                                          FCcutoff = 1,
                                          ylim = c(0,12.5),
                                          xlim = c(-7.5,7.5),
                                          title = NULL,#'DE genes in kidney monocytes',
                                          subtitle = NULL,#'healthy versus fibrotic',
                                          cutoffLineType = 'twodash',
                                          cutoffLineWidth = 0.8,
                                          pointSize = 2,
                                          labSize = 5.0,
                                          colAlpha = 0.6,
                                          legendLabels=c('NS','logFC','p-value',
                                                         'p-value & logFC'),
                                          legendPosition = 'right',
                                          legendLabSize = 10,
                                          legendIconSize = 2.0,
                                          gridlines.major = FALSE,
                                          gridlines.minor = FALSE)
png(paste0(DE.path, "lung_META_monocyte_derived_macrophage_volcano_2.png"), width=1200,height=800,units="px")
print(monocyte_derived_macrophage.plot2)
dev.off()

#########################################
### no.2: tissue resident Macrophage  ###
#########################################

## prepare data ##
# subset data
pseudobulk.cell <- lung.pseudobulk[,label2==lung.pseudobulk$meta.celltype]

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
summary(decideTests(res)) # 2149 down and 1550 up
top.genes <- topTags(res, n=10472)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_lung_META_tissue_resident_macrophage.csv"))

# save results in data frame 
tissue_resident_macrophage.genes <- data.frame(gene = rownames(top.genes$table),
                                     p.adj = top.genes$table$FDR,
                                     logFC = top.genes$table$logFC)

####################
### volcano plot ###
####################
# prepare data #
p.adj.cutoff <- 0.05
logFC.cutoff <- 1.0

tissue_resident_macrophage.genes$significance <- "NS"
tissue_resident_macrophage.genes$significance[(abs(tissue_resident_macrophage.genes$logFC) > logFC.cutoff)] <- "FC"
tissue_resident_macrophage.genes$significance[(tissue_resident_macrophage.genes$p.adj < p.adj.cutoff)] <- "FDR"
tissue_resident_macrophage.genes$significance[(tissue_resident_macrophage.genes$p.adj < p.adj.cutoff) & (abs(tissue_resident_macrophage.genes$logFC) > logFC.cutoff)] <- "FC_FDR"
table(tissue_resident_macrophage.genes$significance)
tissue_resident_macrophage.genes$significance <- factor(tissue_resident_macrophage.genes$significance, levels=c("NS", "FC", "FDR", "FC_FDR"))
# sort genes by ordered padj
tissue_resident_macrophage.genes_ordered <- tissue_resident_macrophage.genes[order(tissue_resident_macrophage.genes$p.adj), ]
# create a column to indicate which genes to label
tissue_resident_macrophage.genes_ordered$gene.labels <- ""
tissue_resident_macrophage.genes_ordered$gene.labels <- tissue_resident_macrophage.genes_ordered$gene %in% tissue_resident_macrophage.genes_ordered$gene[1:100] 

# volcano plot
tissue_resident_macrophage.plot1 <- ggplot(tissue_resident_macrophage.genes_ordered, aes(x=logFC, y=-log10(p.adj))) +
  geom_point(aes(color=factor(significance)), alpha=1/2, size=2.0) +
  geom_text_repel(aes(x=logFC, y=-log10(p.adj), label = ifelse(gene.labels == T, paste0(gene),""))) +
  scale_color_manual( values=c(NS="grey30", FC="forestgreen", FDR="royalblue", FC_FDR="red2"), 
                      labels=c(NS="NS", 
                               FC="logFC", 
                               FDR="p-value (adj)", 
                               FC_FDR="logFC & p-value (adj)")) +
  scale_x_continuous(limits = c(-7.5,7.5)) +
  scale_y_continuous(limits = c(0,10)) +
  ggtitle("Lung tissue resident macrophages DE genes\nhealthy vs. fibrotic") +
  xlab(bquote(~Log[2]~ "fold change")) +
  ylab(bquote(~-Log[10]~adjusted~italic(P))) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.position="top",            #Moves the legend to the top of the plot
        legend.key=element_blank(),       #removes the border
        legend.key.size=unit(0.5, "cm"),  #Sets overall area/size of the legend
        legend.text=element_text(size=12), #Text size
        legend.title=element_blank()) +
  geom_vline(xintercept=c(-logFC.cutoff, logFC.cutoff), linetype="longdash", colour="black", size=0.4) +
  geom_hline(yintercept=-log10(p.adj.cutoff), linetype="longdash", colour="black", size=0.4)
png(paste0(DE.path, "lung_META_tissue_resident_macrophage_volcano_1.png"), width=1200,height=800,units="px")
print(tissue_resident_macrophage.plot1)
dev.off()

#############################
### enhanced volcano plot ###
#############################
tissue_resident_macrophage.plot2 <- EnhancedVolcano(top.genes$table,
                                          lab = rownames(top.genes$table),
                                          x = 'logFC',
                                          y = 'FDR',
                                          pCutoff = 0.05,
                                          FCcutoff = 1,
                                          xlim = c(-7.5,7.5),
                                          ylim = c(0,10),
                                          title = NULL,#'DE genes in lung tissue resident macrophages',
                                          subtitle = NULL,#'healthy versus fibrotic',
                                          cutoffLineType = 'twodash',
                                          cutoffLineWidth = 0.8,
                                          pointSize = 2,
                                          labSize = 5.0,
                                          colAlpha = 0.6,
                                          legendLabels=c('NS','logFC','p-value',
                                                         'p-value & logFC'),
                                          legendPosition = 'right',
                                          legendLabSize = 10,
                                          legendIconSize = 2.0,
                                          gridlines.major = FALSE,
                                          gridlines.minor = FALSE)
png(paste0(DE.path, "lung_META_tissue_resident_macrophage_volcano_2.png"), width=1200,height=800,units="px")
print(tissue_resident_macrophage.plot2)
dev.off()

######################
### no.3: monocyte ###
######################

## prepare data ##
# subset data
pseudobulk.cell <- lung.pseudobulk[,label3==lung.pseudobulk$meta.celltype]

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
top.genes <- topTags(res, n=5291)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_lung_META_monocyte.csv"))

# save results in data frame 
monocyte.genes <- data.frame(gene = rownames(top.genes$table),
                        p.adj = top.genes$table$FDR,
                        logFC = top.genes$table$logFC)

####################
### volcano plot ###
####################
# prepare data #
p.adj.cutoff <- 0.05
logFC.cutoff <- 1.0

monocyte.genes$significance <- "NS"
monocyte.genes$significance[(abs(monocyte.genes$logFC) > logFC.cutoff)] <- "FC"
monocyte.genes$significance[(monocyte.genes$p.adj < p.adj.cutoff)] <- "FDR"
monocyte.genes$significance[(monocyte.genes$p.adj < p.adj.cutoff) & (abs(monocyte.genes$logFC) > logFC.cutoff)] <- "FC_FDR"
table(monocyte.genes$significance)
monocyte.genes$significance <- factor(monocyte.genes$significance, levels=c("NS", "FC", "FDR", "FC_FDR"))
# sort genes by ordered padj
monocyte.genes_ordered <- monocyte.genes[order(monocyte.genes$p.adj), ]
# create a column to indicate which genes to label
monocyte.genes_ordered$gene.labels <- ""
monocyte.genes_ordered$gene.labels <- monocyte.genes_ordered$gene %in% monocyte.genes_ordered$gene[1:100] 

# volcano plot
monocyte.plot1 <- ggplot(monocyte.genes_ordered, aes(x=logFC, y=-log10(p.adj))) +
  geom_point(aes(color=factor(significance)), alpha=1/2, size=2.0) +
  geom_text_repel(aes(x=logFC, y=-log10(p.adj), label = ifelse(gene.labels == T, paste0(gene),""))) +
  scale_color_manual( values=c(NS="grey30", FC="forestgreen", FDR="royalblue", FC_FDR="red2"), 
                      labels=c(NS="NS", 
                               FC="logFC", 
                               FDR="p-value (adj)", 
                               FC_FDR="logFC & p-value (adj)")) +
  scale_x_continuous(limits = c(-5,5)) +
  scale_y_continuous(limits = c(0,10)) +
  ggtitle("Lung monocyte DE genes\nhealthy vs. fibrotic") +
  xlab(bquote(~Log[2]~ "fold change")) +
  ylab(bquote(~-Log[10]~adjusted~italic(P))) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.position="top",            #Moves the legend to the top of the plot
        legend.key=element_blank(),       #removes the border
        legend.key.size=unit(0.5, "cm"),  #Sets overall area/size of the legend
        legend.text=element_text(size=12), #Text size
        legend.title=element_blank()) +
  geom_vline(xintercept=c(-logFC.cutoff, logFC.cutoff), linetype="longdash", colour="black", size=0.4) +
  geom_hline(yintercept=-log10(p.adj.cutoff), linetype="longdash", colour="black", size=0.4)
png(paste0(DE.path, "lung_META_monocyte_volcano_1.png"), width=1200,height=800,units="px")
print(monocyte.plot1)
dev.off()

#############################
### enhanced volcano plot ###
#############################
monocyte.plot2 <- EnhancedVolcano(top.genes$table,
                             lab = rownames(top.genes$table),
                             x = 'logFC',
                             y = 'FDR',
                             pCutoff = 0.05,
                             FCcutoff = 1,
                             xlim = c(-5,5),
                             ylim = c(0,10),
                             title = NULL,#'DE genes in lung monocytes',
                             subtitle = NULL,#'healthy versus fibrotic',
                             cutoffLineType = 'twodash',
                             cutoffLineWidth = 0.8,
                             pointSize = 2,
                             labSize = 5.0,
                             colAlpha = 0.6,
                             legendLabels=c('NS','logFC','p-value',
                                            'p-value & logFC'),
                             legendPosition = 'right',
                             legendLabSize = 10,
                             legendIconSize = 2.0,
                             gridlines.major = FALSE,
                             gridlines.minor = FALSE)
png(paste0(DE.path, "lung_META_monocyte_volcano_2.png"), width=1200,height=800,units="px")
print(monocyte.plot2)
dev.off()

#################
### no.4: cDC ###
#################

## prepare data ##
# subset data
pseudobulk.cell <- lung.pseudobulk[,label4==lung.pseudobulk$meta.celltype]

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
top.genes <- topTags(res, n=5204)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_lung_META_cDC.csv"))

# save results in data frame 
cDC.genes <- data.frame(gene = rownames(top.genes$table),
                        p.adj = top.genes$table$FDR,
                        logFC = top.genes$table$logFC)

####################
### volcano plot ###
####################
# prepare data #
p.adj.cutoff <- 0.05
logFC.cutoff <- 1.0

cDC.genes$significance <- "NS"
cDC.genes$significance[(abs(cDC.genes$logFC) > logFC.cutoff)] <- "FC"
cDC.genes$significance[(cDC.genes$p.adj < p.adj.cutoff)] <- "FDR"
cDC.genes$significance[(cDC.genes$p.adj < p.adj.cutoff) & (abs(cDC.genes$logFC) > logFC.cutoff)] <- "FC_FDR"
table(cDC.genes$significance)
cDC.genes$significance <- factor(cDC.genes$significance, levels=c("NS", "FC", "FDR", "FC_FDR"))
# sort genes by ordered padj
cDC.genes_ordered <- cDC.genes[order(cDC.genes$p.adj), ]
# create a column to indicate which genes to label
cDC.genes_ordered$gene.labels <- ""
cDC.genes_ordered$gene.labels <- cDC.genes_ordered$gene %in% cDC.genes_ordered$gene[1:100] 

# volcano plot
cDC.plot1 <- ggplot(cDC.genes_ordered, aes(x=logFC, y=-log10(p.adj))) +
  geom_point(aes(color=factor(significance)), alpha=1/2, size=2.0) +
  geom_text_repel(aes(x=logFC, y=-log10(p.adj), label = ifelse(gene.labels == T, paste0(gene),""))) +
  scale_color_manual( values=c(NS="grey30", FC="forestgreen", FDR="royalblue", FC_FDR="red2"), 
                      labels=c(NS="NS", 
                               FC="logFC", 
                               FDR="p-value (adj)", 
                               FC_FDR="logFC & p-value (adj)")) +
  scale_x_continuous(limits = c(-5,5)) +
  scale_y_continuous(limits = c(0,12.5)) +
  ggtitle("lung cDC DE genes\nhealthy vs. fibrotic") +
  xlab(bquote(~Log[2]~ "fold change")) +
  ylab(bquote(~-Log[10]~adjusted~italic(P))) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.position="top",            #Moves the legend to the top of the plot
        legend.key=element_blank(),       #removes the border
        legend.key.size=unit(0.5, "cm"),  #Sets overall area/size of the legend
        legend.text=element_text(size=12), #Text size
        legend.title=element_blank()) +
  geom_vline(xintercept=c(-logFC.cutoff, logFC.cutoff), linetype="longdash", colour="black", size=0.4) +
  geom_hline(yintercept=-log10(p.adj.cutoff), linetype="longdash", colour="black", size=0.4)
png(paste0(DE.path, "lung_META_cDC_volcano_1.png"), width=1200,height=800,units="px")
print(cDC.plot1)
dev.off()

#############################
### enhanced volcano plot ###
#############################
cDC.plot2 <- EnhancedVolcano(top.genes$table,
                             lab = rownames(top.genes$table),
                             x = 'logFC',
                             y = 'FDR',
                             pCutoff = 0.05,
                             FCcutoff = 1,
                             xlim = c(-5,5),
                             ylim = c(0,12.5),
                             title = NULL,#'DE genes in lung monocytes',
                             subtitle = NULL,#'healthy versus fibrotic',
                             cutoffLineType = 'twodash',
                             cutoffLineWidth = 0.8,
                             pointSize = 2,
                             labSize = 5.0,
                             colAlpha = 0.6,
                             legendLabels=c('NS','logFC','p-value',
                                            'p-value & logFC'),
                             legendPosition = 'right',
                             legendLabSize = 10,
                             legendIconSize = 2.0,
                             gridlines.major = FALSE,
                             gridlines.minor = FALSE)
png(paste0(DE.path, "lung_META_cDC_volcano_2.png"), width=1200,height=800,units="px")
print(cDC.plot2)
dev.off()
