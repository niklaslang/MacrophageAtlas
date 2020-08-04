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
kidney.path <- "/home/s1987963/ds_group/Niklas/combined_organs/combined_kidney_annotated.rds"
DE.path <- "/home/s1987963/ds_group/Niklas/all_kidney/DE/"

### read files ###
kidney <- readRDS(kidney.path)

### reorder condition factor ###
kidney$condition <- factor(kidney$condition)
kidney$condition <- factor(kidney$condition, levels(kidney$condition)[c(2,1)])

### convert to sce ###
kidney.sce <- as.SingleCellExperiment(kidney)

### overview ###
table(kidney.sce$celltype_integrated, kidney.sce$condition)

# creating pseudo-bulk samples
kidney.pseudobulk <- aggregateAcrossCells(kidney.sce, id=colData(kidney.sce)[,c("celltype_integrated", "patient.ID")])

# picking one cell type
label1 <- "Combined Kidney Monocyte"
label2 <- "Combined Kidney Macrophage"

######################
### no.1: Monocyte ###
######################

## prepare data ##
# subset data
pseudobulk.cell <- kidney.pseudobulk[,label1==kidney.pseudobulk$celltype_integrated]

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
summary(decideTests(res)) # 965 down and 1057 up
top.genes <- topTags(res, n=2572)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_kidney_monocyte.csv"))

# save results in data frame 
monocyte.genes <- data.frame(gene = rownames(top.genes$table),
                                p.adj = top.genes$table$FDR,
                                logFC = top.genes$table$logFC)

####################
### volcano plot ###
####################
# prepare data #
p.adj.cutoff <- 0.05
logFC.cutoff <- 2.0

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
monocyte.genes_ordered$gene.labels <- monocyte.genes_ordered$gene %in% monocyte.genes_ordered$gene[1:60] 

# volcano plot
monocyte.plot <- ggplot(monocyte.genes_ordered, aes(x=logFC, y=-log10(p.adj))) +
  geom_point(aes(color=factor(significance)), alpha=1/2, size=2.0) +
  geom_text_repel(aes(x=logFC, y=-log10(p.adj), label = ifelse(gene.labels == T, paste0(gene),""))) +
  scale_color_manual( values=c(NS="grey30", FC="forestgreen", FDR="royalblue", FC_FDR="red2"), 
                      labels=c(NS="NS", 
                               FC="logFC", 
                               FDR="p-value (adj)", 
                               FC_FDR="logFC & p-value (adj)")) +
  scale_x_continuous(limits = c(-10,10)) +
  scale_y_continuous(limits = c(0,20)) +
  ggtitle("Kidney monocyte DE genes\nhealthy vs. fibrotic") +
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
png(paste0(DE.path, "kidney_monocyte_volcano_1.png"), width=1200,height=800,units="px")
print(monocyte.plot)
dev.off()

#############################
### enhanced volcano plot ###
#############################
monocyte.plot2 <- EnhancedVolcano(top.genes$table,
                                    lab = rownames(top.genes$table),
                                    x = 'logFC',
                                    y = 'FDR',
                                    pCutoff = 0.05,
                                    FCcutoff = 2,
                                    xlim = c(-10,10),
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
png(paste0(DE.path, "kidney_monocyte_volcano_2.png"), width=1200,height=800,units="px")
print(monocyte.plot2)
dev.off()

#########################
### no.2: Macrophages ###
#########################

## prepare data ##
# subset data
pseudobulk.cell <- kidney.pseudobulk[,label2==kidney.pseudobulk$celltype_integrated]

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
summary(decideTests(res)) # 859 down and 788 up
top.genes <- topTags(res, n=2031)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_kidney_macrophage.csv"))

# save results in data frame 
macrophages.genes <- data.frame(gene = rownames(top.genes$table),
                                p.adj = top.genes$table$FDR,
                                logFC = top.genes$table$logFC)

####################
### volcano plot ###
####################
# prepare data #
p.adj.cutoff <- 0.05
logFC.cutoff <- 2.0

macrophages.genes$significance <- "NS"
macrophages.genes$significance[(abs(macrophages.genes$logFC) > logFC.cutoff)] <- "FC"
macrophages.genes$significance[(macrophages.genes$p.adj < p.adj.cutoff)] <- "FDR"
macrophages.genes$significance[(macrophages.genes$p.adj < p.adj.cutoff) & (abs(macrophages.genes$logFC) > logFC.cutoff)] <- "FC_FDR"
table(macrophages.genes$significance)
macrophages.genes$significance <- factor(macrophages.genes$significance, levels=c("NS", "FC", "FDR", "FC_FDR"))
# sort genes by ordered padj
macrophages.genes_ordered <- macrophages.genes[order(macrophages.genes$p.adj), ]
# create a column to indicate which genes to label
macrophages.genes_ordered$gene.labels <- ""
macrophages.genes_ordered$gene.labels <- macrophages.genes_ordered$gene %in% macrophages.genes_ordered$gene[1:60] 

# volcano plot
macrophage.plot <- ggplot(macrophages.genes_ordered, aes(x=logFC, y=-log10(p.adj))) +
  geom_point(aes(color=factor(significance)), alpha=1/2, size=2.0) +
  geom_text_repel(aes(x=logFC, y=-log10(p.adj), label = ifelse(gene.labels == T, paste0(gene),""))) +
  scale_color_manual( values=c(NS="grey30", FC="forestgreen", FDR="royalblue", FC_FDR="red2"), 
                      labels=c(NS="NS", 
                               FC="logFC", 
                               FDR="p-value (adj)", 
                               FC_FDR="logFC & p-value (adj)")) +
  scale_x_continuous(limits = c(-10,10)) +
  scale_y_continuous(limits = c(0,15)) +
  ggtitle("Kidney macrophages DE genes\nhealthy vs. fibrotic") +
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
png(paste0(DE.path, "kidney_macrophages_volcano_1.png"), width=1200,height=800,units="px")
print(macrophage.plot)
dev.off()

#############################
### enhanced volcano plot ###
#############################
macrophage.plot2 <- EnhancedVolcano(top.genes$table,
                lab = rownames(top.genes$table),
                x = 'logFC',
                y = 'FDR',
                pCutoff = 0.05,
                FCcutoff = 2,
                xlim = c(-10,10),
                title = NULL,#'DE genes in kidney macrophages',
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
png(paste0(DE.path, "kidney_macrophages_volcano_2.png"), width=1200,height=800,units="px")
print(macrophage.plot2)
dev.off()