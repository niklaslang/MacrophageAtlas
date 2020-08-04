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
blood.path <- "/home/s1987963/ds_group/Niklas/combined_organs/combined_blood_annotated.rds"
DE.path <- "/home/s1987963/ds_group/Niklas/all_blood/DE/"

### read files ###
blood <- readRDS(blood.path)

### create meta cell types ###
cell.data <- data.table(barcode = colnames(blood),
                        celltype = blood$celltype_integrated)
cell.data[, meta.celltype := ifelse(celltype %in% c("Combined Blood CD14+ Monocyte 1", "Combined Blood CD14+ Monocyte 2", "Combined Blood CD14+ Monocyte 3"),
                                    "Combined Blood CD14+ Monocyte", paste0(celltype))]
cell.data[, meta.celltype := ifelse(celltype %in% c("Combined Blood cDC1", "Combined Blood cDC2"),
                                    "Combined Blood cDC", paste0(meta.celltype))] 
cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$barcode <- NULL
cell.data$celltype <- NULL
blood <- AddMetaData(blood, cell.data, col.name = "meta.celltype")

### convert to sce ###
blood.sce <- as.SingleCellExperiment(blood)

### overview ###
table(blood.sce$meta.celltype, blood.sce$condition)

# creating pseudo-bulk samples
blood.pseudobulk <- aggregateAcrossCells(blood.sce, id=colData(blood.sce)[,c("meta.celltype", "patient.ID")])

# label meta cell types
label1 <- "Combined Blood CD14+ Monocyte"
label2 <- "Combined Blood CD16+ Monocyte"
label3 <- "Combined Blood cDC"

##############################
### no.1: CD14+ Monocytes  ###
##############################

## prepare data ##
# subset data
pseudobulk.cell <- blood.pseudobulk[,label1==blood.pseudobulk$meta.celltype]

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
summary(decideTests(res)) #  1507 down and 805 up
top.genes <- topTags(res, n=10664)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_blood_META_cd14+_monocytes.csv"))

# save results in data frame 
cd14pos_monocyte.genes <- data.frame(gene = rownames(top.genes$table),
                                     p.adj = top.genes$table$FDR,
                                     logFC = top.genes$table$logFC)

####################
### volcano plot ###
####################
# prepare data #
p.adj.cutoff <- 0.05
logFC.cutoff <- 1.0

cd14pos_monocyte.genes$significance <- "NS"
cd14pos_monocyte.genes$significance[(abs(cd14pos_monocyte.genes$logFC) > logFC.cutoff)] <- "FC"
cd14pos_monocyte.genes$significance[(cd14pos_monocyte.genes$p.adj < p.adj.cutoff)] <- "FDR"
cd14pos_monocyte.genes$significance[(cd14pos_monocyte.genes$p.adj < p.adj.cutoff) & (abs(cd14pos_monocyte.genes$logFC) > logFC.cutoff)] <- "FC_FDR"
table(cd14pos_monocyte.genes$significance)
cd14pos_monocyte.genes$significance <- factor(cd14pos_monocyte.genes$significance, levels=c("NS", "FC", "FDR", "FC_FDR"))
# sort genes by ordered padj
cd14pos_monocyte.genes_ordered <- cd14pos_monocyte.genes[order(cd14pos_monocyte.genes$p.adj), ]
# create a column to indicate which genes to label
cd14pos_monocyte.genes_ordered$gene.labels <- ""
cd14pos_monocyte.genes_ordered$gene.labels <- cd14pos_monocyte.genes_ordered$gene %in% cd14pos_monocyte.genes_ordered$gene[1:60] 

# volcano plot
cd14pos_monocyte.plot1 <- ggplot(cd14pos_monocyte.genes_ordered, aes(x=logFC, y=-log10(p.adj))) +
  geom_point(aes(color=factor(significance)), alpha=1/2, size=2.0) +
  geom_text_repel(aes(x=logFC, y=-log10(p.adj), label = ifelse(gene.labels == T, paste0(gene),""))) +
  scale_color_manual( values=c(NS="grey30", FC="forestgreen", FDR="royalblue", FC_FDR="red2"), 
                      labels=c(NS="NS", 
                               FC="logFC", 
                               FDR="p-value (adj)", 
                               FC_FDR="logFC & p-value (adj)")) +
  scale_x_continuous(limits = c(-20,20)) +
  scale_y_continuous(limits = c(0,20)) +
  ggtitle("Blood CD14+ DE genes\nhealthy vs. fibrotic") +
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
png(paste0(DE.path, "blood_META_cd14pos_monocyte_volcano_1.png"), width=1200,height=800,units="px")
print(cd14pos_monocyte.plot1)
dev.off()

#############################
### enhanced volcano plot ###
#############################
cd14pos_monocyte.plot2 <- EnhancedVolcano(top.genes$table,
                                          lab = rownames(top.genes$table),
                                          x = 'logFC',
                                          y = 'FDR',
                                          pCutoff = 0.05,
                                          FCcutoff = 1,
                                          xlim = c(-20,20),
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
png(paste0(DE.path, "blood_META_cd14pos_monocyte_volcano_2.png"), width=1200,height=800,units="px")
print(cd14pos_monocyte.plot2)
dev.off()

##############################
### no.2: CD16+ Monocytes  ###
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
summary(decideTests(res)) # 605 down and 509 up
top.genes <- topTags(res, n=6749)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_blood_META_cd16+_monocyte.csv"))

# save results in data frame 
cd16pos_monocyte.genes <- data.frame(gene = rownames(top.genes$table),
                                     p.adj = top.genes$table$FDR,
                                     logFC = top.genes$table$logFC)

####################
### volcano plot ###
####################
# prepare data #
p.adj.cutoff <- 0.05
logFC.cutoff <- 1.0

cd16pos_monocyte.genes$significance <- "NS"
cd16pos_monocyte.genes$significance[(abs(cd16pos_monocyte.genes$logFC) > logFC.cutoff)] <- "FC"
cd16pos_monocyte.genes$significance[(cd16pos_monocyte.genes$p.adj < p.adj.cutoff)] <- "FDR"
cd16pos_monocyte.genes$significance[(cd16pos_monocyte.genes$p.adj < p.adj.cutoff) & (abs(cd16pos_monocyte.genes$logFC) > logFC.cutoff)] <- "FC_FDR"
table(cd16pos_monocyte.genes$significance)
cd16pos_monocyte.genes$significance <- factor(cd16pos_monocyte.genes$significance, levels=c("NS", "FC", "FDR", "FC_FDR"))
# sort genes by ordered padj
cd16pos_monocyte.genes_ordered <- cd16pos_monocyte.genes[order(cd16pos_monocyte.genes$p.adj), ]
# create a column to indicate which genes to label
cd16pos_monocyte.genes_ordered$gene.labels <- ""
cd16pos_monocyte.genes_ordered$gene.labels <- cd16pos_monocyte.genes_ordered$gene %in% cd16pos_monocyte.genes_ordered$gene[1:60] 

# volcano plot
cd16pos_monocyte.plot1 <- ggplot(cd16pos_monocyte.genes_ordered, aes(x=logFC, y=-log10(p.adj))) +
  geom_point(aes(color=factor(significance)), alpha=1/2, size=2.0) +
  geom_text_repel(aes(x=logFC, y=-log10(p.adj), label = ifelse(gene.labels == T, paste0(gene),""))) +
  scale_color_manual( values=c(NS="grey30", FC="forestgreen", FDR="royalblue", FC_FDR="red2"), 
                      labels=c(NS="NS", 
                               FC="logFC", 
                               FDR="p-value (adj)", 
                               FC_FDR="logFC & p-value (adj)")) +
  scale_x_continuous(limits = c(-10,10)) +
  scale_y_continuous(limits = c(0,25)) +
  ggtitle("Blood CD16+ monocyte DE genes\nhealthy vs. fibrotic") +
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
png(paste0(DE.path, "blood_META_cd16pos_monocyte_volcano_1.png"), width=1200,height=800,units="px")
print(cd16pos_monocyte.plot1)
dev.off()

#############################
### enhanced volcano plot ###
#############################
cd16pos_monocyte.plot2 <- EnhancedVolcano(top.genes$table,
                                          lab = rownames(top.genes$table),
                                          x = 'logFC',
                                          y = 'FDR',
                                          pCutoff = 0.05,
                                          FCcutoff = 1,
                                          xlim = c(-10,10),
                                          title = NULL,#'DE genes in blood monocytes',
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
png(paste0(DE.path, "blood_META_cd16pos_monocyte_volcano_2.png"), width=1200,height=800,units="px")
print(cd16pos_monocyte.plot2)
dev.off()

##################
### no.3: cDCs ###
##################

## prepare data ##
# subset data
pseudobulk.cell <- blood.pseudobulk[,label3==blood.pseudobulk$meta.celltype]

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
summary(decideTests(res)) # 600 down and 583 up
top.genes <- topTags(res, n=7838)
write.csv(top.genes$table, file = paste0(DE.path, "DE_combined_blood_META_cDC.csv"))

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
cDC.genes_ordered$gene.labels <- cDC.genes_ordered$gene %in% cDC.genes_ordered$gene[1:60] 

# volcano plot
cDC.plot1 <- ggplot(cDC.genes_ordered, aes(x=logFC, y=-log10(p.adj))) +
  geom_point(aes(color=factor(significance)), alpha=1/2, size=2.0) +
  geom_text_repel(aes(x=logFC, y=-log10(p.adj), label = ifelse(gene.labels == T, paste0(gene),""))) +
  scale_color_manual( values=c(NS="grey30", FC="forestgreen", FDR="royalblue", FC_FDR="red2"), 
                      labels=c(NS="NS", 
                               FC="logFC", 
                               FDR="p-value (adj)", 
                               FC_FDR="logFC & p-value (adj)")) +
  scale_x_continuous(limits = c(-15,15)) +
  scale_y_continuous(limits = c(0,45)) +
  ggtitle("Blood cDC DE genes\nhealthy vs. fibrotic") +
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
png(paste0(DE.path, "blood_META_cDC_volcano_1.png"), width=1200,height=800,units="px")
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
                             xlim = c(-15,15),
                             title = NULL,#'DE genes in blood monocytes',
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
png(paste0(DE.path, "blood_META_cDC_volcano_2.png"), width=1200,height=800,units="px")
print(cDC.plot2)
dev.off()
