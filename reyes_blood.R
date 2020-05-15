library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(umap)
library(reticulate)

### load blood data ###
blood.dir <- "/home/s1987963/processed_data/reyes_blood/healthy/spc_gex_SeuratObj.rds"
blood <- readRDS(blood.dir)

### QC ###

## Calculate missing metrics ##
# compute fraction of mitochondrial gene counts
blood[["percent.mt"]] <- PercentageFeatureSet(blood, pattern = "^MT-")

## QC plots before filtering ##

# violin plots
QC.violin <- function(data){
  violin.plots <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  return(violin.plots)
}
QC.violin.before <- QC.violin(blood)
QC.violin.before

# scatter plots
QC.scatter <- function(data){
  scatter.plots <- list()
  scatter1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
  scatter.plots[[1]] <- scatter1
  scatter2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  scatter.plots[[2]] <- scatter2
  return(scatter.plots)
}
QC.scatter.before <- QC.scatter(blood)
QC.scatter.before[[1]] + QC.scatter.before[[2]]

# histograms

QC.histograms <- function(data){
  histograms <- list()
  
  # distribution of genes per cell
  hist1 <- qplot(x =data[["nFeature_RNA"]]$nFeature_RNA , fill=..count.., geom="histogram", binwidth = 100,
                 xlab = "Unique genes per cell",
                 ylab = "Frequency",
                 main = "Gene Count Distribution")+scale_fill_gradient(low="lightblue", high="darkblue")
  histograms[[1]] <- hist1
  
  # distribution of count depth
  hist2 <- qplot(x =data[["nCount_RNA"]]$nCount_RNA, fill=..count.., geom="histogram", binwidth = 100,
                 xlab = "Unique transcripts per cell",
                 ylab = "Frequency",
                 main = "Transcript Count Distribution")+scale_fill_gradient(low="orange", high="red")
  histograms[[2]] <- hist2
  
  # distribution of mitochondrial gene fraction
  hist3 <- qplot(x =data[["percent.mt"]]$percent.mt, fill=..count.., geom="histogram", binwidth = 0.1,
                 xlab = "Mitochondrial fraction per cell",
                 ylab = "Frequency",
                 main = "Mitochondrial Reads Fraction")+scale_fill_gradient(low="lightgreen", high="darkgreen")
  histograms[[3]] <- hist3
  return(histograms)
}
QC.histograms.before <- QC.histograms(blood)
QC.histograms.before[[1]]+QC.histograms.before[[2]]+QC.histograms.before[[3]]

### QC conclusion: QC already performed ###

### normalisation ###
blood <- NormalizeData(blood, normalization.method = "LogNormalize", scale.factor = 10000)

### feature selection ###
# selecting the most highly variable genes
blood <- FindVariableFeatures(blood, selection.method = "vst", nfeatures = 2000)

# 10 most highly variable genes
top10 <- head(VariableFeatures(blood), 10)

# plot variable features with and without labels
HVGs.plot <- VariableFeaturePlot(blood)
HVGs.plot <- LabelPoints(plot = HVGs.plot, points = top10, repel = TRUE)
HVGs.plot

### scale data ###
blood <- ScaleData(blood, vars.to.regress = "percent.mt")

### save data ###
saveRDS(blood, file = "/home/s1987963/MacrophageAtlas/reyes_blood.rds")



