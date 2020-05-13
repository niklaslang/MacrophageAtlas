library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(umap)
library(reticulate)

### load lung data ###
lung.data.dir <- c("/home/s1987963/processed_data/raredon_lung/healthy/GSM3926545_Hum1/", "/home/s1987963/processed_data/raredon_lung/healthy/GSM3926546_Hum2/")
lung.data <- Read10X(data.dir = lung.data.dir)

# initialize the Seurat object
lung <- CreateSeuratObject(counts = lung.data, project = "lung_raredon", min.cells = 3, min.features = 200)

### QC ###

## Calculate missing metrics ##
# compute fraction of mitochondrial gene counts
lung[["percent.mt"]] <- PercentageFeatureSet(lung, pattern = "^MT-")

## QC plots before filtering ##

# violin plots
QC.violin <- function(data){
  violin.plots <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  return(violin.plots)
}
QC.violin.before <- QC.violin(lung)
QC.violin.before

#plot1 <- VlnPlot(lung, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1

# scatter plots
QC.scatter <- function(data){
  scatter.plots <- list()
  scatter1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
  scatter.plots[[1]] <- scatter1
  scatter2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  scatter.plots[[2]] <- scatter2
  return(scatter.plots)
}
QC.scatter.before <- QC.scatter(lung)
QC.scatter.before[[1]] + QC.scatter.before[[2]]

#plot2 <- FeatureScatter(lung, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot3 <- FeatureScatter(lung, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot2 + plot3

# histograms

QC.histograms <- function(data){
  histograms <- list()
  
  # distribution of genes per cell
  hist1 <- qplot(x =data[["nFeature_RNA"]]$nFeature_RNA , fill=..count.., geom="histogram", binwidth = 100,
                 xlab = "Genes",
                 ylab = "Frequency",
                 main = "Genes per cell")+scale_fill_gradient(low="lightblue", high="darkblue")
  histograms[[1]] <- hist1
  
  # distribution of count depth
  hist2 <- qplot(x =data[["nCount_RNA"]]$nCount_RNA, fill=..count.., geom="histogram", binwidth = 1000,
                 xlab = "Count depth",
                 ylab = "Frequency",
                 main = "Count depth per cell")+scale_fill_gradient(low="orange", high="red")
  histograms[[2]] <- hist2
  
  # distribution of mitochondrial gene fraction
  hist3 <- qplot(x =data[["percent.mt"]]$percent.mt, fill=..count.., geom="histogram", binwidth = 0.1,
                    xlab = "Fraction Mitochondrial Counts",
                    ylab = "Frequency",
                    main = "Mitochondrial Reads Fraction")+scale_fill_gradient(low="lightgreen", high="darkgreen")
  histograms[[3]] <- hist3
  return(histograms)
}
QC.histograms.before <- QC.histograms(lung)
QC.histograms.before[[1]]+QC.histograms.before[[2]]+QC.histograms.before[[3]]

## distribution of count depth
#plot4 <- qplot(x =lung[["nCount_RNA"]]$nCount_RNA, fill=..count.., geom="histogram", binwidth = 250,
#               xlab = "Count depth",
#               ylab = "Frequency",
#               main = "Count depth per cell")+scale_fill_gradient(low="blue", high="red")
#
## distribution of mitochondrial gene fraction
#plot5 <- qplot(x =lung[["percent.mt"]]$percent.mt, fill=..count.., geom="histogram", binwidth = 0.1,
#          xlab = "Fraction Mitochondrial Counts",
#          ylab = "Frequency",
#          main = "Mitochondrial Reads Fraction")+scale_fill_gradient(low="yellow", high="green")
#plot4 + plot5

## QC plots with proposed cutoff ##
# plot cut-offs
QC.hist1.cutoff <- QC.histograms(lung)[[1]] + geom_vline(aes(xintercept=7500),color="black", linetype="dashed", size=.5) + geom_vline(aes(xintercept=200),color="black", linetype="dashed", size=.5)
QC.hist2.cutoff <- QC.histograms(lung)[[2]] + geom_vline(aes(xintercept=500),color="black", linetype="dashed", size=.5)
QC.hist3.cutoff <- QC.histograms(lung)[[3]] + geom_vline(aes(xintercept=7),color="black", linetype="dashed", size=.5)
QC.hist1.cutoff + QC.hist2.cutoff + QC.hist3.cutoff

## filter cells ##
lung <- subset(lung, subset = nFeature_RNA > 200)

## QC plots after filtering ##

# violin plots
QC.violin.after <- QC.violin(lung)
QC.violin.before / QC.violin.after

# scatter plots
QC.scatter.after <- QC.scatter(lung)
QC.scatter.before[[1]] + QC.scatter.after[[1]]
QC.scatter.before[[2]] + QC.scatter.after[[2]]

# histograms
QC.histograms.after <- QC.histograms(lung)
QC.histograms.before[[1]] + QC.histograms.after[[1]]
QC.histograms.before[[2]] + QC.histograms.after[[2]]
QC.histograms.before[[3]] + QC.histograms.after[[3]]

### normalisation ###
lung <- NormalizeData(lung, normalization.method = "LogNormalize", scale.factor = 10000)

# feature selection: selecting the most highly variable genes
lung <- FindVariableFeatures(lung, selection.method = "vst", nfeatures = 2000)

# 10 most highly variable genes
top10 <- head(VariableFeatures(lung), 10)

# plot variable features with and without labels
HVGs.plot <- VariableFeaturePlot(lung)
HVGs.plot <- LabelPoints(plot = HVGs.plot, points = top10, repel = TRUE)
HVGs.plot

# scale data
all.genes <- rownames(lung)
lung <- ScaleData(lung, features = all.genes)

# dimensionality reduction: PCA
lung <- RunPCA(lung, features = VariableFeatures(object = lung))

