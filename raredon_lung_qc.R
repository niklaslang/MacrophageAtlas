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
