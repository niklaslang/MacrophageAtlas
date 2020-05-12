library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(umap)
library(reticulate)

# load lung data
lung.data.dir <- c("/home/s1987963/processed_data/raredon_lung/healthy/GSM3926545_Hum1/", "/home/s1987963/processed_data/raredon_lung/healthy/GSM3926546_Hum2/")

lung.data <- Read10X(data.dir = lung.data.dir)

# initialize the Seurat object
lung <- CreateSeuratObject(counts = lung.data, project = "lung_raredon", min.cells = 3, min.features = 200)

# compute fraction of mitochondrial gene counts
lung[["percent.mt"]] <- PercentageFeatureSet(lung, pattern = "^MT-")

# violin plots
VlnPlot(lung, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# scatter plots
plot1 <- FeatureScatter(lung, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(lung, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# histograms
# distribution of count depth
plot3 <- qplot(x =lung[["nCount_RNA"]]$nCount_RNA, fill=..count.., geom="histogram", binwidth = 250,
               xlab = "Count depth",
               ylab = "Frequency",
               main = "Count depth per cell")+scale_fill_gradient(low="blue", high="red")

# distribution of mitochondrial gene fraction
plot4 <- qplot(x =lung[["percent.mt"]]$percent.mt, fill=..count.., geom="histogram", binwidth = 0.1,
          xlab = "Fraction Mitochondrial Counts",
          ylab = "Frequency",
          main = "Mitochondrial reads fraction per cell")+scale_fill_gradient(low="yellow", high="green")

plot3 + plot4

# filter cells




