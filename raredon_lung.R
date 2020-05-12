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
          main = "Mitochondrial Reads Fraction")+scale_fill_gradient(low="yellow", high="green")
plot3 + plot4

# filter cells
# plot cut-offs
plot5 <- plot3 + geom_vline(aes(xintercept=22000),color="black", linetype="dashed", size=.5) + geom_vline(aes(xintercept=200),color="black", linetype="dashed", size=.5)
plot6 <- plot4 + geom_vline(aes(xintercept=5),color="black", linetype="dashed", size=.5)
plot5 + plot6

# apply cut-offs
lung <- subset(lung, subset = nFeature_RNA > 200 & nFeature_RNA < 22000 & percent.mt < 5)

# normalize count matrix
lung <- NormalizeData(lung, normalization.method = "LogNormalize", scale.factor = 10000)

# feature selection: selecting the most highly variable genes
lung <- FindVariableFeatures(lung, selection.method = "vst", nfeatures = 2000)

# 10 most highly variable genes
top10 <- head(VariableFeatures(lung), 10)

# plot variable features with and without labels
plot7 <- VariableFeaturePlot(lung)
plot8 <- LabelPoints(plot = plot7, points = top10, repel = TRUE)
plot7
plot8

# scale data
all.genes <- rownames(lung)
lung <- ScaleData(lung, features = all.genes)

# dimensionality reduction: PCA
lung <- RunPCA(lung, features = VariableFeatures(object = lung))
