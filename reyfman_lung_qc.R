library(Seurat)
library(dplyr)
library(umap)
library(reticulate)
library(ggplot2)
library(patchwork)
library(celda)

### load raw data ###
lung.raw <- readRDS("/home/s1987963/ds_group/Niklas/reyfman_lung/reyfman_lung_raw.rds")
reyfman_lung.path <- "/home/s1987963/ds_group/Niklas/reyfman_lung/"

### QC metrics ###
## QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(lung.raw, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(reyfman_lung.path,"QC.raw.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(lung.raw, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(reyfman_lung.path,"QC.raw.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(lung.raw, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(reyfman_lung.path,"QC.raw.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(lung.raw, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(reyfman_lung.path,"QC.raw.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(lung.raw, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(reyfman_lung.path,"QC.raw.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(lung.raw, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(reyfman_lung.path,"QC.raw.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

## violin plots ##
QC.violin <- function(data){
  violin.plots <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = NULL)
  return(violin.plots)
}
QC.violin.raw <- QC.violin(lung.raw)
png(paste0(reyfman_lung.path,"QC.raw.violin.png"), width=1500,height=500,units="px")
print(QC.violin.raw)
dev.off()

## scatter plots ##
QC.scatter <- function(data){
  scatter.plots <- list()
  scatter1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
  scatter.plots[[1]] <- scatter1
  scatter2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  scatter.plots[[2]] <- scatter2
  return(scatter.plots)
}
QC.scatter.raw <- QC.scatter(lung.raw)
QC.scatter.raw[[1]] + QC.scatter.raw[[2]]

## histograms ##
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
QC.histograms.raw <- QC.histograms(lung.raw)
QC.all.histograms.raw <- QC.histograms.raw[[1]]+QC.histograms.raw[[2]]+QC.histograms.raw[[3]]
png(paste0(reyfman_lung.path,"QC.raw.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.raw)
dev.off()

### load filtered data ###
lung.filtered <- readRDS("/home/s1987963/ds_group/Niklas/reyfman_lung/reyfman_lung_filtered.rds")

### QC metrics ###
## QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(lung.filtered, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(reyfman_lung.path,"QC.filtered.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(lung.filtered, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(reyfman_lung.path,"QC.filtered.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(lung.filtered, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(reyfman_lung.path,"QC.filtered.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(lung.filtered, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(reyfman_lung.path,"QC.filtered.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(lung.filtered, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(reyfman_lung.path,"QC.filtered.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(lung.filtered, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(reyfman_lung.path,"QC.filtered.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

## QC histograms ##
QC.histograms.filtered <- QC.histograms(lung.filtered)
QC.all.histograms.filtered <- QC.histograms.filtered[[1]]+QC.histograms.filtered[[2]]+QC.histograms.filtered[[3]]
png(paste0(reyfman_lung.path,"QC.filtered.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.filtered)
dev.off()