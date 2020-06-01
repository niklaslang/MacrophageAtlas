library(Seurat)
library(dplyr)
library(umap)
library(reticulate)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# path to data
reyfman_lung.path <- "/home/s1987963/ds_group/Niklas/reyfman_lung/"

### load raw data ###
lung.raw <- readRDS("/home/s1987963/ds_group/Niklas/reyfman_lung/reyfman_lung_raw.rds")

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

## scatter plot ##
QC.scatter <- function(data){
  scatter.plot <- ggplot(data[[]], aes( x = nFeature_RNA, y = nCount_RNA)) + 
    geom_point(aes(colour = percent.mt), size = 0.1) + 
    coord_cartesian(xlim = c(0.0 , 10000), ylim = c(0.0 , 100000)) +
    labs(title = "Overall QC", x  ="Count depth", y = "Unique Genes") + 
    theme(
      plot.title = element_text(color = "black", size = 20 , face = "bold"),
      axis.title.x = element_text(color = "black", size = 20, face = "bold"),
      axis.title.y = element_text(color = "black", size = 20, face = "bold"),
      legend.title = element_text(color = "black", size = 16, face = "bold", angle = 90)
    ) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")),
                           guide = guide_colourbar("Mitochondrial fraction", title.position = "right", title.vjust = 0.5, title.hjust = 0.5, barwidth = 1.0, barheight = 60))
    
  return(scatter.plot)
}
QC.scatter.raw <- QC.scatter(lung.raw)
png(paste0(reyfman_lung.path,"QC.raw.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.raw)
dev.off()

#guides(fill = guide_colourbar("Mitochondrial fraction", title.position = "right", barwidth = 1.0, barheight = 10))

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
                 xlab = "Count depth per cell",
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

## QC scatter plot ##
QC.scatter.filtered <- QC.scatter(lung.filtered)
png(paste0(reyfman_lung.path,"QC.filtered.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.filtered)
dev.off()

### apply QC cut-off no. 1 ###
lung.cutoff1 <- subset(lung.filtered, subset = nCount_RNA > 1200 & percent.mt < 20)

## post cut-off QC histograms ##
QC.histograms.cutoff1 <- QC.histograms(lung.cutoff1)
QC.all.histograms.cutoff1 <- QC.histograms.cutoff1[[1]]+QC.histograms.cutoff1[[2]]+QC.histograms.cutoff1[[3]]
png(paste0(reyfman_lung.path,"QC.cutoff1.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.cutoff1)
dev.off()

## post cut-off QC scatter plot ##
QC.scatter.cutoff1 <- QC.scatter(lung.cutoff1)
png(paste0(reyfman_lung.path,"QC.cutoff1.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.cutoff1)
dev.off()

## post cuf-off QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(lung.cutoff1, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(reyfman_lung.path,"QC.cutoff1.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(lung.cutoff1, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(reyfman_lung.path,"QC.cutoff1.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(lung.cutoff1, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(reyfman_lung.path,"QC.cutoff1.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(lung.cutoff1, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(reyfman_lung.path,"QC.cutoff1.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(lung.cutoff1, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(reyfman_lung.path,"QC.cutoff1.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(lung.cutoff1, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(reyfman_lung.path,"QC.cutoff1.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()




