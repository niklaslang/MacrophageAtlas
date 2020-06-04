library(Seurat)
library(dplyr)
library(umap)
library(reticulate)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# path to data
reyfman_lung.path <- "/home/s1987963/ds_group/Niklas/reyfman_lung/"

### load filtered data ###
lung.filtered <- readRDS(paste0(reyfman_lung.path,"reyfman_lung_filtered.rds"))

### functions ###
# scatter plots
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
png(paste0(reyfman_lung.path,"QC.filtered.scatter.png"), width=1600, height=1000,units="px")
print(QC.scatter.filtered)
dev.off()

### how many doublets were detected by scrublet? ###
table(lung.filtered[[]]$scrublet_auto)

### remove doublets ###
lung.scrublet <- subset(lung.filtered, subset = scrublet_auto == FALSE)

### post-doublet removal QC plots ###
## histograms ## 
QC.histograms.scrublet <- QC.histograms(lung.scrublet)
QC.all.histograms.scrublet <- QC.histograms.scrublet[[1]]+QC.histograms.scrublet[[2]]+QC.histograms.scrublet[[3]]
png(paste0(reyfman_lung.path,"QC.scrublet.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.scrublet)
dev.off()

## post-doublet removal QC scatter plot ##
QC.scatter.scrublet <- QC.scatter(lung.scrublet)
png(paste0(reyfman_lung.path,"QC.scrublet.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.scrublet)
dev.off()

## post-doublet removal QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(lung.scrublet, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(reyfman_lung.path,"QC.scrublet.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(lung.scrublet, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(reyfman_lung.path,"QC.scrublet.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(lung.scrublet, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(reyfman_lung.path,"QC.scrublet.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(lung.scrublet, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(reyfman_lung.path,"QC.scrublet.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(lung.scrublet, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(reyfman_lung.path,"QC.scrublet.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(lung.scrublet, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(reyfman_lung.path,"QC.scrublet.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

### remove low quality cells ###
# cells with fewer than 1500 transcripts
# cells with mitochondrial fraction higher than 15
lung.qc1 <- subset(lung.scrublet, subset = nCount_RNA > 1500 & percent.mt < 15)

### post cut-off QC plots ###
## histograms ## 
QC.histograms.qc1 <- QC.histograms(lung.qc1)
QC.all.histograms.qc1 <- QC.histograms.qc1[[1]]+QC.histograms.qc1[[2]]+QC.histograms.qc1[[3]]
png(paste0(reyfman_lung.path,"QC.qc1.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.qc1)
dev.off()

## post cut-off QC scatter plot ##
QC.scatter.qc1 <- QC.scatter(lung.qc1)
png(paste0(reyfman_lung.path,"QC.qc1.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.qc1)
dev.off()

## post cuf-off QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(lung.qc1, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(reyfman_lung.path,"QC.qc1.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(lung.qc1, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(reyfman_lung.path,"QC.qc1.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(lung.qc1, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(reyfman_lung.path,"QC.qc1.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(lung.qc1, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(reyfman_lung.path,"QC.qc1.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(lung.qc1, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(reyfman_lung.path,"QC.qc1.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(lung.qc1, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(reyfman_lung.path,"QC.qc1.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

### split data by condition ###
lung.healthy <- subset(lung.qc1, subset = condition == "healthy")
lung.fibrotic <- subset(lung.qc1, subset = condition == "fibrotic")

### save data ###
saveRDS(lung.healthy, paste0(reyfman_lung.path, "reyfman_lung_healthy.rds"))
saveRDS(lung.fibrotic, paste0(reyfman_lung.path, "reyfman_lung_fibrotic.rds"))

