library(Seurat)
library(dplyr)
library(umap)
library(reticulate)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# path to data
reyfman_lung.path <- "/home/s1987963/ds_group/Niklas/reyfman_lung/"
lung.healthy.path <- "/home/s1987963/ds_group/Niklas/reyfman_lung/healthy/"
lung.fibrotic.path <- "/home/s1987963/ds_group/Niklas/reyfman_lung/fibrotic/"

### load filtered data ###
lung.all <- readRDS(paste0(reyfman_lung.path,"reyfman_lung_filtered.rds")) # 79175cells

### split data ###
lung.healthy <- subset(lung.all, subset = condition == "healthy") # 43627 cells
lung.fibrotic <- subset(lung.all, subset = condition == "fibrotic") # 35548 cells

### save data ###
# healthy cells
saveRDS(lung.healthy, paste0(lung.healthy.path, "habermann_lung_healthy.rds"))
# fibrotic cells 
saveRDS(lung.fibrotic, paste0(lung.fibrotic.path, "habermann_lung_fibrotic.rds"))

### perform QC ###
### functions ###
# scatter plots
QC.scatter <- function(data){
  scatter.plot <- ggplot(data[[]], aes( x = nCount_RNA, y = nFeature_RNA)) + 
    geom_point(aes(colour = percent.mt), size = 1) + 
    coord_cartesian(xlim = c(0.0 , 50000), ylim = c(0.0 , 10000)) +
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
  hist2 <- qplot(x =data[["nCount_RNA"]]$nCount_RNA, fill=..count.., geom="histogram", binwidth = 500,
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

########################
#####  ALL DATA    #####
########################

### QC metrics ###
## QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(lung.all, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(reyfman_lung.path,"QC.filtered.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(lung.all, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(reyfman_lung.path,"QC.filtered.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(lung.all, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(reyfman_lung.path,"QC.filtered.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(lung.all, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(reyfman_lung.path,"QC.filtered.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(lung.all, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(reyfman_lung.path,"QC.filtered.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(lung.all, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(reyfman_lung.path,"QC.filtered.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

## QC histograms ##
QC.histograms.filtered <- QC.histograms(lung.all)
QC.all.histograms.filtered <- QC.histograms.filtered[[1]]+QC.histograms.filtered[[2]]+QC.histograms.filtered[[3]]
png(paste0(reyfman_lung.path,"QC.filtered.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.filtered)
dev.off()

## QC scatter plot ##
QC.scatter.filtered <- QC.scatter(lung.all)
png(paste0(reyfman_lung.path,"QC.filtered.scatter.png"), width=1600, height=1000,units="px")
print(QC.scatter.filtered)
dev.off()

########################
####  HEALTHY DATA  ####
########################

### QC metrics ###
## QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(lung.healthy, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.healthy.path,"QC.filtered.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(lung.healthy, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.healthy.path,"QC.filtered.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(lung.healthy, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.healthy.path,"QC.filtered.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(lung.healthy, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.healthy.path,"QC.filtered.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(lung.healthy, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.healthy.path,"QC.filtered.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(lung.healthy, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.healthy.path,"QC.filtered.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

## QC histograms ##
QC.histograms.filtered <- QC.histograms(lung.healthy)
QC.all.histograms.filtered <- QC.histograms.filtered[[1]]+QC.histograms.filtered[[2]]+QC.histograms.filtered[[3]]
png(paste0(lung.healthy.path,"QC.filtered.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.filtered)
dev.off()

## QC scatter plot ##
QC.scatter.filtered <- QC.scatter(lung.healthy)
png(paste0(lung.healthy.path,"QC.filtered.scatter.png"), width=1600, height=1000,units="px")
print(QC.scatter.filtered)
dev.off()

### scrublet: detected 1149 doublets ###
table(lung.healthy[[]]$scrublet_auto)

### remove doublets ###
lung.healthy <- subset(lung.healthy, subset = scrublet_auto == FALSE)

### post-doublet removal QC plots ###
## histograms ## 
QC.histograms.scrublet <- QC.histograms(lung.healthy)
QC.all.histograms.scrublet <- QC.histograms.scrublet[[1]]+QC.histograms.scrublet[[2]]+QC.histograms.scrublet[[3]]
png(paste0(lung.healthy.path,"QC.scrublet.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.scrublet)
dev.off()

## post-doublet removal QC scatter plot ##
QC.scatter.scrublet <- QC.scatter(lung.healthy)
png(paste0(lung.healthy.path,"QC.scrublet.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.scrublet)
dev.off()

## post-doublet removal QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(lung.healthy, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.healthy.path,"QC.scrublet.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(lung.healthy, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.healthy.path,"QC.scrublet.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(lung.healthy, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.healthy.path,"QC.scrublet.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(lung.healthy, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.healthy.path,"QC.scrublet.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(lung.healthy, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.healthy.path,"QC.scrublet.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(lung.healthy, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.healthy.path,"QC.scrublet.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

### remove low quality cells ###
# cells with mitochondrial fraction higher than 15
# cells with fewer than 500 and more than 6,000 genes
lung.healthy <- subset(lung.healthy, subset = nFeature_RNA > 500 & percent.mt < 15)

### post cut-off QC plots ###
## histograms ## 
QC.histograms.healthy <- QC.histograms(lung.healthy)
QC.all.histograms.healthy <- QC.histograms.healthy[[1]]+QC.histograms.healthy[[2]]+QC.histograms.healthy[[3]]
png(paste0(lung.healthy.path,"QC.healthy.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.healthy)
dev.off()

## post cut-off QC scatter plot ##
QC.scatter.healthy <- QC.scatter(lung.healthy)
png(paste0(lung.healthy.path,"QC.healthy.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.healthy)
dev.off()

## post cuf-off QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(lung.healthy, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.healthy.path,"QC.healthy.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(lung.healthy, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.healthy.path,"QC.healthy.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(lung.healthy, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.healthy.path,"QC.healthy.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(lung.healthy, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.healthy.path,"QC.healthy.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(lung.healthy, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.healthy.path,"QC.healthy.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(lung.healthy, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.healthy.path,"QC.healthy.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

### save healthy data ###
saveRDS(lung.healthy, paste0(lung.healthy.path, "reyfman_lung_healthy_filtered.rds")) # 42264 cells

########################
####  FIBROTIC DATA  ###
########################

### QC at patient level ###
# unique genes per cell
nFeature.plot1 <- VlnPlot(lung.fibrotic, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.fibrotic.path,"QC.raw.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(lung.fibrotic, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.fibrotic.path,"QC.raw.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(lung.fibrotic, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.fibrotic.path,"QC.raw.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(lung.fibrotic, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.fibrotic.path,"QC.raw.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(lung.fibrotic, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.fibrotic.path,"QC.raw.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(lung.fibrotic, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.fibrotic.path,"QC.raw.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

## QC histograms ##
QC.histograms.raw <- QC.histograms(lung.fibrotic)
QC.all.histograms.raw <- QC.histograms.raw[[1]]+QC.histograms.raw[[2]]+QC.histograms.raw[[3]]
png(paste0(lung.fibrotic.path,"QC.raw.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.raw)
dev.off()

## QC scatter plot ##
QC.scatter.raw <- QC.scatter(lung.fibrotic)
png(paste0(lung.fibrotic.path,"QC.raw.scatter.png"), width=1600, height=1000,units="px")
print(QC.scatter.raw)
dev.off()

### doublet detection with scrublet ###
# perform on personal machine #

### scrublet results ###
table(lung.fibrotic[[]]$scrublet_auto)

### remove 859 doublets ###
lung.fibrotic <- subset(lung.fibrotic, subset = scrublet_auto == FALSE)

### post-doublet removal QC plots ###
## histograms ## 
QC.histograms.scrublet <- QC.histograms(lung.fibrotic)
QC.all.histograms.scrublet <- QC.histograms.scrublet[[1]]+QC.histograms.scrublet[[2]]+QC.histograms.scrublet[[3]]
png(paste0(lung.fibrotic.path,"QC.scrublet.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.scrublet)
dev.off()

## post-doublet removal QC scatter plot ##
QC.scatter.scrublet <- QC.scatter(lung.fibrotic)
png(paste0(lung.fibrotic.path,"QC.scrublet.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.scrublet)
dev.off()

## post-doublet removal QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(lung.fibrotic, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.fibrotic.path,"QC.scrublet.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(lung.fibrotic, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.fibrotic.path,"QC.scrublet.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(lung.fibrotic, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.fibrotic.path,"QC.scrublet.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(lung.fibrotic, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.fibrotic.path,"QC.scrublet.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(lung.fibrotic, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.fibrotic.path,"QC.scrublet.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(lung.fibrotic, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.fibrotic.path,"QC.scrublet.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

### remove low quality cells ###
### cells with < 500 features
### cells with mitochondrial fraction < 15%
lung.fibrotic.filtered <- subset(lung.fibrotic, subset = nFeature_RNA > 500 & percent.mt < 15)

### post filtering QC plots ###
## histograms ## 
QC.histograms.filtered <- QC.histograms(lung.fibrotic.filtered)
QC.all.histograms.filtered <- QC.histograms.filtered[[1]]+QC.histograms.filtered[[2]]+QC.histograms.filtered[[3]]
png(paste0(lung.fibrotic.path,"QC.filtered.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.filtered)
dev.off()

## QC scatter plot ##
QC.scatter.filtered <- QC.scatter(lung.fibrotic.filtered)
png(paste0(lung.fibrotic.path,"QC.filtered.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.filtered)
dev.off()

## post filtering QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(lung.fibrotic.filtered, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.fibrotic.path,"QC.filtered.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(lung.fibrotic.filtered, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.fibrotic.path,"QC.filtered.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(lung.fibrotic.filtered, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.fibrotic.path,"QC.filtered.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(lung.fibrotic.filtered, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.fibrotic.path,"QC.filtered.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(lung.fibrotic.filtered, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.fibrotic.path,"QC.filtered.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(lung.fibrotic.filtered, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.fibrotic.path,"QC.filtered.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

### save data: 31888 healthy cells ###
saveRDS(lung.fibrotic.filtered, paste0(lung.fibrotic.path, "reyfman_lung_fibrotic_filtered.rds"))