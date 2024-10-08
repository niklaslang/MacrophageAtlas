library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(umap)
library(reticulate)
library(RColorBrewer)

### path.variables ###
blood.path <- "/home/s1987963/ds_group/Niklas/ramacha_blood/"

### read data ###
blood <- readRDS(paste0(blood.path, "ramacha_blood_fibrotic.rds"))

### add meta data ###
blood$organ <- "blood"
blood$condition <- "healthy"
blood$percent.mt <- blood$percent.mito
blood$percent.mito <- NULL
blood$patient.ID <- blood$sample
blood$sample <- NULL
blood$source <- NULL
blood$study <- "ramachandran_blood"
blood$cohort <- "Edinburgh"

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

### QC metrics ###
## QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(blood, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(blood.path,"QC.raw.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(blood, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(blood.path,"QC.raw.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(blood, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(blood.path,"QC.raw.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(blood, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(blood.path,"QC.raw.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(blood, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(blood.path,"QC.raw.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(blood, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(blood.path,"QC.raw.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

## QC histograms ##
QC.histograms.raw <- QC.histograms(blood)
QC.all.histograms.raw <- QC.histograms.raw[[1]]+QC.histograms.raw[[2]]+QC.histograms.raw[[3]]
png(paste0(blood.path,"QC.raw.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.raw)
dev.off()

## QC scatter plot ##
QC.scatter.raw <- QC.scatter(blood)
png(paste0(blood.path,"QC.raw.scatter.png"), width=1600, height=1000,units="px")
print(QC.scatter.raw)
dev.off()

### save data ###
saveRDS(blood, paste0(blood.path, "ramacha_blood_fibrotic.rds"))

### perform doublet detection with scrublet ###
## on personal machine ##
# load file if necessary #
blood <- readRDS(paste0(blood.path, "ramacha_blood_fibrotic.rds"))

# show scrublet results #
table(blood$scrublet_auto)

# remove 1605 doublets #
blood <- subset(blood, subset = scrublet_auto == FALSE)

### post doublet removal QC metrics ###
## QC histograms ##
QC.histograms.scrublet <- QC.histograms(blood)
QC.all.histograms.scrublet <- QC.histograms.scrublet[[1]]+QC.histograms.scrublet[[2]]+QC.histograms.scrublet[[3]]
png(paste0(blood.path,"QC.scrublet.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.scrublet)
dev.off()

## QC scatter plot ##
QC.scatter.scrublet <- QC.scatter(blood)
png(paste0(blood.path,"QC.scrublet.scatter.png"), width=1600, height=1000,units="px")
print(QC.scatter.scrublet)
dev.off()

## QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(blood, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(blood.path,"QC.scrublet.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(blood, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(blood.path,"QC.scrublet.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(blood, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(blood.path,"QC.scrublet.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(blood, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(blood.path,"QC.scrublet.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochrial fraction per cell
percent.mt.plot1 <- VlnPlot(blood, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(blood.path,"QC.scrublet.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(blood, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(blood.path,"QC.scrublet.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

### remove 2041 low quality cells ###
### cells with mitochondrial fraction < 10%
### cells with < 500 features
blood.filtered <- subset(blood, subset = nFeature_RNA > 500 & percent.mt < 10)

### post filtering QC plots ###
## histograms ## 
QC.histograms.filtered <- QC.histograms(blood.filtered)
QC.all.histograms.filtered <- QC.histograms.filtered[[1]]+QC.histograms.filtered[[2]]+QC.histograms.filtered[[3]]
png(paste0(blood.path,"QC.filtered.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.filtered)
dev.off()

## post filtering QC scatter plot ##
QC.scatter.filtered <- QC.scatter(blood.filtered)
png(paste0(blood.path,"QC.filtered.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.filtered)
dev.off()

## post filtering QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(blood.filtered, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(blood.path,"QC.filtered.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(blood.filtered, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(blood.path,"QC.filtered.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(blood.filtered, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(blood.path,"QC.filtered.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(blood.filtered, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(blood.path,"QC.filtered.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(blood.filtered, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(blood.path,"QC.filtered.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(blood.filtered, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(blood.path,"QC.filtered.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

### save healthy blood data: 45145 cells ###
saveRDS(blood.filtered, paste0(blood.path, "ramacha_blood_fibrotic_filtered.rds"))
