library(dplyr)
library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(umap)
library(reticulate)
library(RColorBrewer)

### path.variables ###
liver.path <- "/home/s1987963/ds_group/Niklas/ramacha_liver/"
liver.healthy.path <- "/home/s1987963/ds_group/Niklas/ramacha_liver/healthy/"
liver.fibrotic.path <- "/home/s1987963/ds_group/Niklas/ramacha_liver/fibrotic/"

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

####################
### HEALTHY DATA ###
####################

### load data ###
liver.healthy <- readRDS(paste0(liver.path, "Healthy_Livers_MergeSeuratObj.rds"))

### remove 4962 blood cells ###
liver.healthy <- subset(liver.healthy, subset = sample %in% c("Healthy_1", "Healthy_2", "Healthy_3", "Healthy_4", "Healthy_5", "Healthy_6", "Healthy_7"))

### add meta data ###
liver.healthy$organ <- "liver"
liver.healthy$condition <- "healthy"
liver.healthy$percent.mt <- liver.healthy$percent.mito
liver.healthy$percent.mito <- NULL
liver.healthy$patient.ID <- liver.healthy$sample
liver.healthy$sample <- NULL
liver.healthy$source <- NULL

### QC metrics ###
## QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(liver.healthy, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.healthy.path,"QC.raw.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(liver.healthy, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.healthy.path,"QC.raw.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(liver.healthy, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.healthy.path,"QC.raw.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(liver.healthy, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.healthy.path,"QC.raw.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(liver.healthy, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.healthy.path,"QC.raw.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(liver.healthy, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.healthy.path,"QC.raw.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

## QC histograms ##
QC.histograms.raw <- QC.histograms(liver.healthy)
QC.all.histograms.raw <- QC.histograms.raw[[1]]+QC.histograms.raw[[2]]+QC.histograms.raw[[3]]
png(paste0(liver.healthy.path,"QC.raw.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.raw)
dev.off()

## QC scatter plot ##
QC.scatter.raw <- QC.scatter(liver.healthy)
png(paste0(liver.healthy.path,"QC.raw.scatter.png"), width=1600, height=1000,units="px")
print(QC.scatter.raw)
dev.off()

### save healthy liver data: 54362 cells ###
saveRDS(liver.healthy, paste0(liver.healthy.path, "ramachandran_liver_healthy.rds"))

### perform doublet detection with scrublet ###
## on personal machine ##
# load file if necessary #
liver.healthy <- readRDS(paste0(liver.healthy.path, "ramachandran_liver_healthy.rds"))

# show scrublet results #
table(liver.healthy$scrublet_auto)

# remove 967 doublets #
liver.healthy <- subset(liver.healthy, subset = scrublet_auto == FALSE)

### post doublet removal QC metrics ###
## QC histograms ##
QC.histograms.scrublet <- QC.histograms(liver.healthy)
QC.all.histograms.scrublet <- QC.histograms.scrublet[[1]]+QC.histograms.scrublet[[2]]+QC.histograms.scrublet[[3]]
png(paste0(liver.healthy.path,"QC.scrublet.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.scrublet)
dev.off()

## QC scatter plot ##
QC.scatter.scrublet <- QC.scatter(liver.healthy)
png(paste0(liver.healthy.path,"QC.scrublet.scatter.png"), width=1600, height=1000,units="px")
print(QC.scatter.scrublet)
dev.off()

## QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(liver.healthy, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.healthy.path,"QC.scrublet.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(liver.healthy, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.healthy.path,"QC.scrublet.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(liver.healthy, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.healthy.path,"QC.scrublet.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(liver.healthy, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.healthy.path,"QC.scrublet.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochrial fraction per cell
percent.mt.plot1 <- VlnPlot(liver.healthy, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.healthy.path,"QC.scrublet.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(liver.healthy, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.healthy.path,"QC.scrublet.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

### remove 8250 low quality cells ###
### cells with mitochondrial fraction < 30%
### cells with < 500 features
liver.healthy.filtered <- subset(liver.healthy, subset = nFeature_RNA > 500 & percent.mt < 30)

### post filtering QC plots ###
## histograms ## 
QC.histograms.filtered <- QC.histograms(liver.healthy.filtered)
QC.all.histograms.filtered <- QC.histograms.filtered[[1]]+QC.histograms.filtered[[2]]+QC.histograms.filtered[[3]]
png(paste0(liver.healthy.path,"QC.filtered.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.filtered)
dev.off()

## post filtering QC scatter plot ##
QC.scatter.filtered <- QC.scatter(liver.healthy.filtered)
png(paste0(liver.healthy.path,"QC.filtered.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.filtered)
dev.off()

## post filtering QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(liver.healthy.filtered, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.healthy.path,"QC.filtered.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(liver.healthy.filtered, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.healthy.path,"QC.filtered.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(liver.healthy.filtered, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.healthy.path,"QC.filtered.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(liver.healthy.filtered, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.healthy.path,"QC.filtered.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(liver.healthy.filtered, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.healthy.path,"QC.filtered.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(liver.healthy.filtered, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.healthy.path,"QC.filtered.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

### save healthy liver data: 45145 cells ###
saveRDS(liver.filtered, paste0(liver.healthy.path, "ramacha_liver_healthy_filtered.rds"))

#####################
### FIBROTIC DATA ###
#####################

### load data ###
liver.fibrotic <- readRDS(paste0(liver.path, "cirrhotic_human_preprocessed.rds"))

### subset 32481 blood cells ###
blood.fibrotic <- subset(liver.fibrotic, subset = source == "Blood")
saveRDS(blood.fibrotic, "/home/s1987963/ds_group/Niklas/ramacha_blood/ramacha_blood_fibrotic.rds")
### subset 42745 fibrotic cells ###
liver.fibrotic <- subset(liver.fibrotic, subset = source == "Tissue")

### add meta data ###
liver.fibrotic$organ <- "liver"
liver.fibrotic$condition <- "fibrotic"
liver.fibrotic$percent.mt <- liver.fibrotic$percent.mito
liver.fibrotic$percent.mito <- NULL
liver.fibrotic$patient.ID <- liver.fibrotic$sample
liver.fibrotic$sample <- NULL
liver.fibrotic$source <- NULL
liver.fibrotic$study <- "ramachandran_liver"
liver.fibrotic$cohort <- "Edinburgh"

### QC metrics ###
## QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(liver.fibrotic, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.fibrotic.path,"QC.raw.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(liver.fibrotic, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.fibrotic.path,"QC.raw.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(liver.fibrotic, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.fibrotic.path,"QC.raw.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(liver.fibrotic, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.fibrotic.path,"QC.raw.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(liver.fibrotic, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.fibrotic.path,"QC.raw.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(liver.fibrotic, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.fibrotic.path,"QC.raw.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

## QC histograms ##
QC.histograms.raw <- QC.histograms(liver.fibrotic)
QC.all.histograms.raw <- QC.histograms.raw[[1]]+QC.histograms.raw[[2]]+QC.histograms.raw[[3]]
png(paste0(liver.fibrotic.path,"QC.raw.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.raw)
dev.off()

## QC scatter plot ##
QC.scatter.raw <- QC.scatter(liver.fibrotic)
png(paste0(liver.fibrotic.path,"QC.raw.scatter.png"), width=1600, height=1000,units="px")
print(QC.scatter.raw)
dev.off()

### save fibrotic liver data: 42745 cells ###
saveRDS(liver.fibrotic, paste0(liver.fibrotic.path, "ramachandran_liver_fibrotic.rds"))

### perform doublet detection with scrublet ###
## on personal machine ##
# load file if necessary #
liver.fibrotic <- readRDS(paste0(liver.fibrotic.path, "ramachandran_liver_fibrotic.rds"))

# show scrublet results #
table(liver.fibrotic$scrublet_auto)

# remove 682 doublets #
liver.fibrotic <- subset(liver.fibrotic, subset = scrublet_auto == FALSE)

### post doublet removal QC metrics ###
## QC histograms ##
QC.histograms.scrublet <- QC.histograms(liver.fibrotic)
QC.all.histograms.scrublet <- QC.histograms.scrublet[[1]]+QC.histograms.scrublet[[2]]+QC.histograms.scrublet[[3]]
png(paste0(liver.fibrotic.path,"QC.scrublet.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.scrublet)
dev.off()

## QC scatter plot ##
QC.scatter.scrublet <- QC.scatter(liver.fibrotic)
png(paste0(liver.fibrotic.path,"QC.scrublet.scatter.png"), width=1600, height=1000,units="px")
print(QC.scatter.scrublet)
dev.off()

## QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(liver.fibrotic, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.fibrotic.path,"QC.scrublet.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(liver.fibrotic, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.fibrotic.path,"QC.scrublet.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(liver.fibrotic, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.fibrotic.path,"QC.scrublet.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(liver.fibrotic, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.fibrotic.path,"QC.scrublet.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochrial fraction per cell
percent.mt.plot1 <- VlnPlot(liver.fibrotic, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.fibrotic.path,"QC.scrublet.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(liver.fibrotic, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.fibrotic.path,"QC.scrublet.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

### remove 5674 low quality cells ###
### cells with mitochondrial fraction < 30%
### cells with < 500 features
liver.fibrotic.filtered <- subset(liver.fibrotic, subset = nFeature_RNA > 500 & percent.mt < 30)

### post filtering QC plots ###
## histograms ## 
QC.histograms.filtered <- QC.histograms(liver.fibrotic.filtered)
QC.all.histograms.filtered <- QC.histograms.filtered[[1]]+QC.histograms.filtered[[2]]+QC.histograms.filtered[[3]]
png(paste0(liver.fibrotic.path,"QC.filtered.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.filtered)
dev.off()

## post filtering QC scatter plot ##
QC.scatter.filtered <- QC.scatter(liver.fibrotic.filtered)
png(paste0(liver.fibrotic.path,"QC.filtered.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.filtered)
dev.off()

## post filtering QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(liver.fibrotic.filtered, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.fibrotic.path,"QC.filtered.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(liver.fibrotic.filtered, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.fibrotic.path,"QC.filtered.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(liver.fibrotic.filtered, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.fibrotic.path,"QC.filtered.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(liver.fibrotic.filtered, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.fibrotic.path,"QC.filtered.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(liver.fibrotic.filtered, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.fibrotic.path,"QC.filtered.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(liver.fibrotic.filtered, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.fibrotic.path,"QC.filtered.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

### save healthy liver data: 37071 cells ###
saveRDS(liver.fibrotic.filtered, paste0(liver.fibrotic.path, "ramacha_liver_fibrotic_filtered.rds"))

