library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(umap)
library(reticulate)

### path.variables ###
lung.path <- "/home/s1987963/ds_group/Niklas/habermann_lung/"
lung.healthy.path <- "/home/s1987963/ds_group/Niklas/habermann_lung/healthy/"
lung.fibrotic.path <- "/home/s1987963/ds_group/Niklas/habermann_lung/fibrotic/"
metadata.path <- "/home/s1987963/ds_group/Niklas/habermann_lung/GSE135893_IPF_metadata.csv"

### load lung data ###
lung.counts <- Read10X(lung.path, gene.column=1)

### create seurat object: 212566 cells ###
lung <- CreateSeuratObject(counts = lung.counts, project = "habermann_lung", min.cells = 3, min.features = 200)

### add meta data ###
# calculate fraction of mitochondrial counts #
lung$percent.mt <- PercentageFeatureSet(lung, pattern = "^MT-")
lung$organ <- "lung"

### load meta data ###
lung.metadata.dt <- fread(metadata.path, stringsAsFactors = TRUE, header = TRUE)
lung.metadata.dt <- data.table(lung.metadata.dt)

# clean filtered meta data #
lung.metadata.dt$V1 <- NULL
lung.metadata.dt$nCount_RNA <- NULL
lung.metadata.dt$nFeature_RNA <- NULL
lung.metadata.dt$nCount_SCT <- NULL
lung.metadata.dt$nFeature_SCT <- NULL
lung.metadata.dt$percent.mt <- NULL
lung.metadata.dt$seurat_clusters <- NULL
lung.metadata.dt$celltype <- NULL
lung.metadata.dt$population <- NULL
lung.metadata.dt <- unique(lung.metadata.dt)
lung.metadata.dt[, condition := ifelse(Status == "Disease", "fibrotic", "healthy")]

# reconstruct complete meta data #
metadata.dt <- data.table(barcode=colnames(lung))
metadata.dt[, sample := strsplit(barcode, "_")[[1]][1], by=barcode]
lung.metadata.dt <- merge(metadata.dt, lung.metadata.dt, by.x = "sample", by.y = "orig.ident")

# clean complete meta data #
lung.metadata.dt$Sample_Source <- NULL

# convert complete meta data #
lung.metadata <- data.frame(lung.metadata.dt, row.names = lung.metadata.dt$barcode) # set rownames equal to cell labels
lung.metadata$barcode <- NULL # remove barcode column

# add meta data to seurat object #
lung <- AddMetaData( object = lung, metadata = lung.metadata$condition, col.name = "condition")
lung <- AddMetaData( object = lung, metadata = lung.metadata$Sample_Name, col.name = "patient.ID")

### split data by condition: healthy and fibrotic ###
lung.healthy <- subset(lung, subset = condition == "healthy")
lung.fibrotic <- subset(lung, subset = condition == "fibrotic")
rm(lung)

### save data ###
# 56222 healthy cells #
saveRDS(lung.healthy, paste0(lung.path, "habermann_lung_healthy.rds"))
# 156344 fibrotic cells #4
saveRDS(lung.fibrotic, paste0(lung.path, "habermann_lung_fibrotic.rds"))

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

### QC metrics ###
## QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(lung.healthy, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.healthy.path,"QC.raw.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(lung.healthy, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.healthy.path,"QC.raw.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(lung.healthy, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.healthy.path,"QC.raw.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(lung.healthy, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.healthy.path,"QC.raw.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(lung.healthy, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.healthy.path,"QC.raw.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(lung.healthy, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.healthy.path,"QC.raw.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

## QC histograms ##
QC.histograms.raw <- QC.histograms(lung.healthy)
QC.all.histograms.raw <- QC.histograms.raw[[1]]+QC.histograms.raw[[2]]+QC.histograms.raw[[3]]
png(paste0(lung.healthy.path,"QC.raw.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.raw)
dev.off()

## QC scatter plot ##
QC.scatter.raw <- QC.scatter(lung.healthy)
png(paste0(lung.healthy.path,"QC.raw.scatter.png"), width=1600, height=1000,units="px")
print(QC.scatter.raw)
dev.off()

### doublet detection with scrublet ###
# perform on personal machine #

# load scrublet output #
lung.healthy <- readRDS(paste0(lung.path, "habermann_lung_healthy.rds"))

### how many doublets were detected by scrublet? ###
table(lung.healthy[[]]$scrublet_auto)

### remove doublets ###
lung.scrublet <- subset(lung.healthy, subset = scrublet_auto == FALSE)

### post-doublet removal QC plots ###
## histograms ## 
QC.histograms.scrublet <- QC.histograms(lung.scrublet)
QC.all.histograms.scrublet <- QC.histograms.scrublet[[1]]+QC.histograms.scrublet[[2]]+QC.histograms.scrublet[[3]]
png(paste0(lung.healthy.path,"QC.scrublet.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.scrublet)
dev.off()

## post-doublet removal QC scatter plot ##
QC.scatter.scrublet <- QC.scatter(lung.scrublet)
png(paste0(lung.healthy.path,"QC.scrublet.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.scrublet)
dev.off()

## post-doublet removal QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(lung.scrublet, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.healthy.path,"QC.scrublet.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(lung.scrublet, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.healthy.path,"QC.scrublet.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(lung.scrublet, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.healthy.path,"QC.scrublet.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(lung.scrublet, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.healthy.path,"QC.scrublet.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(lung.scrublet, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.healthy.path,"QC.scrublet.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(lung.scrublet, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.healthy.path,"QC.scrublet.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

### remove low quality cells ###
### cells with < 500 features
### cells with mitochondrial fraction < 30%
lung.filtered <- subset(lung.scrublet, subset = nFeature_RNA > 500 & percent.mt < 25)

### post filtering QC plots ###
## histograms ## 
QC.histograms.filtered <- QC.histograms(lung.filtered)
QC.all.histograms.filtered <- QC.histograms.filtered[[1]]+QC.histograms.filtered[[2]]+QC.histograms.filtered[[3]]
png(paste0(lung.healthy.path,"QC.filtered.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.filtered)
dev.off()

## QC scatter plot ##
QC.scatter.filtered <- QC.scatter(lung.filtered)
png(paste0(lung.healthy.path,"QC.filtered.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.filtered)
dev.off()

## post filtering QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(lung.filtered, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.healthy.path,"QC.filtered.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(lung.filtered, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.healthy.path,"QC.filtered.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(lung.filtered, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.healthy.path,"QC.filtered.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(lung.filtered, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.healthy.path,"QC.filtered.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(lung.filtered, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(lung.healthy.path,"QC.filtered.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(lung.filtered, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(lung.healthy.path,"QC.filtered.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

### save data: 42922 healthy cells ###
saveRDS(lung.filtered, paste0(lung.healthy.path, "reyfman_lung_healthy.rds"))

















