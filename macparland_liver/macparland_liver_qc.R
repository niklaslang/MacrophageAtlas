library(dplyr)
library(ggplot2)
library(data.table)
library(Seurat)
library(patchwork)
library(umap)
library(reticulate)

### path.variables ###
liver.path <- "/home/s1987963/ds_group/Niklas/macparland_liver/"

### load liver data ###
liver <- readRDS(paste0(liver.path, "macparland_liver.rds"))

### add meta data ###
liver$organ <- "liver"
liver$condition <- "healthy"
liver$patient.ID <- liver$orig.ident
liver$percent.mt <- liver$percent.mito
liver$percent.mito <- NULL

### functions ###
# scatter plots
QC.scatter <- function(data){
  scatter.plot <- ggplot(data[[]], aes( x = nCount_RNA, y = nFeature_RNA)) + 
    geom_point(aes(colour = percent.mt), size = 1) + 
    coord_cartesian(xlim = c(0.0 , 10000), ylim = c(0.0 , 10000)) +
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
  hist2 <- qplot(x =data[["nCount_RNA"]]$nCount_RNA, fill=..count.., geom="histogram", binwidth = 200,
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
nFeature.plot1 <- VlnPlot(liver, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.path,"QC.raw.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(liver, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.path,"QC.raw.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(liver, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.path,"QC.raw.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(liver, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.path,"QC.raw.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochrial fraction per cell
percent.mt.plot1 <- VlnPlot(liver, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.path,"QC.raw.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(liver, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.path,"QC.raw.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

## QC histograms ##
QC.histograms.raw <- QC.histograms(liver)
QC.all.histograms.raw <- QC.histograms.raw[[1]]+QC.histograms.raw[[2]]+QC.histograms.raw[[3]]
png(paste0(liver.path,"QC.raw.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.raw)
dev.off()

## QC scatter plot ##
QC.scatter.raw <- QC.scatter(liver)
png(paste0(liver.path,"QC.raw.scatter.png"), width=1600, height=1000,units="px")
print(QC.scatter.raw)
dev.off()

### remove doublets ###
# run scrublet script on personal machine #

### check scrublet results ###
table(liver$scrublet_auto == FALSE)

### remove 118 doublet ###
liver <- subset(liver, subset = scrublet_auto == FALSE)

### post-doublet removal QC plots ###
## histograms ## 
QC.histograms.scrublet <- QC.histograms(liver)
QC.all.histograms.scrublet <- QC.histograms.scrublet[[1]]+QC.histograms.scrublet[[2]]+QC.histograms.scrublet[[3]]
png(paste0(liver.path,"QC.scrublet.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.scrublet)
dev.off()

## post-doublet removal QC scatter plot ##
QC.scatter.scrublet <- QC.scatter(liver)
png(paste0(liver.path,"QC.scrublet.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.scrublet)
dev.off()

## post-doublet removal QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(liver, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.path,"QC.scrublet.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(liver, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.path,"QC.scrublet.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(liver, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.path,"QC.scrublet.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(liver, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.path,"QC.scrublet.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(liver, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.path,"QC.scrublet.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(liver, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.path,"QC.scrublet.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

### remove low quality cells ###
### cells with < 500 features
liver.filtered <- subset(liver, subset = nFeature_RNA > 500)

### post filtering QC plots ###
## histograms ## 
QC.histograms.filtered <- QC.histograms(liver.filtered)
QC.all.histograms.filtered <- QC.histograms.filtered[[1]]+QC.histograms.filtered[[2]]+QC.histograms.filtered[[3]]
png(paste0(liver.path,"QC.filtered.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.filtered)
dev.off()

## QC scatter plot ##
QC.scatter.filtered <- QC.scatter(liver.filtered)
png(paste0(liver.path,"QC.filtered.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.filtered)
dev.off()

## post filtering QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(liver.filtered, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.path,"QC.filtered.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(liver.filtered, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.path,"QC.filtered.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(liver.filtered, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.path,"QC.filtered.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(liver.filtered, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.path,"QC.filtered.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(liver.filtered, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(liver.path,"QC.filtered.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(liver.filtered, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(liver.path,"QC.filtered.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

### save data ###
saveRDS(liver.filtered, paste0(liver.path, "macparland_liver_filtered.rds"))
