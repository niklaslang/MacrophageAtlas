library(dplyr)
library(ggplot2)
library(data.table)
library(Seurat)
library(patchwork)
library(umap)
library(reticulate)

### path.variables ###
blood.path <- "/home/s1987963/ds_group/Niklas/reyes_blood/"
metadata.path <- "/home/s1987963/ds_group/Niklas/reyes_blood/scp_meta.txt"

### load data ###
blood <- readRDS(paste0(blood.path, "reyes_blood.rds"))

### load meta data ###
blood.metadata <- fread(metadata.path, stringsAsFactors = TRUE, header = TRUE)

### clean meta data ###
blood.metadata <- blood.metadata[-1] # remove sub-header
blood.metadata$NAME <- gsub(x = blood.metadata$NAME, pattern = "-", replacement = ".") # replace "-"
blood.metadata <- blood.metadata[NAME %in% colnames(blood)] # remove all cells that are not in count matrix

# convert meta data #
blood.metadata <- data.frame(blood.metadata, row.names = blood.metadata$NAME) # set rownames equal to cell labels
blood.metadata$NAME <- NULL # remove cell label column

### add meta-data ###
blood <- AddMetaData( object = blood, metadata = blood.metadata$Cohort, col.name = "condition")
blood <- AddMetaData( object = blood, metadata = blood.metadata$Cell_Type, col.name = "cell_type")
blood <- AddMetaData( object = blood, metadata = blood.metadata$Patient, col.name = "patient.ID")
blood[["percent.mt"]] <- PercentageFeatureSet(blood, pattern = "^MT-")

### subset data ###
blood.healthy <- subset(blood, subset = condition == "Control")

### save data ###
saveRDS(blood.healthy, file = paste0(blood.path, "reyes_blood_healthy.rds"))

### QC ###

### functions ###
# scatter plots
QC.scatter <- function(data){
  scatter.plot <- ggplot(data[[]], aes( x = nFeature_RNA, y = nCount_RNA)) + 
    geom_point(aes(colour = percent.mt), size = 0.1) + 
    coord_cartesian(xlim = c(0.0 , 5000), ylim = c(0.0 , 10000)) +
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
nFeature.plot1 <- VlnPlot(blood.healthy, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(blood.path,"QC.filtered.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(blood.healthy, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(blood.path,"QC.filtered.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(blood.healthy, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(blood.path,"QC.filtered.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(blood.healthy, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(blood.path,"QC.filtered.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(blood.healthy, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(blood.path,"QC.filtered.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(blood.healthy, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(blood.path,"QC.filtered.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

## QC histograms ##
QC.histograms.filtered <- QC.histograms(blood.healthy)
QC.all.histograms.filtered <- QC.histograms.filtered[[1]]+QC.histograms.filtered[[2]]+QC.histograms.filtered[[3]]
png(paste0(blood.path,"QC.filtered.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.filtered)
dev.off()

## QC scatter plot ##
QC.scatter.filtered <- QC.scatter(blood.healthy)
png(paste0(blood.path,"QC.filtered.scatter.png"), width=1600, height=1000,units="px")
print(QC.scatter.filtered)
dev.off()

### remove doublets with scrublet ###
# see reyes_blood_scrublet.R - run on personal machine #
### how many doublets were detected by scrublet? ###
table(blood.healthy[[]]$scrublet_auto)

### remove doublets ###
blood.healthy <- subset(blood.healthy, subset = scrublet_auto == FALSE)
saveRDS(blood.healthy, file = paste0(blood.path, "reyes_blood_healthy.rds"))

### post-doublet removal QC plots ###
## histograms ## 
QC.histograms.scrublet <- QC.histograms(blood.healthy)
QC.all.histograms.scrublet <- QC.histograms.scrublet[[1]]+QC.histograms.scrublet[[2]]+QC.histograms.scrublet[[3]]
png(paste0(blood.path,"QC.scrublet.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.scrublet)
dev.off()

## post-doublet removal QC scatter plot ##
QC.scatter.scrublet <- QC.scatter(blood.healthy)
png(paste0(blood.path,"QC.scrublet.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.scrublet)
dev.off()

## post-doublet removal QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(blood.healthy, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(blood.path,"QC.scrublet.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(blood.healthy, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(blood.path,"QC.scrublet.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(blood.healthy, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(blood.path,"QC.scrublet.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(blood.healthy, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(blood.path,"QC.scrublet.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(blood.healthy, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(blood.path,"QC.scrublet.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(blood.healthy, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(blood.path,"QC.scrublet.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()
                          
### remove low quality cells ###
# cells with mitochondrial fraction higher than 15
# cells with fewer than 500 genes
blood.qc1 <- subset(blood.healthy, subset = nFeature_RNA > 500)

### post cut-off QC plots ###
## histograms ## 
QC.histograms.qc1 <- QC.histograms(blood.qc1)
QC.all.histograms.qc1 <- QC.histograms.qc1[[1]]+QC.histograms.qc1[[2]]+QC.histograms.qc1[[3]]
png(paste0(blood.path,"QC.qc1.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.qc1)
dev.off()

## post-doublet removal QC scatter plot ##
QC.scatter.qc1 <- QC.scatter(blood.qc1)
png(paste0(blood.path,"QC.qc1.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.qc1)
dev.off()

## post-doublet removal QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(blood.qc1, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(blood.path,"QC.qc1.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(blood.qc1, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(blood.path,"QC.qc1.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(blood.qc1, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(blood.path,"QC.qc1.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(blood.qc1, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(blood.path,"QC.qc1.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(blood.qc1, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(blood.path,"QC.qc1.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(blood.qc1, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(blood.path,"QC.qc1.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()
