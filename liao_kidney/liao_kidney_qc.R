library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(umap)
library(reticulate)

### path.variables ###
kidney.path <- "/home/s1987963/ds_group/Niklas/liao_kidney/"
patient.ID <- c("GSM4145204", "GSM4145205", "GSM4145206")

### load data ###
# load samples samples separately #
for(i in 1:length(patient.ID)){
  assign(paste0("kidney.",patient.ID[i],".data"), Read10X(data.dir = paste0(kidney.path, patient.ID[i] , "/")))
}

# initialize Seurat objects #
for(i in 1:length(patient.ID)){
  data <- paste0("kidney.",patient.ID[i])
  count_data <- paste0("kidney.",patient.ID[i],".data")
  eval(parse(text=paste0(data," <- CreateSeuratObject(counts = ", count_data,", project = \"stewart_kidney\", min.cells=3, min.features=200)")))
  eval(parse(text=paste0("rm(", count_data, ")")))
  
  # add metadata: study/cohort
  eval(parse(text=paste0(data,"$study <- \"liao_kidney\"")))
  eval(parse(text=paste0(data,"$cohort <- \"Guangxi\"")))
  
  # add metadata: organ
  eval(parse(text=paste0(data,"$organ <- \"kidney\"")))
  
  # add metadata: condition
  eval(parse(text=paste0(data,"$condition <- \"healthy\"")))
  
  # add metadata: patient ID
  eval(parse(text=paste0(data,"$patient.ID <- \"", patient.ID[i], "\"")))
  
  # add meta data: calculate missing metrics
  eval(parse(text=paste0(data,"[[\"percent.mt\"]] <- PercentageFeatureSet(", data, ", pattern = \"^MT-\")")))
}

# merge seurat objects: 25279 cells #
kidney <- merge(kidney.GSM4145204, c(kidney.GSM4145205, kidney.GSM4145206))

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

### QC at patient level ###
# unique genes per cell
nFeature.plot1 <- VlnPlot(kidney, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.path,"QC.raw.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(kidney, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.path,"QC.raw.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(kidney, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.path,"QC.raw.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(kidney, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.path,"QC.raw.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(kidney, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.path,"QC.raw.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(kidney, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.path,"QC.raw.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

## QC histograms ##
QC.histograms.raw <- QC.histograms(kidney)
QC.all.histograms.raw <- QC.histograms.raw[[1]]+QC.histograms.raw[[2]]+QC.histograms.raw[[3]]
png(paste0(kidney.path,"QC.raw.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.raw)
dev.off()

## QC scatter plot ##
QC.scatter.raw <- QC.scatter(kidney)
png(paste0(kidney.path,"QC.raw.scatter.png"), width=1600, height=1000,units="px")
print(QC.scatter.raw)
dev.off()

### save data: 25279 cells ###
saveRDS(kidney, paste0(kidney.path, "liao_kidney.rds"))

### doublet detection with scrublet ###
# perform on personal machine #

# load scrublet output #
kidney.scrublet <- readRDS("/home/s1987963/ds_group/Niklas/liao_kidney/liao_kidney.rds")

### how many doublets were detected by scrublet? ###
table(kidney.scrublet[[]]$scrublet_auto)

### remove  doublets ###
kidney.scrublet <- subset(kidney.scrublet, subset = scrublet_auto == FALSE)

### post-doublet removal QC plots ###
## histograms ## 
QC.histograms.scrublet <- QC.histograms(kidney.scrublet)
QC.all.histograms.scrublet <- QC.histograms.scrublet[[1]]+QC.histograms.scrublet[[2]]+QC.histograms.scrublet[[3]]
png(paste0(kidney.path,"QC.scrublet.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.scrublet)
dev.off()

## post-doublet removal QC scatter plot ##
QC.scatter.scrublet <- QC.scatter(kidney.scrublet)
png(paste0(kidney.path,"QC.scrublet.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.scrublet)
dev.off()

## post-doublet removal QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(kidney.scrublet, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.path,"QC.scrublet.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(kidney.scrublet, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.path,"QC.scrublet.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(kidney.scrublet, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.path,"QC.scrublet.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(kidney.scrublet, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.path,"QC.scrublet.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(kidney.scrublet, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.path,"QC.scrublet.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(kidney.scrublet, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.path,"QC.scrublet.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

### remove low quality cells ###
### cells with < 500 features
### cells with mitochondrial fraction < 20%
kidney.filtered <- subset(kidney.scrublet, subset = nFeature_RNA > 500 & percent.mt < 20)

### post filtering QC plots ###
## histograms ## 
QC.histograms.filtered <- QC.histograms(kidney.filtered)
QC.all.histograms.filtered <- QC.histograms.filtered[[1]]+QC.histograms.filtered[[2]]+QC.histograms.filtered[[3]]
png(paste0(kidney.path,"QC.filtered.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.filtered)
dev.off()

## QC scatter plot ##
QC.scatter.filtered <- QC.scatter(kidney.filtered)
png(paste0(kidney.path,"QC.filtered.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.filtered)
dev.off()

## post filtering QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(kidney.filtered, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.path,"QC.filtered.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(kidney.filtered, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.path,"QC.filtered.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(kidney.filtered, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.path,"QC.filtered.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(kidney.filtered, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.path,"QC.filtered.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(kidney.filtered, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.path,"QC.filtered.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(kidney.filtered, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.path,"QC.filtered.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

### save data: 16317 healthy cells ###
saveRDS(kidney.filtered, paste0(kidney.path, "liao_kidney_filtered.rds"))