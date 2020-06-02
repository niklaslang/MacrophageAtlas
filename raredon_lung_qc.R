library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

### load lung data ###
lung.data.dir <- "/home/s1987963/processed_data/raredon_lung/healthy/"
lung.data.patients <- c("GSM3926545_Hum1","GSM3926546_Hum2", "GSM4050097_Hum3",
                        "GSM4050099_Hum4","GSM4050101_Hum5","GSM4050103_Hum6",
                        "GSM4050105_Hum7","GSM4050107_Hum8","GSM4050109_Hum9",
                        "GSM4050111_Hum10","GSM4050112_Hum11","GSM4050113_Hum12",
                        "GSM4050114_Hum13","GSM4050115_Hum14")
lung.patient.ID <- sprintf("%02d", 1:length(lung.data.patients))

samples.path <- "/home/s1987963/ds_group/Niklas/raredon_lung/samples/"
merge.path <- "/home/s1987963/ds_group/Niklas/raredon_lung/raredon_lung.rds"
raredon_lung.path <- "/home/s1987963/ds_group/Niklas/raredon_lung/"

# load lung samples separately #
for(i in 1:length(lung.patient.ID)){
  assign(paste0("lung.patient",lung.patient.ID[i],".data"), Read10X(data.dir = paste0(lung.data.dir, lung.data.patients[i])))
}

# initialize Seurat objects #
for(i in 1:length(lung.patient.ID)){
  data <- paste0("lung.patient",lung.patient.ID[i])
  count_data <- paste0("lung.patient",lung.patient.ID[i],".data")
  pat.ID <- paste0("patient",lung.patient.ID[i])
  eval(parse(text=paste0(data," <- CreateSeuratObject(counts = ", count_data,", project = \"raredon_lung\", min.cells=3, min.features=200)")))
  eval(parse(text=paste0("rm(", count_data, ")")))
  
  # add metadata: author_dataset
  eval(parse(text=paste0(data,"$author <- \"raredon_lung\"")))
  
  # add metadata: organ
  eval(parse(text=paste0(data,"$organ <- \"lung\"")))
  
  # add metadata: condition
  eval(parse(text=paste0(data,"$condition <- \"healthy\"")))
  
  # add metadata: patient ID
  eval(parse(text=paste0(data,"$patient.ID <- \"", pat.ID, "\"")))
  
  # add meta data: calculate missing metrics
  eval(parse(text=paste0(data,"[[\"percent.mt\"]] <- PercentageFeatureSet(", data, ", pattern = \"^MT-\")")))
}

# save samples separately #
for(i in 1:length(lung.patient.ID)){
  data <- paste0("lung.patient",lung.patient.ID[i])
  eval(parse(text=paste0("saveRDS(", data, ", file = \"", samples.path, data, ".rds\")")))
}

### scrublet ###
## perform doublet removal on personal machine ##

### load scrublet output ###
lung <- readRDS(merge.path)

### QC metrics ###
## QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(lung, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(raredon_lung.path,"QC.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(lung, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(raredon_lung.path,"QC.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(lung, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(raredon_lung.path,"QC.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(lung, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(raredon_lung.path,"QC.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(lung, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(raredon_lung.path,"QC.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(lung, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(raredon_lung.path,"QC.percent.mt.2.png"), width=1500,height=500,units="px")
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
QC.scatter <- QC.scatter(lung)
png(paste0(raredon_lung.path,"QC.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter)
dev.off()

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
QC.histograms <- QC.histograms(lung)
QC.all.histograms <- QC.histograms[[1]]+QC.histograms[[2]]+QC.histograms[[3]]
png(paste0(raredon_lung.path,"QC.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms)
dev.off()

### apply cut-offs ###