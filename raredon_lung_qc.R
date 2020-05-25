library(Seurat)
library(SoupX)
library(dplyr)
library(umap)
library(reticulate)
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

# load lung samples separately
for(i in 1:length(lung.patient.ID)){
  assign(paste0("lung.patient",lung.patient.ID[i],".data"), Read10X(data.dir = paste0(lung.data.dir, lung.data.patients[i])))
}

# initialize Seurat objects
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

## merge Seurat objects
lung <- merge(lung.patient01, c(lung.patient02, lung.patient03, lung.patient04,
                                   lung.patient05, lung.patient06, lung.patient07,
                                   lung.patient08, lung.patient09, lung.patient10,
                                   lung.patient11, lung.patient12, lung.patient13,
                                   lung.patient14))

### save merged data ###
saveRDS(lung, file = "/home/s1987963/MacrophageAtlas/raredon_lung.rds")

### QC metrics ###
## QC at patient level ##
# unique genes per cell
VlnPlot(lung, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
VlnPlot(lung, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)

# unique transcripts per cell
VlnPlot(lung, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
VlnPlot(lung, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)

# mitochondrial fraction per cell
VlnPlot(lung, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
VlnPlot(lung, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)

## overall QC ##
# violin plots
QC.violin <- function(data){
  violin.plots <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = NULL)
  return(violin.plots)
}
QC.violin.before <- QC.violin(lung)
QC.violin.before

# scatter plots
QC.scatter <- function(data){
  scatter.plots <- list()
  scatter1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
  scatter.plots[[1]] <- scatter1
  scatter2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  scatter.plots[[2]] <- scatter2
  return(scatter.plots)
}
QC.scatter.before <- QC.scatter(lung)
QC.scatter.before[[1]] + QC.scatter.before[[2]]

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
                 xlab = "Unique transcripts per cell",
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
QC.histograms.before <- QC.histograms(lung)
QC.histograms.before[[1]]+QC.histograms.before[[2]]+QC.histograms.before[[3]]

### conclusion ###