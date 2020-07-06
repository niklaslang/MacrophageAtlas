library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(umap)
library(reticulate)

### path.variables ###
kidney.path <- "/home/s1987963/ds_group/Niklas/wilson_kidney/"
kidney.healthy.path <- "/home/s1987963/ds_group/Niklas/wilson_kidney/healthy/"
kidney.diabetes.path <- "/home/s1987963/ds_group/Niklas/wilson_kidney/diabetes/"

# patient IDs
patient.ID <- c("GSM3823939_control", "GSM3823940_control", "GSM3823941_control",
                "GSM3823942_diabetes", "GSM3823943_diabetes", "GSM3823944_diabetes")

### load data ###
# load samples samples separately #
for(i in 1:length(patient.ID)){

  data <- paste0("kidney.", patient.ID[i])
  assign(paste0("kidney.",patient.ID[i]), readRDS(paste0(kidney.path, "kidney.", patient.ID[i], ".rds")))
  
  # add metadata: study/cohort
  eval(parse(text=paste0(data,"$study <- \"wilson_kidney\"")))
  eval(parse(text=paste0(data,"$cohort <- \"Washington/Boston\"")))
  
  # add metadata: organ
  eval(parse(text=paste0(data,"$organ <- \"kidney\"")))
  
  # add metadata: condition
  if(substr(patient.ID[i],12,14) == "con"){
    eval(parse(text=paste0(data,"$condition <- \"healthy\"")))
  } else if (substr(patient.ID[i],12,14) == "dia"){
    eval(parse(text=paste0(data,"$condition <- \"diabetes\"")))
  }
  
  # add metadata: patient ID
  eval(parse(text=paste0(data,"$patient.ID <- \"", patient.ID[i], "\"")))
}

### merge data ###
kidney.healthy <- merge(kidney.GSM3823939_control, c(kidney.GSM3823940_control, kidney.GSM3823941_control))
kidney.diabetes <- merge(kidney.GSM3823942_diabetes, c(kidney.GSM3823943_diabetes, kidney.GSM3823944_diabetes))

### add meta data ###
kidney.healthy$percent.mt <- PercentageFeatureSet(kidney.healthy, pattern = "^MT-")
kidney.diabetes$percent.mt <- PercentageFeatureSet(kidney.diabetes, pattern = "^MT-")

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
####  HEALTHY DATA  ####
########################

### QC at patient level ###
# unique genes per cell
nFeature.plot1 <- VlnPlot(kidney.healthy, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.healthy.path,"QC.raw.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(kidney.healthy, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.healthy.path,"QC.raw.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(kidney.healthy, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.healthy.path,"QC.raw.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(kidney.healthy, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.healthy.path,"QC.raw.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(kidney.healthy, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.healthy.path,"QC.raw.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(kidney.healthy, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.healthy.path,"QC.raw.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

## QC histograms ##
QC.histograms.raw <- QC.histograms(kidney.healthy)
QC.all.histograms.raw <- QC.histograms.raw[[1]]+QC.histograms.raw[[2]]+QC.histograms.raw[[3]]
png(paste0(kidney.healthy.path,"QC.raw.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.raw)
dev.off()

## QC scatter plot ##
QC.scatter.raw <- QC.scatter(kidney.healthy)
png(paste0(kidney.healthy.path,"QC.raw.scatter.png"), width=1600, height=1000,units="px")
print(QC.scatter.raw)
dev.off()

### save data ###
saveRDS(kidney.healthy, paste0(kidney.healthy.path, "wilson_kidney_healthy.rds"))

### doublet detection with scrublet ###
# perform on personal machine #

# load scrublet output #
kidney.healthy <- readRDS(paste0(kidney.healthy.path, "wilson_kidney_healthy.rds"))

### how many doublets were detected by scrublet? ###
table(kidney.healthy[[]]$scrublet_auto)

### remove 574 doublets ###
kidney.healthy.scrublet <- subset(kidney.healthy, subset = scrublet_auto == FALSE)

### post-doublet removal QC plots ###
## histograms ## 
QC.histograms.scrublet <- QC.histograms(kidney.healthy.scrublet)
QC.all.histograms.scrublet <- QC.histograms.scrublet[[1]]+QC.histograms.scrublet[[2]]+QC.histograms.scrublet[[3]]
png(paste0(kidney.healthy.path,"QC.scrublet.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.scrublet)
dev.off()

## post-doublet removal QC scatter plot ##
QC.scatter.scrublet <- QC.scatter(kidney.healthy.scrublet)
png(paste0(kidney.healthy.path,"QC.scrublet.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.scrublet)
dev.off()

## post-doublet removal QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(kidney.healthy.scrublet, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.healthy.path,"QC.scrublet.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(kidney.healthy.scrublet, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.healthy.path,"QC.scrublet.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(kidney.healthy.scrublet, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.healthy.path,"QC.scrublet.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(kidney.healthy.scrublet, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.healthy.path,"QC.scrublet.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(kidney.healthy.scrublet, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.healthy.path,"QC.scrublet.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(kidney.healthy.scrublet, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.healthy.path,"QC.scrublet.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

### remove low quality cells ###
### cells with < 500 features
### cells with mitochondrial fraction < 25%
kidney.healthy.filtered <- subset(kidney.healthy.scrublet, subset = nFeature_RNA > 500 & percent.mt < 25)

### post filtering QC plots ###
## histograms ## 
QC.histograms.filtered <- QC.histograms(kidney.healthy.filtered)
QC.all.histograms.filtered <- QC.histograms.filtered[[1]]+QC.histograms.filtered[[2]]+QC.histograms.filtered[[3]]
png(paste0(kidney.healthy.path,"QC.filtered.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.filtered)
dev.off()

## QC scatter plot ##
QC.scatter.filtered <- QC.scatter(kidney.healthy.filtered)
png(paste0(kidney.healthy.path,"QC.filtered.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.filtered)
dev.off()

## post filtering QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(kidney.healthy.filtered, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.healthy.path,"QC.filtered.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(kidney.healthy.filtered, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.healthy.path,"QC.filtered.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(kidney.healthy.filtered, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.healthy.path,"QC.filtered.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(kidney.healthy.filtered, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.healthy.path,"QC.filtered.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(kidney.healthy.filtered, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.healthy.path,"QC.filtered.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(kidney.healthy.filtered, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.healthy.path,"QC.filtered.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

### save data: 13207 healthy cells ###
saveRDS(kidney.healthy.filtered, paste0(kidney.healthy.path, "wilson_kidney_healthy_filtered.rds"))

#########################
####  FIBROTIC DATA  ####
#########################

### QC at patient level ###
# unique genes per cell
nFeature.plot1 <- VlnPlot(kidney.diabetes, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.diabetes.path,"QC.raw.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(kidney.diabetes, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.diabetes.path,"QC.raw.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(kidney.diabetes, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.diabetes.path,"QC.raw.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(kidney.diabetes, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.diabetes.path,"QC.raw.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(kidney.diabetes, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.diabetes.path,"QC.raw.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(kidney.diabetes, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.diabetes.path,"QC.raw.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

## QC histograms ##
QC.histograms.raw <- QC.histograms(kidney.diabetes)
QC.all.histograms.raw <- QC.histograms.raw[[1]]+QC.histograms.raw[[2]]+QC.histograms.raw[[3]]
png(paste0(kidney.diabetes.path,"QC.raw.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.raw)
dev.off()

## QC scatter plot ##
QC.scatter.raw <- QC.scatter(kidney.diabetes)
png(paste0(kidney.diabetes.path,"QC.raw.scatter.png"), width=1600, height=1000,units="px")
print(QC.scatter.raw)
dev.off()

### save data ###
saveRDS(kidney.diabetes, paste0(kidney.diabetes.path, "wilson_kidney_diabetes.rds"))

### doublet detection with scrublet ###
# perform on personal machine #

# load scrublet output #
kidney.diabetes <- readRDS(paste0(kidney.diabetes.path, "wilson_kidney_diabetes.rds"))

### how many doublets were detected by scrublet? ###
table(kidney.diabetes[[]]$scrublet_auto)

### remove 127 doublets ###
kidney.diabetes.scrublet <- subset(kidney.diabetes, subset = scrublet_auto == FALSE)

### post-doublet removal QC plots ###
## histograms ## 
QC.histograms.scrublet <- QC.histograms(kidney.diabetes.scrublet)
QC.all.histograms.scrublet <- QC.histograms.scrublet[[1]]+QC.histograms.scrublet[[2]]+QC.histograms.scrublet[[3]]
png(paste0(kidney.diabetes.path,"QC.scrublet.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.scrublet)
dev.off()

## post-doublet removal QC scatter plot ##
QC.scatter.scrublet <- QC.scatter(kidney.diabetes.scrublet)
png(paste0(kidney.diabetes.path,"QC.scrublet.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.scrublet)
dev.off()

## post-doublet removal QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(kidney.diabetes.scrublet, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.diabetes.path,"QC.scrublet.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(kidney.diabetes.scrublet, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.diabetes.path,"QC.scrublet.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(kidney.diabetes.scrublet, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.diabetes.path,"QC.scrublet.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(kidney.diabetes.scrublet, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.diabetes.path,"QC.scrublet.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(kidney.diabetes.scrublet, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.diabetes.path,"QC.scrublet.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(kidney.diabetes.scrublet, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.diabetes.path,"QC.scrublet.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

### remove low quality cells ###
### cells with < 500 features
### cells with mitochondrial fraction < 30%
kidney.diabetes.filtered <- subset(kidney.diabetes.scrublet, subset = nFeature_RNA > 500 & percent.mt < 25)

### post filtering QC plots ###
## histograms ## 
QC.histograms.filtered <- QC.histograms(kidney.diabetes.filtered)
QC.all.histograms.filtered <- QC.histograms.filtered[[1]]+QC.histograms.filtered[[2]]+QC.histograms.filtered[[3]]
png(paste0(kidney.diabetes.path,"QC.filtered.histogram.png"), width=1500,height=500,units="px")
print(QC.all.histograms.filtered)
dev.off()

## QC scatter plot ##
QC.scatter.filtered <- QC.scatter(kidney.diabetes.filtered)
png(paste0(kidney.diabetes.path,"QC.filtered.scatter.png"), width=1600,height=1000,units="px")
print(QC.scatter.filtered)
dev.off()

## post filtering QC at patient level ##
# unique genes per cell
nFeature.plot1 <- VlnPlot(kidney.diabetes.filtered, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.diabetes.path,"QC.filtered.nFeature_RNA.1.png"), width=1500,height=500,units="px")
print(nFeature.plot1)
dev.off()

nFeature.plot2 <- VlnPlot(kidney.diabetes.filtered, features = c("nFeature_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.diabetes.path,"QC.filtered.nFeature_RNA.2.png"), width=1500,height=500,units="px")
print(nFeature.plot2)
dev.off()

# unique transcripts per cell
nCount.plot1 <- VlnPlot(kidney.diabetes.filtered, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.diabetes.path,"QC.filtered.nCount_RNA.1.png"), width=1500,height=500,units="px")
print(nCount.plot1)
dev.off()

nCount.plot2 <- VlnPlot(kidney.diabetes.filtered, features = c("nCount_RNA"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.diabetes.path,"QC.filtered.nCount_RNA.2.png"), width=1500,height=500,units="px")
print(nCount.plot2)
dev.off()

# mitochondrial fraction per cell
percent.mt.plot1 <- VlnPlot(kidney.diabetes.filtered, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0.5)
png(paste0(kidney.diabetes.path,"QC.filtered.percent.mt.1.png"), width=1500,height=500,units="px")
print(percent.mt.plot1)
dev.off()

percent.mt.plot2 <- VlnPlot(kidney.diabetes.filtered, features = c("percent.mt"), group.by = "patient.ID", pt.size = 0)
png(paste0(kidney.diabetes.path,"QC.filtered.percent.mt.2.png"), width=1500,height=500,units="px")
print(percent.mt.plot2)
dev.off()

### save data: 9945 diabetic cells ###
saveRDS(kidney.diabetes.filtered, paste0(kidney.diabetes.path, "wilson_kidney_diabetes_filtered.rds"))


