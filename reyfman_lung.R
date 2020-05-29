library(Seurat)
library(dplyr)
library(umap)
library(reticulate)
library(ggplot2)
library(patchwork)
library(celda)

### read raw data ###
raw.dir <- "reyfman_raw/"
patient.ID <- paste0(c(rep("healthy_",8),rep("fibrotic_",8)), rep(sprintf("%02d", 1:8),2))

# load raw lung samples separately #
for(i in 1:length(patient.ID)){
  #print(paste0(raw.dir, "healthy_", patient.ID[i], "_raw.h5"))
  #print(paste0(raw.dir, "fibrotic_", patient.ID[i], "_raw.h5"))
  assign(paste0("lung.", patient.ID[i],".raw.data"), Read10X_h5(paste0(raw.dir, patient.ID[i], "_raw.h5")))
}

# initialize seurat objects #
for(i in 1:length(patient.ID)){
  
  # name of seurat objects
  data <- paste0("lung.", patient.ID[i],".raw")
  
  # name of count data
  count.data <- paste0("lung.", patient.ID[i],".raw.data")
  
  # initialize seurat object
  eval(parse(text=paste0(data," <- CreateSeuratObject(counts = ", count.data,", project = \"reyfman_lung\", min.cells=3, min.features=200)")))
  eval(parse(text=paste0("rm(", count.data, ")")))
  
  # add metadata: author_dataset
  eval(parse(text=paste0(data,"$author <- \"refyman_lung\"")))
  
  # add metadata: organ
  eval(parse(text=paste0(data,"$organ <- \"lung\"")))
  
  # add metadata: condition
  if(substr(patient.ID[i],1,1) == "h"){
    eval(parse(text=paste0(data,"$condition <- \"healthy\"")))
  } else if (substr(patient.ID[i],1,1) == "f"){
    eval(parse(text=paste0(data,"$condition <- \"fibrotic\"")))
  }
  
  # add metadata: patient ID
  eval(parse(text=paste0(data,"$patient.ID <- \"", patient.ID, "\"")))
  
  # add meta data: calculate missing metrics
  eval(parse(text=paste0(data,"[[\"percent.mt\"]] <- PercentageFeatureSet(", data, ", pattern = \"^MT-\")")))
}

### merge raw data ###
## merge Seurat objects
lung.raw <- merge(lung.healthy_01.raw, c(lung.healthy_02.raw, lung.healthy_03.raw, lung.healthy_04.raw,
                                         lung.healthy_05.raw, lung.healthy_06.raw, lung.healthy_07.raw,
                                         lung.healthy_08.raw, lung.fibrotic_01.raw, lung.fibrotic_02.raw,
                                         lung.fibrotic_03.raw, lung.fibrotic_04.raw, lung.fibrotic_05.raw,
                                         lung.fibrotic_06.raw, lung.fibrotic_07.raw, lung.fibrotic_08.raw))

### save merged raw data ###
saveRDS(lung.raw, file = "reyfman_lung_raw.rds")

### read filtered data ###
filtered.dir <- "/Users/Niklas/Documents/Bioinformatics/MacrophageAtlas/reyfman_filtered/"
patient.ID <- paste0(c(rep("healthy_",8),rep("fibrotic_",8)), rep(sprintf("%02d", 1:8),2))

# load filtered lung samples separately #
for(i in 1:length(patient.ID)){
  assign(paste0("lung.", patient.ID[i],".filtered.data"), Read10X_h5(paste0(filtered.dir, patient.ID[i], "_filtered.h5")))
}

# initialize seurat objects #
for(i in 1:length(patient.ID)){
  
  # name of seurat objects
  data <- paste0("lung.", patient.ID[i],".filtered")
  
  # name of count data
  count.data <- paste0("lung.", patient.ID[i],".filtered.data")
  
  # initialize seurat object
  eval(parse(text=paste0(data," <- CreateSeuratObject(counts = ", count.data,", project = \"reyfman_lung\", min.cells=3, min.features=200)")))
  eval(parse(text=paste0("rm(", count.data, ")")))
  
  # add metadata: author_dataset
  eval(parse(text=paste0(data,"$author <- \"refyman_lung\"")))
  
  # add metadata: organ
  eval(parse(text=paste0(data,"$organ <- \"lung\"")))
  
  # add metadata: condition
  if(substr(patient.ID[i],1,1) == "h"){
    eval(parse(text=paste0(data,"$condition <- \"healthy\"")))
  } else if (substr(patient.ID[i],1,1) == "f"){
    eval(parse(text=paste0(data,"$condition <- \"fibrotic\"")))
  }
  
  # add metadata: patient ID
  eval(parse(text=paste0(data,"$patient.ID <- \"", patient.ID, "\"")))
  
  # add meta data: calculate missing metrics
  eval(parse(text=paste0(data,"[[\"percent.mt\"]] <- PercentageFeatureSet(", data, ", pattern = \"^MT-\")")))
}

### merge raw data ###
## merge Seurat objects
lung.filtered <- merge(lung.healthy_01.filtered, c(lung.healthy_02.filtered, lung.healthy_03.filtered, lung.healthy_04.filtered,
                                              lung.healthy_05.filtered, lung.healthy_06.filtered, lung.healthy_07.filtered,
                                              lung.healthy_08.filtered, lung.fibrotic_01.filtered, lung.fibrotic_02.filtered,
                                              lung.fibrotic_03.filtered, lung.fibrotic_04.filtered, lung.fibrotic_05.filtered,
                                              lung.fibrotic_06.filtered, lung.fibrotic_07.filtered, lung.fibrotic_08.filtered))

### save merged filtered data ###
saveRDS(lung.filtered, file = "reyfman_lung_filtered.rds")