library(reticulate)
use_condaenv("r-reticulate", required = TRUE)
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

### load python packages ###
scr <- import("scrublet")
plt <- import("matplotlib.pyplot")

### path variables ###
habermann_lung.path <- "/Users/Niklas/Documents/MASTER/scrublet/habermann_lung/"

### load data ###
lung <- readRDS(paste0(habermann_lung.path, "habermann_lung_fibrotic.rds"))

### patient/sample identifiers ###
patient.ID <- names(table(lung$patient.ID))[-c(1,2,3,11,12,13,14,15,16,17)]

### split data for scrublet ###
for(i in 1:length(patient.ID)){
  patient.split <- paste0("lung.",patient.ID[i])
  eval(parse(text=paste0(patient.split, "<- subset(lung, subset = patient.ID == \"", patient.ID[[i]], "\")")))
}

# remove complete seurat object
rm(lung)

### run scrublet on each sample separately ###
for(i in 1:length(patient.ID)){
  data <- paste0("lung.",patient.ID[i])
  
  ## remove doublets ##
  eval(parse(text=paste0("scrubbed <- scr$Scrublet(t(as.matrix(", data, "@assays$RNA@data)), expected_doublet_rate=0.1, sim_doublet_ratio=5)")))
  #scrubbed <- scr$Scrublet(t(as.matrix(lung@assays$RNA@data)), expected_doublet_rate=0.1, sim_doublet_ratio=5)
  tmp <- scrubbed$scrub_doublets()
  
  # compare doublet score distribution for simulated doublets and observed cells #
  scrubbed$plot_histogram()
  plt$savefig(paste0(habermann_lung.path, data, ".hist.png"))
  
  # add auto doublet score cut-off to seurat-object #
  eval(parse(text=paste0(data, "$scrublet_auto <- tmp[[2]]")))
  
  # add manual cut-off as control ###
  #eval(parse(text=paste0(data, "$scrublet_manual <- tmp[[1]] > 0.35")))
}

### merge samples ###
lung <- merge(lung.TILD001, c( lung.TILD006, lung.TILD010, lung.TILD015, lung.TILD019, 
                               lung.TILD028, lung.TILD030, lung.VUILD48, lung.VUILD53, 
                               lung.VUILD54, lung.VUILD55, lung.VUILD57, lung.VUILD58, 
                               lung.VUILD59, lung.VUILD60, lung.VUILD61, lung.VUILD62, 
                               lung.VUILD63, lung.VUILD64, lung.VUILD65))

### show results: 1013 doublets detected ###
table(lung$scrublet_auto)

### save data ###
saveRDS(lung, "/Users/Niklas/Documents/MASTER/scrublet/habermann_lung/habermann_lung_fibrotic.rds")
