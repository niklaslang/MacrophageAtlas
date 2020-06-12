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
macparland_liver.path <- "/Users/Niklas/Documents/MASTER/scrublet/macparland_liver/"

### load data ###
liver <- readRDS(paste0(macparland_liver.path, "macparland_liver.rds"))

### patient/sample identifiers ###
patient.ID <- levels(liver[[]]$patient.ID)

### split data for scrublet ###
for(i in 1:length(patient.ID)){
  patient.split <- paste0("liver.",patient.ID[i])
  eval(parse(text=paste0(patient.split, "<- subset(liver, subset = patient.ID == \"", patient.ID[[i]], "\")")))
}

### run scrublet on each sample separately ###
for(i in 1:length(patient.ID)){
  data <- paste0("liver.",patient.ID[i])
  
  ## remove doublets ##
  eval(parse(text=paste0("scrubbed <- scr$Scrublet(t(as.matrix(", data, "@assays$RNA@data)), expected_doublet_rate=0.1, sim_doublet_ratio=5)")))
  #scrubbed <- scr$Scrublet(t(as.matrix(liver@assays$RNA@data)), expected_doublet_rate=0.1, sim_doublet_ratio=5)
  tmp <- scrubbed$scrub_doublets()
  
  # compare doublet score distribution for simulated doublets and observed cells #
  scrubbed$plot_histogram()
  plt$savefig(paste0(macparland_liver.path, data, ".hist.png"))
  
  # add auto doublet score cut-off to seurat-object #
  eval(parse(text=paste0(data, "$scrublet_auto <- tmp[[2]]")))
  
  # add manual cut-off as control ###
  #eval(parse(text=paste0(data, "$scrublet_manual <- tmp[[1]] > 0.35")))
}

### merge samples ###
liver <- merge(liver.P1TLH, c(liver.P2TLH, liver.P3TLH, liver.P4TLH, liver.P5TLH))

### save data ###
saveRDS(liver, "/Users/Niklas/Documents/MASTER/scrublet/macparland_liver/macparland_liver.rds")
