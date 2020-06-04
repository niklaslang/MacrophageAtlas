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
reyfman_lung.path <- "/Users/Niklas/Documents/MASTER/scrublet/reyfman_lung/"

### patient/sample identifiers ###
patient.ID <- paste0(c(rep("healthy_",8),rep("fibrotic_",8)), rep(sprintf("%02d", 1:8),2))

### load data ###
lung <- readRDS(paste0(reyfman_lung.path, "reyfman_lung_filtered.rds"))

### split data for scrublet ###
for(i in 1:length(patient.ID)){
  patient.split <- paste0("lung.",patient.ID[i])
  eval(parse(text=paste0(patient.split, "<- subset(lung, subset = patient.ID == \"", patient.ID[[i]], "\")")))
}

### run scrublet on each sample separately ###
for(i in 1:length(patient.ID)){
  data <- paste0("lung.",patient.ID[i])
  
  ## remove doublets ##
  eval(parse(text=paste0("scrubbed <- scr$Scrublet(t(as.matrix(", data, "@assays$RNA@data)), expected_doublet_rate=0.1, sim_doublet_ratio=5)")))
  #scrubbed <- scr$Scrublet(t(as.matrix(lung@assays$RNA@data)), expected_doublet_rate=0.1, sim_doublet_ratio=5)
  tmp <- scrubbed$scrub_doublets()
  
  # compare doublet score distribution for simulated doublets and observed cells #
  scrubbed$plot_histogram()
  plt$savefig(paste0(reyfman_lung.path, data, ".hist.png"))
  
  # add auto doublet score cut-off to seurat-object #
  eval(parse(text=paste0(data, "$scrublet_auto <- tmp[[2]]")))
  
  # add manual cut-off as control ###
  #eval(parse(text=paste0(data, "$scrublet_manual <- tmp[[1]] > 0.35")))
}

### merge samples ###
lung <- merge(lung.healthy_01, c(lung.healthy_02, lung.healthy_03, lung.healthy_04,
                                 lung.healthy_05, lung.healthy_06, lung.healthy_07,
                                 lung.healthy_08, lung.fibrotic_01, lung.fibrotic_02,
                                 lung.fibrotic_03, lung.fibrotic_04, lung.fibrotic_05,
                                 lung.fibrotic_06, lung.fibrotic_07, lung.fibrotic_08))

### save data ###
saveRDS(lung, "/Users/Niklas/Documents/MASTER/scrublet/reyfman_lung/reyfman_lung_filtered.rds")
