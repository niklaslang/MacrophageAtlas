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
ramacha_liver.path <- "/Users/Niklas/Documents/MASTER/scrublet/ramacha_liver/"

### load data ###
liver <- readRDS(paste0(ramacha_liver.path, "ramachandran_liver_healthy.rds"))

### patient/sample identifiers ###
patient.ID <- c("Healthy_1", "Healthy_2", "Healthy_3", "Healthy_4", "Healthy_5", "Healthy_6", "Healthy_7")

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
  plt$savefig(paste0(ramacha_liver.path, data, ".hist.png"))
  
  # add auto doublet score cut-off to seurat-object #
  eval(parse(text=paste0(data, "$scrublet_auto <- tmp[[2]]")))
  
  # add manual cut-off as control ###
  #eval(parse(text=paste0(data, "$scrublet_manual <- tmp[[1]] > 0.35")))
}

### merge samples ###
liver <- merge(liver.Healthy_1, c(liver.Healthy_2, liver.Healthy_3, liver.Healthy_4, liver.Healthy_5))

### show results: 732 doublets detected ###
table(liver$scrublet_auto)

### save data ###
saveRDS(liver, "/Users/Niklas/Documents/MASTER/scrublet/ramacha_liver/ramachandran_liver_healthy.rds")
