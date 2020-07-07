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
ramacha_blood.path <- "/Users/Niklas/Documents/MASTER/scrublet/ramacha_blood/"

### load data ###
blood <- readRDS(paste0(ramacha_blood.path, "ramacha_blood_fibrotic.rds"))

### patient/sample identifiers ###
patient.ID <- names(table(blood$patient.ID))

### split data for scrublet ###
for(i in 1:length(patient.ID)){
  patient.split <- paste0("blood.",patient.ID[i])
  eval(parse(text=paste0(patient.split, "<- subset(blood, subset = patient.ID == \"", patient.ID[[i]], "\")")))
}

### run scrublet on each sample separately ###
for(i in 1:length(patient.ID)){
  data <- paste0("blood.",patient.ID[i])
  
  ## remove doublets ##
  eval(parse(text=paste0("scrubbed <- scr$Scrublet(t(as.matrix(", data, "@assays$RNA@data)), expected_doublet_rate=0.1, sim_doublet_ratio=5)")))
  #scrubbed <- scr$Scrublet(t(as.matrix(blood@assays$RNA@data)), expected_doublet_rate=0.1, sim_doublet_ratio=5)
  tmp <- scrubbed$scrub_doublets()
  
  # compare doublet score distribution for simulated doublets and observed cells #
  scrubbed$plot_histogram()
  plt$savefig(paste0(ramacha_blood.path, data, ".hist.png"))
  
  # add auto doublet score cut-off to seurat-object #
  eval(parse(text=paste0(data, "$scrublet_auto <- tmp[[2]]")))
  
  # add manual cut-off as control ###
  #eval(parse(text=paste0(data, "$scrublet_manual <- tmp[[1]] > 0.35")))
}

### merge samples ###
blood <- merge(blood.Blood_1, c(blood.Blood_2, blood.Blood_3, blood.Blood_4))

### show results: 1605 doublets detected ###
table(blood$scrublet_auto)

### save data ###
saveRDS(blood, "/Users/Niklas/Documents/MASTER/scrublet/ramacha_blood/ramacha_blood_fibrotic.rds")
