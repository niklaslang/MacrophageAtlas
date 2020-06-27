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
wilson_kidney.path <- "/Users/Niklas/Documents/MASTER/scrublet/wilson_kidney/"

### load data ###
kidney <- readRDS(paste0(wilson_kidney.path, "wilson_kidney_diabetes.rds"))

### patient/sample identifiers ###
patient.ID <- c("GSM3823942_diabetes", "GSM3823943_diabetes", "GSM3823944_diabetes")

### split data for scrublet ###
for(i in 1:length(patient.ID)){
  patient.split <- paste0("kidney.",patient.ID[i])
  eval(parse(text=paste0(patient.split, "<- subset(kidney, subset = patient.ID == \"", patient.ID[[i]], "\")")))
}

# remove complete seurat object
rm(kidney)

### run scrublet on each sample separately ###
for(i in 1:length(patient.ID)){
  data <- paste0("kidney.",patient.ID[i])
  
  ## remove doublets ##
  eval(parse(text=paste0("scrubbed <- scr$Scrublet(t(as.matrix(", data, "@assays$RNA@data)), expected_doublet_rate=0.1, sim_doublet_ratio=5)")))
  #scrubbed <- scr$Scrublet(t(as.matrix(kidney@assays$RNA@data)), expected_doublet_rate=0.1, sim_doublet_ratio=5)
  tmp <- scrubbed$scrub_doublets()
  
  # compare doublet score distribution for simulated doublets and observed cells #
  scrubbed$plot_histogram()
  plt$savefig(paste0(wilson_kidney.path, data, ".hist.png"))
  
  # add auto doublet score cut-off to seurat-object #
  eval(parse(text=paste0(data, "$scrublet_auto <- tmp[[2]]")))
  
  # add manual cut-off as control ###
  #eval(parse(text=paste0(data, "$scrublet_manual <- tmp[[1]] > 0.4")))
}

### merge samples ###
kidney <- merge(kidney.GSM3823942_diabetes, c(kidney.GSM3823943_diabetes, kidney.GSM3823944_diabetes))

### show results: 127 doublets detected ###
table(kidney$scrublet_auto)
#table(kidney$scrublet_manual)

### save data ###
saveRDS(kidney, "/Users/Niklas/Documents/MASTER/scrublet/wilson_kidney/wilson_kidney_diabetes.rds")
