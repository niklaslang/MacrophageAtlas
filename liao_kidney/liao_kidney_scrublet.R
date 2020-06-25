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
liao_kidney.path <- "/Users/Niklas/Documents/MASTER/scrublet/liao_kidney/"


### load data ###
kidney <- readRDS(paste0(liao_kidney.path, "liao_kidney.rds"))

### patient/sample identifiers ###
patient.ID <- c("GSM4145204", "GSM4145205", "GSM4145206")

### split data for scrublet ###
for(i in 1:length(patient.ID)){
  patient.split <- paste0("kidney.",patient.ID[i])
  eval(parse(text=paste0(patient.split, "<- subset(kidney, subset = patient.ID == \"", patient.ID[[i]], "\")")))
}

# remove old seurat object
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
  plt$savefig(paste0(liao_kidney.path, data, ".hist.png"))
  
  # add auto doublet score cut-off to seurat-object #
  eval(parse(text=paste0(data, "$scrublet_auto <- tmp[[2]]")))
  
  # add manual cut-off as control ###
  #eval(parse(text=paste0(data, "$scrublet_manual <- tmp[[1]] > 0.35")))
}

### merge samples ###
kidney <- merge(kidney.GSM4145204, c(kidney.GSM4145205, kidney.GSM4145206))

### show results: 14 doublets detected ###
table(kidney$scrublet_auto)

### save data ###
saveRDS(kidney, "/Users/Niklas/Documents/MASTER/scrublet/liao_kidney/liao_kidney.rds")
