### DO NOT BIOINFMSC4.MED.ED.AC.UK ###
### RUN ON PERSONAL MACHINE ###

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
raredon_lung.path <- "/Users/Niklas/Documents/MASTER/scrublet/raredon_lung/samples/"

### patient/sample identifiers ###
lung.patient.ID <- sprintf("%02d", 1:14)

### load data for scrublet ###
for(i in 1:length(lung.patient.ID)){
  data <- paste0("lung.patient",lung.patient.ID[i])
  
  # load sample
  eval(parse(text=paste0(data, "<- readRDS(\"", raredon_lung.path, data, ".rds\")")))
}

### run scrublet on each sample separately ###
for(i in 1:length(lung.patient.ID)){
  data <- paste0("lung.patient",lung.patient.ID[i])
  
  ## remove doublets ##
  eval(parse(text=paste0("scrubbed <- scr$Scrublet(t(as.matrix(", data, "@assays$RNA@data)), expected_doublet_rate=0.1, sim_doublet_ratio=5)")))
  #scrubbed <- scr$Scrublet(t(as.matrix(lung@assays$RNA@data)), expected_doublet_rate=0.1, sim_doublet_ratio=5)
  tmp <- scrubbed$scrub_doublets()
  
  # compare doublet score distribution for simulated doublets and observed cells #
  scrubbed$plot_histogram()
  plt$savefig(paste0(raredon_lung.path, data, ".hist.png"))
  
  # add auto doublet score cut-off to seurat-object #
  eval(parse(text=paste0(data, "$scrublet_auto <- tmp[[2]]")))
}

### merge samples ###
lung <- merge(lung.patient01, c(lung.patient02, lung.patient03, lung.patient04,
                                lung.patient05, lung.patient06, lung.patient07,
                                lung.patient08, lung.patient09, lung.patient10,
                                lung.patient11, lung.patient12, lung.patient13,
                                lung.patient14))

### save data ###
saveRDS(lung, "/Users/Niklas/Documents/MASTER/scrublet/raredon_lung/raredon_lung.rds")
