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
blood.path <- "/Users/Niklas/Documents/MASTER/scrublet/reyes_blood/"

### patient/sample identifiers ###
patient.ID <- c("C2P01H", "C2P05F", "C2P10H", "C2P13F", "C2P15H", 
                "C2P16H", "C2P19H", "P02H", "P04H", "P06F",
                "P07H", "P08H", "P09H", "P13H", "P15F",
                "P17H", "P18F", "P20H")

patient.no <- sprintf("%02d", 1:18)
blood <- readRDS(paste0(blood.path, "reyes_blood_healthy.rds"))

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
  plt$savefig(paste0(blood.path, data, ".hist.png"))
  
  # add auto doublet score cut-off to seurat-object #
  eval(parse(text=paste0(data, "$scrublet_auto <- tmp[[2]]")))
  
  # add manual cut-off as control ###
  #eval(parse(text=paste0(data, "$scrublet_manual <- tmp[[1]] > 0.35")))
}

### merge samples ###
blood <- merge(blood.C2P01H, c(blood.C2P05F, blood.C2P10H, blood.C2P13F,
                               blood.C2P15H, blood.C2P16H, blood.C2P19H,
                               blood.P02H, blood.P04H, blood.P06F,
                               blood.P07H, blood.P08H, blood.P09H,
                               blood.P13H, blood.P15F, blood.P17H,
                               blood.P18F, blood.P20H))

### save data ###
saveRDS(blood, "/Users/Niklas/Documents/MASTER/scrublet/reyes_blood/reyes_blood_healthy.rds")
