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
menon_kidney.path <- "/Users/Niklas/Documents/MASTER/scrublet/menon_kidney/"

### load data ###
kidney <- readRDS(paste0(menon_kidney.path, "menon_kidney.rds"))

### patient/sample identifiers ###
patient.ID <- paste0(rep("GSM41919", 24), seq(41,64,1))

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
  plt$savefig(paste0(menon_kidney.path, data, ".hist.png"))
  
  # add auto doublet score cut-off to seurat-object #
  eval(parse(text=paste0(data, "$scrublet_auto <- tmp[[2]]")))
  
  # add manual cut-off as control ###
  #eval(parse(text=paste0(data, "$scrublet_manual <- tmp[[1]] > 0.35")))
}

### merge samples ###
kidney <- merge(kidney.GSM4191941, c(kidney.GSM4191942, kidney.GSM4191943, kidney.GSM4191944,
                                     kidney.GSM4191945, kidney.GSM4191946, kidney.GSM4191947,
                                     kidney.GSM4191948, kidney.GSM4191949, kidney.GSM4191950,
                                     kidney.GSM4191951, kidney.GSM4191952, kidney.GSM4191953,
                                     kidney.GSM4191954, kidney.GSM4191955, kidney.GSM4191956,
                                     kidney.GSM4191957, kidney.GSM4191958, kidney.GSM4191959,
                                     kidney.GSM4191960, kidney.GSM4191961, kidney.GSM4191962,
                                     kidney.GSM4191963, kidney.GSM4191964))
### show results:  176 doublets detected ###
table(kidney$scrublet_auto)

### save data ###
saveRDS(kidney, "/Users/Niklas/Documents/MASTER/scrublet/menon_kidney/menon_kidney.rds")