library(dplyr)
library(ggplot2)
library(data.table)
library(Seurat)
library(patchwork)
library(umap)
library(reticulate)

### path.variables ###
blood.path <- "/home/s1987963/ds_group/Niklas/reyes_blood/"
metadata.path <- "/home/s1987963/ds_group/Niklas/reyes_blood/scp_meta.txt"

### load data ###
blood <- readRDS(paste0(blood.path, "reyes_blood.rds"))

### load meta data ###
blood.metadata <- fread(metadata.path, stringsAsFactors = TRUE, header = TRUE)

### clean meta data ###
blood.metadata <- blood.metadata[-1] # remove sub-header
blood.metadata$NAME <- gsub(x = blood.metadata$NAME, pattern = "-", replacement = ".") # replace "-"
blood.metadata <- blood.metadata[NAME %in% colnames(blood)] # remove all cells that are not in count matrix

# convert meta data #
blood.metadata <- data.frame(blood.metadata, row.names = blood.metadata$NAME) # set rownames equal to cell labels
blood.metadata$NAME <- NULL # remove cell label column

### add meta-data ###
blood <- AddMetaData( object = blood, metadata = blood.metadata$Cohort, col.name = "condition")
blood <- AddMetaData( object = blood, metadata = blood.metadata$Cell_Type, col.name = "cell_type")
blood <- AddMetaData( object = blood, metadata = blood.metadata$Patient, col.name = "patient.ID")

#stopifnot(colnames(x = blood) == blood.metadata$NAME)
#blood[["condition"]] <- blood.metadata$Cohort
#blood[["cell_type"]] <- blood.metadata$Cell_Type

### subset data ###
blood.healthy <- subset(blood, subset = condition == "Control")

### save data ###
saveRDS(blood.healthy, file = paste(blood.path, "blood_reyes_healthy.rds"))

### QC ###
