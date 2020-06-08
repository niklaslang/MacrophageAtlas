library(dplyr)
library(ggplot2)
library(Seurat)


### path.variables ###
blood.path <- "/home/s1987963/ds_group/Niklas/reyes_blood/scp_gex_matrix.csv"

# load sepsis data
blood.counts <- read.csv(file= blood.path, sep=",", header = T, row.names = 1)
row.names(blood.counts) <- gsub(x = row.names(blood.counts), pattern = "_", replacement = "-")
names(blood.counts) <- gsub(x = names(blood.counts), pattern = "_", replacement = "-")
blood.counts <- as.sparse(blood.counts)

# create SeuratObject
blood <- CreateSeuratObject(counts = blood.counts, min.cells = 3, min.features=200, names.delim = "-")

# Save the SeuratObject as rds data
saveRDS(blood, paste0("/home/s1987963/ds_group/Niklas/reyes_blood/reyes_blood.rds"))

