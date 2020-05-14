library(dplyr)
library(ggplot2)
library(Seurat)

# load sepsis data
raw_counts <- read.csv(file=paste0("/home/s1987963/processed_data/reyes_blood/healthy/scp_gex_matrix.csv"), sep=",", header = T, row.names = 1)
row.names(raw_counts) <- gsub(x = row.names(raw_counts), pattern = "_", replacement = "-")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "_", replacement = "-")
raw_counts <- as.sparse(raw_counts)

# create SeuratObject
spc_gex <- CreateSeuratObject(counts = raw_counts, min.cells = 3, min.features=200, names.delim = "-")

# Save the SeuratObject as rds data
saveRDS(spc_gex, paste0("~/Desktop/human/spc_gex_SeuratObj.rds"))

