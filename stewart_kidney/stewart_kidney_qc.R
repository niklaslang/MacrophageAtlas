library(Seurat)
library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(umap)
library(reticulate)
library(RColorBrewer)
library(data.table)

### path.variables ###
kidney.path <- "/home/s1987963/ds_group/Niklas/stewart_kidney/"
cell.metadata.path <- "/home/s1987963/ds_group/Niklas/stewart_kidney/cells.tsv.gz"
gene.metadata.path  <- "/home/s1987963/ds_group/Niklas/stewart_kidney/features.tsv.gz"

### load kidney data ###
# load counts 
kidney.counts <- Read10X(data.dir = kidney.path, unique.features = FALSE)

# load meta data
cell.metadata <- read.table(cell.metadata.path, sep = "\t", stringsAsFactors = TRUE, header = TRUE)
cell.metadata <- data.table(cell.metadata)
gene.metadata <- read.table(gene.metadata.path, sep = "\t", stringsAsFactors = TRUE, header = FALSE)
gene.metadata <- data.table(gene.metadata)

# reconstruct meta data
cell.metadata[, unique.barcodes := paste0(barcode, "-", cellkey)]
cell.metadata[, sample := ifelse(donor_organism.development_stage.ontology_label == "Carnegie stage 01", "fetal", "adult")]

# add colnames to count matrix
colnames(kidney.counts) <- cell.metadata$unique.barcodes

# add rownames to count matrix
rownames(kidney.counts) <- gene.metadata$V2
  
### create SeuratObject: 77609 cells ###
kidney <- CreateSeuratObject(counts = kidney.counts, min.cells = 3, min.features=200, names.delim = "-")

# subset meta data
cell.metadata <- cell.metadata[unique.barcodes %in% colnames(kidney)]

#convert meta data
cell.metadata <- data.frame(cell.metadata, row.names = cell.metadata$unique.barcodes)

# add meta data
kidney <- AddMetaData( object = kidney, metadata = cell.metadata$donor_organism.development_stage.ontology_label, col.name = "sample")

# calculate meta data
kidney$percent.mt <- PercentageFeatureSet(kidney, pattern = "^MT-")

