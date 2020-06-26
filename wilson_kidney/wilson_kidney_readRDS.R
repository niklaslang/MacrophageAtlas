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
kidney.path <- "/home/s1987963/ds_group/Niklas/wilson_kidney/"
metadata.path <- "/home/s1987963/ds_group/Niklas/stewart_kidney/genes.tsv.gz"

# patient IDs
patient.ID <- c("GSM3823939_control", "GSM3823940_control", "GSM3823941_control",
                "GSM3823942_diabetes", "GSM3823943_diabetes", "GSM3823944_diabetes")

### load count matrices ###
for(i in 1:length(patient.ID)){
  assign(paste0("kidney.",patient.ID[i],".counts"), readRDS(paste0(kidney.path, patient.ID[i], ".rds")))
}

### load meta data ###
gene.metadata <- read.table(metadata.path, sep = "\t", stringsAsFactors = TRUE, header = TRUE)
gene.metadata <- data.table(gene.metadata)

### set exon/intron matrix 
kidney.GSM3823939_control.counts <- kidney.GSM3823939_control.counts[["umicount"]][["inex"]][["all"]]
### convert count matrix ###
kidney.GSM3823939_control.counts <- as.matrix(kidney.GSM3823939_control.counts)
### clean count data ###
kidney.GSM3823939_control.counts <- subset(kidney.GSM3823939_control.counts, rownames(kidney.GSM3823939_control.counts) %in% gene.metadata$featurekey)
### subset gene names ###
genes.GSM3823939_control <- gene.metadata[featurekey %in% rownames(kidney.GSM3823939_control.counts)]
### update gene names ###
stopifnot(genes.GSM3823939_control$featurekey == rownames(kidney.GSM3823939_control.counts))
### set gene names ###
rownames(kidney.GSM3823939_control.counts) <- genes.GSM3823939_control$featurename 
### create Seurat Object ###
kidney.GSM3823939_control <- CreateSeuratObject( counts = kidney.GSM3823939_control.counts, project = "kidney", min.cells = 3, min.features = 200)
### save Seurat Object ###
saveRDS(kidney.GSM3823939_control, paste0(kidney.path, "kidney.GSM3823939_control.rds"))
### remove data ###
rm(kidney.GSM3823939_control.counts)

### set exon/intron matrix 
kidney.GSM3823940_control.counts <- kidney.GSM3823940_control.counts[["umicount"]][["inex"]][["all"]]
### clean count data ###
idx <- which(rownames(kidney.GSM3823940_control.counts) %in% gene.metadata$featurekey)
kidney.GSM3823940_control.counts <- kidney.GSM3823940_control.counts[idx, , drop=F]
### subset gene names ###
genes.GSM3823940_control <- gene.metadata[featurekey %in% rownames(kidney.GSM3823940_control.counts)]
### update gene names ###
stopifnot(genes.GSM3823940_control$featurekey == rownames(kidney.GSM3823940_control.counts))
### set gene names ###
rownames(kidney.GSM3823940_control.counts) <- genes.GSM3823940_control$featurename 
### create Seurat Object ###
kidney.GSM3823940_control <- CreateSeuratObject( counts = kidney.GSM3823940_control.counts, project = "kidney", min.cells = 3, min.features = 200)
### save Seurat Object ###
saveRDS(kidney.GSM3823940_control, paste0(kidney.path, "kidney.GSM3823940_control.rds"))
### remove data ###
rm(kidney.GSM3823940_control.counts)

### set exon/intron matrix 
kidney.GSM3823941_control.counts <- kidney.GSM3823941_control.counts[["umicount"]][["inex"]][["all"]]
### convert count matrix ###
kidney.GSM3823941_control.counts <- as.matrix(kidney.GSM3823941_control.counts)
### clean count data ###
kidney.GSM3823941_control.counts <- subset(kidney.GSM3823941_control.counts, rownames(kidney.GSM3823941_control.counts) %in% gene.metadata$featurekey)
### subset gene names ###
genes.GSM3823941_control <- gene.metadata[featurekey %in% rownames(kidney.GSM3823941_control.counts)]
### update gene names ###
stopifnot(genes.GSM3823941_control$featurekey == rownames(kidney.GSM3823941_control.counts))
### set gene names ###
rownames(kidney.GSM3823941_control.counts) <- genes.GSM3823941_control$featurename 
### create Seurat Object ###
kidney.GSM3823941_control <- CreateSeuratObject( counts = kidney.GSM3823941_control.counts, project = "kidney", min.cells = 3, min.features = 200)
### save Seurat Object ###
saveRDS(kidney.GSM3823941_control, paste0(kidney.path, "kidney.GSM3823941_control.rds"))
### remove data ###
rm(kidney.GSM3823941_control.counts)

### set exon/intron matrix 
kidney.GSM3823942_diabetes.counts <- kidney.GSM3823942_diabetes.counts[["umicount"]][["inex"]][["all"]]
### convert count matrix ###
kidney.GSM3823942_diabetes.counts <- as.matrix(kidney.GSM3823942_diabetes.counts)
### clean count data ###
kidney.GSM3823942_diabetes.counts <- subset(kidney.GSM3823942_diabetes.counts, rownames(kidney.GSM3823942_diabetes.counts) %in% gene.metadata$featurekey)
### subset gene names ###
genes.GSM3823942_diabetes <- gene.metadata[featurekey %in% rownames(kidney.GSM3823942_diabetes.counts)]
### update gene names ###
stopifnot(genes.GSM3823942_diabetes$featurekey == rownames(kidney.GSM3823942_diabetes.counts))
### set gene names ###
rownames(kidney.GSM3823942_diabetes.counts) <- genes.GSM3823942_diabetes$featurename 
### create Seurat Object ###
kidney.GSM3823942_diabetes <- CreateSeuratObject( counts = kidney.GSM3823942_diabetes.counts, project = "kidney", min.cells = 3, min.features = 200)
### save Seurat Object ###
saveRDS(kidney.GSM3823942_diabetes, paste0(kidney.path, "kidney.GSM3823942_diabetes.rds"))
### remove data ###
rm(kidney.GSM3823942_diabetes.counts)

### set exon/intron matrix 
kidney.GSM3823943_diabetes.counts <- kidney.GSM3823943_diabetes.counts[["umicount"]][["inex"]][["all"]]
### convert count matrix ###
kidney.GSM3823943_diabetes.counts <- as.matrix(kidney.GSM3823943_diabetes.counts)
### clean count data ###
kidney.GSM3823943_diabetes.counts <- subset(kidney.GSM3823943_diabetes.counts, rownames(kidney.GSM3823943_diabetes.counts) %in% gene.metadata$featurekey)
### subset gene names ###
genes.GSM3823943_diabetes <- gene.metadata[featurekey %in% rownames(kidney.GSM3823943_diabetes.counts)]
### update gene names ###
stopifnot(genes.GSM3823943_diabetes$featurekey == rownames(kidney.GSM3823943_diabetes.counts))
### set gene names ###
rownames(kidney.GSM3823943_diabetes.counts) <- genes.GSM3823943_diabetes$featurename 
### create Seurat Object ###
kidney.GSM3823943_diabetes <- CreateSeuratObject( counts = kidney.GSM3823943_diabetes.counts, project = "kidney", min.cells = 3, min.features = 200)
### save Seurat Object ###
saveRDS(kidney.GSM3823943_diabetes, paste0(kidney.path, "kidney.GSM3823943_diabetes.rds"))
### remove data ###
rm(kidney.GSM3823943_diabetes.counts)

### set exon/intron matrix 
kidney.GSM3823944_diabetes.counts <- kidney.GSM3823944_diabetes.counts[["umicount"]][["inex"]][["all"]]
### convert count matrix ###
kidney.GSM3823944_diabetes.counts <- as.matrix(kidney.GSM3823944_diabetes.counts)
### clean count data ###
kidney.GSM3823944_diabetes.counts <- subset(kidney.GSM3823944_diabetes.counts, rownames(kidney.GSM3823944_diabetes.counts) %in% gene.metadata$featurekey)
### subset gene names ###
genes.GSM3823944_diabetes <- gene.metadata[featurekey %in% rownames(kidney.GSM3823944_diabetes.counts)]
### update gene names ###
stopifnot(genes.GSM3823944_diabetes$featurekey == rownames(kidney.GSM3823944_diabetes.counts))
### set gene names ###
rownames(kidney.GSM3823944_diabetes.counts) <- genes.GSM3823944_diabetes$featurename 
### create Seurat Object ###
kidney.GSM3823944_diabetes <- CreateSeuratObject( counts = kidney.GSM3823944_diabetes.counts, project = "kidney", min.cells = 3, min.features = 200)
### save Seurat Object ###
saveRDS(kidney.GSM3823944_diabetes, paste0(kidney.path, "kidney.GSM3823944_diabetes.rds"))
### remove data ###
rm(kidney.GSM3823944_diabetes.counts)

