library(Seurat)
library(harmony)
library(future)
library(reticulate)
library(umap)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(patchwork)
options(future.globals.maxSize = 4000 * 1024^2)

### path variables ###
ramachandran.liver.path <- "/home/s1987963/ds_group/Niklas/ramacha_liver/ramacha_liver_healthy_filtered.rds"
macparland.liver.path <- "/home/s1987963/ds_group/Niklas/macparland_liver/macparland_liver_filtered.rds"
liver.path <- "/home/s1987963/ds_group/Niklas/healthy_liver/"
harmony.samples.path <- "/home/s1987963/ds_group/Niklas/healthy_liver/harmony/integrate_samples/"
harmony.studies.path <- "/home/s1987963/ds_group/Niklas/healthy_liver/harmony/integrate_studies/"
harmony.samples.studies.path <- "/home/s1987963/ds_group/Niklas/healthy_liver/harmony/integrate_samples_studies/"

### read files ###
ramachandran.liver <- readRDS(ramachandran.liver.path) # 45145 cells
macparland.liver <- readRDS(macparland.liver.path) # 8192 cells

### add meta data ###
ramachandran.liver$study <- "ramachandran_liver"
ramachandran.liver$cohort <- "Edinburgh"
macparland.liver$study <- "macparland_liver"
macparland.liver$cohort <- "Toronto"

### normalization ###
ramachandran.liver <- NormalizeData(ramachandran.liver, normalization.method = "LogNormalize", scale.factor = 10000)
macparland.liver <- NormalizeData(macparland.liver, normalization.method = "LogNormalize", scale.factor = 10000)

### feature selection ###
# approach no. 1: select the top 4000 HVGs from each study and only consider genes that are HVGs in both studies
ramachandran.liver <- FindVariableFeatures(ramachandran.liver, selection.method = "vst", nfeatures = 4000)
macparland.liver <- FindVariableFeatures(macparland.liver, selection.method = "vst", nfeatures = 4000)
liver.HVGs1 <- intersect(VariableFeatures(ramachandran.liver), VariableFeatures(macparland.liver))

# approach no. 2: select the top 4000 HVGs of each sample, consider only genes
# that are HVGs in at least two samples
# that are HVGs in both studies
ramachandran.liver.list <- SplitObject(ramachandran.liver, split.by = "patient.ID")
ramachandran.HVGs <- c()
for (i in 1:length(ramachandran.liver.list)) {
  ramachandran.liver.list[[i]] <- FindVariableFeatures(ramachandran.liver.list[[i]], selection.method = "vst", nfeatures = 4000)
  ramachandran.HVGs <- append(ramachandran.HVGs, VariableFeatures(ramachandran.liver.list[[i]]))
}
ramachandran.HVGs <- unique(ramachandran.HVGs[which(table(ramachandran.HVGs) > 1)])

macparland.liver.list <- SplitObject(macparland.liver, split.by = "patient.ID")
macparland.HVGs <- c()
for (i in 1:length(macparland.liver.list)) {
  macparland.liver.list[[i]] <- FindVariableFeatures(macparland.liver.list[[i]], selection.method = "vst", nfeatures = 4000)
  macparland.HVGs <- append(macparland.HVGs, VariableFeatures(macparland.liver.list[[i]]))
}
macparland.HVGs <- unique(macparland.HVGs[which(table(macparland.HVGs) > 1)])
liver.HVGs2 <- intersect(ramachandran.HVGs, macparland.HVGs)

### merge data ###
liver <- merge(ramachandran.liver, macparland.liver)
VariableFeatures(liver) <- liver.HVGs

### feature selection ###
# approach no.3: select top overall 2000 HVGS
# liver <- FindVariableFeatures(liver, selection.method = "vst", nfeatures = 2000)

### scale data ###
liver <- ScaleData(liver) # uncorrected
# liver <- ScaleData(liver, vars.to.regress = "nFeature_RNA") # adjusted for nFeatures
# liver <- ScaleData(liver, vars.to.regress = c("nCount_RNA", "percent.mt")) # adjusted for nCounts + percent.mt

### dimensionality reduction: PCA ###
liver <- RunPCA(liver, features = VariableFeatures(object = liver))

# elbow plot
pca.elbow.plot <- ElbowPlot(liver, ndims = 50, reduction = "pca")
png(paste0(liver.path,"pca.elbow.plot.png"), width=1000,height=600,units="px")
print(pca.elbow.plot)
dev.off()

# visualize PCs
PC1_2.samples.plot <- DimPlot(object = liver, reduction = "pca", pt.size = .1, group.by = "patient.ID")
png(paste0(liver.path,"PC1_2.samples.png"), width=1000,height=1000,units="px")
print(PC1_2.samples.plot)
dev.off()

PC1_2.studies.plot <- DimPlot(object = liver, reduction = "pca", pt.size = .1, group.by = "study")
png(paste0(liver.path,"PC1_2.studies.png"), width=1000,height=1000,units="px")
print(PC1_2.studies.plot)
dev.off()

# save genes making up the PCs  
sink(paste0(liver.path, "PC_genes.txt"))
print(liver[["pca"]], dims = 1:50, nfeatures = 20)
sink()

### integration with harmony ###
## three integration strategies ##
#liver.harmony <- liver %>% RunHarmony("patient.ID", plot_convergence = TRUE) # harmonize all 12 samples independently
#saveRDS(liver.harmony, paste0(harmony.samples.path, "healthy_liver_harmony.samples.rds"))
#liver.harmony <- liver %>% RunHarmony("study", plot_convergence = TRUE) # harmonize 2 studies independently
#saveRDS(liver.harmony, paste0(harmony.studies.path, "healthy_liver_harmony.studies.rds"))
liver.harmony <- liver %>% RunHarmony(c("study","patient.ID"), plot_convergence = TRUE) # harmonize 2 studies and 12 samples independently
saveRDS(liver.harmony, paste0(harmony.samples.studies.path, "healthy_liver_harmony.samples.studies.rds"))

## harmony elbow plot ##
harmony.elbow.plot <- ElbowPlot(liver.harmony, ndims = 50, reduction = "harmony")
png(paste0(harmony.path,"harmony.elbow.plot.png"), width=1000,height=600,units="px")
print(harmony.elbow.plot)
dev.off()

## explore harmony coordinates ##
## visualize PCs ##
harmony.PC1_2.samples.plot <- DimPlot(object = liver.harmony, reduction = "harmony", pt.size = .1, group.by = "patient.ID")
png(paste0(harmony.path,"harmony.PC1_2.samples.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.samples.plot)
dev.off()

harmony.PC1_2.studies.plot <- DimPlot(object = liver.harmony, reduction = "harmony", pt.size = .1, group.by = "study")
png(paste0(harmony.path,"harmony.PC1_2.studies.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.studies.plot)
dev.off()

### visualize harmony coordinates ###
lower_limit <- seq(1,43,6)
upper_limit <- seq(6,48,6)
for(i in 1:8){
  png(paste0(harmony.path, "harmonyPC_", lower_limit[i], "_", upper_limit[i], ".png"), width=7,height=5,res=300,units="in")
  print(DimHeatmap(liver.harmony, reduction = "harmony", nfeatures = 20, dims = lower_limit[i]:upper_limit[i], cells = 500, balanced = TRUE))
  dev.off()
}

### save genes making up the PCs ### 
sink(paste0(harmony.path, "PC_genes.txt"))
print(liver.harmony[["harmony"]], dims = 1:50, nfeatures = 20)
sink()

### visualize batch effect ###
# dims <- c(8,20,34,46,50) # harmonize samples and studies
# dims <- c(8,14,19,25,34,40) # harmonize studies
dims <- c(8,20,35,42,50) # harmonize samples
for(d in dims){
  liver.harmony <- RunUMAP(liver.harmony, reduction = "harmony", dims = 1:d, seed.use=1)
  sample.batch.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
  study.batch.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "study", pt.size = 0.1)
  umap.batch.plot <- sample.batch.plot + study.batch.plot
  png(paste0(harmony.samples.path, "UMAP_dim", d, ".batch.png"), width=1500, height=600, units="px")
  print(umap.batch.plot)
  dev.off()
}
f