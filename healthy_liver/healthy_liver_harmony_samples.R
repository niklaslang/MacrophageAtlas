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
library(data.table)

### path variables ###
ramachandran.liver.path <- "/home/s1987963/ds_group/Niklas/ramacha_liver/ramacha_liver_healthy_filtered.rds"
macparland.liver.path <- "/home/s1987963/ds_group/Niklas/macparland_liver/macparland_liver_filtered.rds"
liver.path <- "/home/s1987963/ds_group/Niklas/healthy_liver/"
harmony.samples.path <- "/home/s1987963/ds_group/Niklas/healthy_liver/harmonize_samples/"

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

### merge data ###
liver <- merge(ramachandran.liver, macparland.liver)

### feature selection ###
### 4 approaches - RUN ONLY ONE OF THEM! ###
## approach no.1: select top 2000 overall HVGs ##
liver <- FindVariableFeatures(liver, selection.method = "vst", nfeatures = 2000)
# new path variable
harmony.samples.path <- "/home/s1987963/ds_group/Niklas/healthy_liver/harmonize_samples/HVGs_overall/"

# approach no. 2: select the top 4000 HVGs of study and only consider genes that are HVGs in both studies
ramachandran.liver <- FindVariableFeatures(ramachandran.liver, selection.method = "vst", nfeatures = 4000)
macparland.liver <- FindVariableFeatures(macparland.liver, selection.method = "vst", nfeatures = 4000)
liver.HVGs <- intersect(VariableFeatures(ramachandran.liver), VariableFeatures(macparland.liver))
VariableFeatures(liver) <- liver.HVGs
# new path variable
harmony.samples.path <- "/home/s1987963/ds_group/Niklas/healthy_liver/harmonize_samples/HVGs_studies/"

## approach no.3: select the top 4000 HVGs of each sample and consider genes that are HVGs in > 2 samples ##
# create vector to store HVGs
liver.HVGs <- c()
# split liver datasets by sample
ramachandran.liver.list <- SplitObject(ramachandran.liver, split.by = "patient.ID")
macparland.liver.list <- SplitObject(macparland.liver, split.by = "patient.ID")
# calculate HVGs for each sample
for (i in 1:length(ramachandran.liver.list)) {
  ramachandran.liver.list[[i]] <- FindVariableFeatures(ramachandran.liver.list[[i]], selection.method = "vst", nfeatures = 4000)
  # add sample HVGs to HVG vector
  liver.HVGs <- append(liver.HVGs, VariableFeatures(ramachandran.liver.list[[i]]))
}
for (i in 1:length(macparland.liver.list)) {
  macparland.liver.list[[i]] <- FindVariableFeatures(macparland.liver.list[[i]], selection.method = "vst", nfeatures = 4000)
  # add sample HVGs to HVG vector
  liver.HVGs <- append(liver.HVGs, VariableFeatures(macparland.liver.list[[i]]))
}
# select only genes that are highly variable in >2 samples
liver.HVGs <- unique(liver.HVGs[which(table(liver.HVGs) > 2)])
# set HVGs of merged data
VariableFeatures(liver) <- liver.HVGs
# new path variable
harmony.samples.path <- "/home/s1987963/ds_group/Niklas/healthy_liver/harmonize_samples/HVGs_samples/"

## approach no. 4: select the top 4000 HVGs of each sample, consider only genes ##
# that are HVGs in at least two samples AND that are HVGs in both studies
# split ramachandran data by sample
ramachandran.liver.list <- SplitObject(ramachandran.liver, split.by = "patient.ID")
# create empty vector for HVGs
ramachandran.HVGs <- c()
# calculate HVGs for each sample
for (i in 1:length(ramachandran.liver.list)) {
  ramachandran.liver.list[[i]] <- FindVariableFeatures(ramachandran.liver.list[[i]], selection.method = "vst", nfeatures = 4000)
  # add HVGs to dataset specific HVG vector
  ramachandran.HVGs <- append(ramachandran.HVGs, VariableFeatures(ramachandran.liver.list[[i]]))
}
# select only genes that are HVG in > 1 sample
ramachandran.HVGs <- unique(ramachandran.HVGs[which(table(ramachandran.HVGs) > 1)])
# split macparland data by sample
macparland.liver.list <- SplitObject(macparland.liver, split.by = "patient.ID")
# create empty vector for HVGs
macparland.HVGs <- c()
# calculate HVGs for each sample
for (i in 1:length(macparland.liver.list)) {
  macparland.liver.list[[i]] <- FindVariableFeatures(macparland.liver.list[[i]], selection.method = "vst", nfeatures = 4000)
  # add HVGs to dataset specific HVG vector
  macparland.HVGs <- append(macparland.HVGs, VariableFeatures(macparland.liver.list[[i]]))
}
# select only genes that are HVG in > 1 sample
macparland.HVGs <- unique(macparland.HVGs[which(table(macparland.HVGs) > 1)])
# select only genes that are HVG in both studies
liver.HVGs <- intersect(ramachandran.HVGs, macparland.HVGs)
# set HVGs of merged data
VariableFeatures(liver) <- liver.HVGs
# new path variable
harmony.samples.path <- "/home/s1987963/ds_group/Niklas/healthy_liver/harmonize_samples/HVGs_samples_studies/"

### scale data ###
liver <- ScaleData(liver) # uncorrected
# liver <- ScaleData(liver, vars.to.regress = "nFeature_RNA") # adjusted for nFeatures
# liver <- ScaleData(liver, vars.to.regress = c("nCount_RNA", "percent.mt")) # adjusted for nCounts + percent.mt

### dimensionality reduction: PCA ###
liver <- RunPCA(liver, features = VariableFeatures(object = liver))

# elbow plot
pca.elbow.plot <- ElbowPlot(liver, ndims = 50, reduction = "pca")
png(paste0(harmony.samples.path,"pca.elbow.plot.png"), width=1000,height=600,units="px")
print(pca.elbow.plot)
dev.off()

# visualize PCs
PC1_2.samples.plot <- DimPlot(object = liver, reduction = "pca", pt.size = .1, group.by = "patient.ID")
png(paste0(harmony.samples.path,"PC1_2.samples.png"), width=1000,height=1000,units="px")
print(PC1_2.samples.plot)
dev.off()

PC1_2.studies.plot <- DimPlot(object = liver, reduction = "pca", pt.size = .1, group.by = "study")
png(paste0(harmony.samples.path,"PC1_2.studies.png"), width=1000,height=1000,units="px")
print(PC1_2.studies.plot)
dev.off()

# save genes making up the PCs  
sink(paste0(harmony.samples.path, "PC_genes.txt"))
print(liver[["pca"]], dims = 1:50, nfeatures = 20)
sink()

### integration with harmony ###
# harmonize samples
liver.harmony <- liver %>% RunHarmony("patient.ID", theta = 2, reduction.save = "harmony_theta2", plot_convergence = TRUE) # harmonize all 12 samples independently

## harmony elbow plot ##
harmony.elbow.plot <- ElbowPlot(liver.harmony, ndims = 50, reduction = "harmony_theta2")
png(paste0(harmony.samples.path,"harmony_theta2.elbow.plot.png"), width=1000,height=600,units="px")
print(harmony.elbow.plot)
dev.off()

## explore harmony coordinates ##
## visualize PCs ##
harmony.PC1_2.samples.plot <- DimPlot(object = liver.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "patient.ID")
png(paste0(harmony.samples.path,"harmony_theta2.PC1_2.samples.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.samples.plot)
dev.off()

harmony.PC1_2.studies.plot <- DimPlot(object = liver.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "study")
png(paste0(harmony.samples.path,"harmony_theta2.PC1_2.studies.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.studies.plot)
dev.off()

### visualize harmony coordinates ###
lower_limit <- seq(1,43,6)
upper_limit <- seq(6,48,6)
for(i in 1:8){
  png(paste0(harmony.path, "harmonyPC_", lower_limit[i], "_", upper_limit[i], ".png"), width=7,height=5,res=300,units="in")
  print(DimHeatmap(harmony.samples.path, reduction = "harmony", nfeatures = 20, dims = lower_limit[i]:upper_limit[i], cells = 500, balanced = TRUE))
  dev.off()
}

### save genes making up the PCs ### 
sink(paste0(harmony.samples.path, "harmony_PC_genes.txt"))
print(liver.harmony[["harmony_theta2"]], dims = 1:50, nfeatures = 20)
sink()

### explore different numbers of harmony PCs ###
# MNP genes
MNP.genes <- c("CD14", "FCGR3A","CSF1R", "CD68", "LYZ", 
               "CCR2", "CX3CR1", "S100A8", "MARCO", "MNDA", 
               "TIMD4", "CD163", "C1QB", "F13A1",	"MCEMP1",	
               "TREM2",	"CD9",	"VCAN")
# cDC 1 genes
cDC1.genes <- c("XCR1", "CLEC10A", "FCER1A", "CLEC9A")

# cDC 2 genes 
cDC2.genes <- c("CD1C","HLA-DRA",	"HLA-DRB1",	"CCR7",	"ITGAX",
                "CD1E")

# pDC genes 
pDC.genes <- c("LILRA4", "IRF8",	"IRF7",	"CLEC4C")

# T cells
Tcell.genes <- c("CD3D",	"CD4",	"CD8A", "IL7R", "CD2", 
                 "CCR7",	"LTB", "TRAC")

# NK cells
NK.genes <- c("GZMA",	"NKG7",	"XCL1",	"FGFBP2",	"GNLY",	
              "GZMK",	"KLRF1",	"KLRC1")

# B cells
Bcell.genes <- c("CD79A", "CD79B",	"MS4A1")

# plasma cells
plasmacell.genes <- c("IGHG1", "MZB1",	"IGKC",	"JCHAIN",	"JSRP1")

# mast cells
mastcell.genes <- c("TPSB2",	"TPSAB1",	"CPA3",	"MS4A2")

# erys
ery.genes <- c("HBA1", "HBB")

# epithelial cells
epithelial.genes <- c("EPCAM",	"KRT19", "ALB", "TF",	"HNF4A",	"SFTPC")

# endothelial cells
endothelial.genes <- c("KDR",	"CD34", "VWF", "CLDN5", "PECAM1", 
                       "ICAM2",	"PDPN",	"CLEC14A",	"MMRN1")

# mesenchymal cells
mesenchymal.genes <- c("COL1A1",	"COL3A1", "ACTA2", "MYH11",	"PDGFRA", 
                       "PDGFRB", "RGS5",	"MSLN",	"DCN")

# proliferating cells
proliferation.genes <- c("MKI67",	"TOP2A")

# visualize batch effect and marker gene expression
dims <- c(10,16,22,32,41,50) # harmonize samples
for(d in dims){
  
  # create folder
  dir.create(paste0(harmony.samples.path, "dim", d, "_annotation"))
  dim.path <- paste0(harmony.samples.path, "dim", d, "_annotation/")
  
  # run UMAP
  liver.harmony <- RunUMAP(liver.harmony, reduction = "harmony_theta2", dims = 1:d, seed.use=1)
  
  ## visualize batch effect 
  # no.1 overview
  sample.batch.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.1)
  study.batch.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "study", pt.size = 0.1)
  batch1.plot <- sample.batch.plot + study.batch.plot
  png(paste0(dim.path, "UMAP_dim", d, ".batch.png"), width=1500, height=600, units="px")
  print(umap.batch1.plot)
  dev.off()
  
  # no.2 - colored by patients
  batch2.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".patient.batch1.png"), width=1600, height=1000, units="px")
  print(batch2.plot)
  dev.off()
  
  # no.3 - split by patients
  batch3.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "study", split.by = "patient.ID", pt.size = 0.01, ncol = 6)
  png(paste0(dim.path, "UMAP_dim", d, ".patient.batch2.png"), width=1800, height=600, units="px")
  print(batch3.plot)
  dev.off()
  
  # no.4 - colored by study
  batch4.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "study", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".study.batch1.png"), width=1600, height=1000, units="px")
  print(batch4.plot)
  dev.off()
  
  # no.5 - split by study
  batch5.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "study", split.by = "study", pt.size = 0.01, ncol = 2)
  png(paste0(dim.path, "UMAP_dim", d, ".study.batch2.png"), width=1800, height=900, units="px")
  print(batch5.plot)
  dev.off()
  
  ## QC metrics at cluster level ##
  cluster.qc.heatmap <- FeaturePlot(liver.harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "cluster.qc.heatmap.png"), width=1500,height=500,units="px")
  print(cluster.qc.heatmap)
  dev.off()
  
  # feature plot with immune cell marker
  immunecell.markers <- FeaturePlot(liver.harmony, features = c("PTPRC"), pt.size = 0.5, ncol = 1) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "immunecell.markers.png"), width=1000,height=1000,units="px")
  print(immunecell.markers)
  dev.off()
  
  # feature plot with MNP markers
  MNP.markers <- FeaturePlot(liver.harmony, features = MNP.genes, pt.size = 0.5, ncol = 6) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "MNP.markers.png"), width=1800,height=900,units="px")
  print(MNP.markers)
  dev.off()
  
  # feature plot with cDC 1 markers
  cDC1.markers <- FeaturePlot(liver.harmony, features = cDC1.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "cDC1.markers.png"), width=1600,height=400,units="px")
  print(cDC1.markers)
  dev.off()
  
  # feature plot with cDC 2 markers
  cDC2.markers <- FeaturePlot(liver.harmony, features = cDC2.genes, pt.size = 0.5, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path, "UMAP_dim", d, "cDC2.markers.png"), width=1200,height=800,units="px")
  print(cDC2.markers)
  dev.off()
  
  # feature plot with pDC markers
  pDC.markers <- FeaturePlot(liver.harmony, features = pDC.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "pDC.markers.png"), width=1600,height=400,units="px")
  print(pDC.markers)
  dev.off()
  
  # feature plot with T cell markers
  Tcell.markers <- FeaturePlot(liver.harmony, features = Tcell.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "Tcell.markers.png"), width=1600,height=800,units="px")
  print(Tcell.markers)
  dev.off()
  
  # feature plot with NK cell markers
  NK.markers <- FeaturePlot(liver.harmony, features = NK.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "NK.markers.png"), width=1600,height=800,units="px")
  print(NK.markers)
  dev.off()
  
  # feature plot with B cell markers
  Bcell.markers <- FeaturePlot(liver.harmony, features = Bcell.genes, pt.size = 0.5, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "Bcell.markers.png"), width=1200,height=400,units="px")
  print(Bcell.markers)
  dev.off()
  
  # feature plot with plasma cell markers
  plasmacell.markers <- FeaturePlot(liver.harmony, features = plasmacell.genes, pt.size = 0.5, ncol = 5) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "plasmacell.markers.png"), width=1500,height=300,units="px")
  print(plasmacell.markers)
  dev.off()
  
  # feature plot with mast cell markers
  mastcell.markers <- FeaturePlot(liver.harmony, features = mastcell.genes, pt.size = 0.5, ncol = 4) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"mastcell.markers.png"), width=1600,height=400,units="px")
  print(mastcell.markers)
  dev.off()
  
  # feature plot with epithelial cell markers
  epithelial.markers <- FeaturePlot(liver.harmony, features = epithelial.genes, pt.size = 0.5, ncol = 5) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"epithelial.markers.png"), width=1500,height=400,units="px")
  print(epithelial.markers)
  dev.off()
  
  # feature plot with endothelial cell markers
  endothelial.markers <- FeaturePlot(liver.harmony, features = endothelial.genes, pt.size = 0.5, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"endothelial.markers.png"), width=1200,height=1200,units="px")
  print(endothelial.markers)
  dev.off()
  
  # feature plot with mesenchymal cell markers
  mesenchymal.markers <- FeaturePlot(liver.harmony, features = mesenchymal.genes, pt.size = 0.5, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"mesenchymal.markers.png"), width=1200,height=1200,units="px")
  print(mesenchymal.markers)
  dev.off()
  
  # feature plot with erythroid cell markers
  ery.markers <- FeaturePlot(liver.harmony, features = ery.genes, pt.size = 0.5, ncol = 2) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"ery.markers.png"), width=800,height=400,units="px")
  print(ery.markers)
  dev.off()
  
  # proliferating cells
  proliferation.markers <- FeaturePlot(liver.harmony, features = proliferation.genes, pt.size = 0.5, ncol = 2) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d,"proliferation.markers.png"), width=800,height=400,units="px")
  print(proliferation.markers)
  dev.off()
  
}

### clustering ###
dims <- c(50)
res <- seq(0.1,2.0,0.1)
for(d in dims){
  
  # set path variable
  dim.path <- paste0(harmony.samples.path, "dim", d, "_annotation/")
  
  # calculate UMAP
  liver.harmony <- RunUMAP(liver.harmony, reduction = "harmony_theta2", dims=1:d, seed.use=1)
  
  # construct kNN graph
  liver.harmony <- FindNeighbors(liver.harmony, reduction = "harmony_theta2", dims = 1:d)
  
  # plot batch effect
  sample.batch.plot <- DimPlot(liver.harmony, reduction = "umap", group.by = "patient.ID", pt.size = 0.01)
  
  for (r in res) {
    # compute clusters
    liver.harmony <- FindClusters(liver.harmony, reduction = "harmony_theta2", resolution = r)
    
    # plot clusters
    cluster.plot <- DimPlot(liver.harmony, reduction = "umap", label = F, pt.size = 0.1)
    
    # create eval plot
    eval.plot <- cluster.plot + sample.batch.plot
    png(paste0(dim.path, "UMAP_dim", d, "_res", r, ".png"), width=1500, height=600, units="px")
    print(eval.plot)
    dev.off()
  }
}

### preliminary clustering ###
liver.harmony <- FindNeighbors(liver.harmony, reduction = "harmony_theta2", dims = 1:50)
liver.harmony <- FindClusters(liver.harmony, reduction = "harmony_theta2", resolution = 1.9)

# compare PTPRC expression across clusters
umap.plot <- DimPlot(liver.harmony, reduction = "umap", label = T, label.size = 6, pt.size = 0.1)
ptprc.plot <- VlnPlot(object = liver.harmony, features = c("PTPRC"), group.by = "seurat_clusters", pt.size = 0.1) + NoLegend()
immuneclusters.plot <- umap.plot + immunecell.markers - ptprc.plot + plot_layout(ncol=1, widths=c(2,1))
png(paste0(harmony.samples.path, "dim50_annotation/immuneclusters.png"), width=1800,height=1200,units="px")
print(immuneclusters.plot)
dev.off()

## compute cluster marker genes ###
liver.markers <- FindAllMarkers(liver.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
liver.top50.markers <- liver.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
write.csv(liver.markers, file = paste0(harmony.samples.path, "ALL_marker_genes.csv"))
write.csv(liver.top50.markers, file = paste0(harmony.samples.path, "top50_marker_genes.csv"))

### cell type annotation ###
cluster.annotation <- c("Healthy Liver NK1", "Healthy Liver KC 1","Healthy Liver LSEndo 1", "Healthy Liver CD8+ T cell 1",
                     "Healthy Liver CD8+ T cell 2", "Healthy Liver CD4+ T cell", "Healthy Liver NK2", "Healthy Liver PVEndo",
                     "Healthy Liver Hepatocyte 1", "Healthy Liver Hepatocyte 2", "Healthy Liver Monocyte 2", "Healthy Liver CVEndo",
                     "Healthy Liver Mesenchyme 1", "Healthy Liver Monocyte 1", "Healthy Liver HAEndo", "Healthy Liver cDC2",
                     "Healthy Liver Cholangiocyte", "Healthy Liver CD8+ T cell 3", "Healthy Liver Proliferating 1", "Healthy Liver Macrophage",
                     "Healthy Liver B cell", "Healthy Liver LSEndo 2", "Healthy Liver LSEndo 3", "Healthy Liver NK3",
                     "Healthy Liver Plasma cell 1", "Healthy Liver KC 2", "Healthy Liver Plasma cell 2", "Healthy Liver LymphEndo",
                     "Healthy Liver cDC1", "Healthy Liver LSEndo 4", "Healthy Liver pDC", "RBC",
                     "Healthy Liver Hepatocyte 3", "Healthy Liver Mesenchyme 2", "RBC/T cell doublets", "Healthy Liver NK4",
                     "Healthy Liver Proliferating 2", "Endo/mac doublet", "Healthy Liver Mast cell", "Healthy Liver CCR7+ DC")
names(cluster.annotation) <- levels(liver.harmony)
liver.harmony <- RenameIdents(liver.harmony, cluster.annotation)

# save annotated UMAP
annotated.umap.plot <- DimPlot(liver.harmony, reduction = "umap", label = T, label.size = 5, pt.size = 0.1)
png(paste0(harmony.samples.path, "dim50_annotation/UMAP_annotated.png"), width=1800,height=1200,units="px")
print(annotated.umap.plot)
dev.off()

### cell lineage annotation ###
cell.data <- data.table(barcode = colnames(liver.harmony),
                           celltype = Idents(liver.harmony))

lineage.annotation <- c("NK cell","MP","Endothelia","T cell","T cell","T cell","NK cell","Endothelia",
                        "Epithelia","Epithelia","MP","Endothelia","Mesenchyme","MP","Endothelia","MP",
                        "Epithelia","T cell","Proliferating","MP","B cell","Endothelia","Endothelia","NK cell",
                        "Plasma cell","MP","Plasma cell","Endothelia","MP","Endothelia","MP","Red Blood Cell",
                        "Epithelia","Mesenchyme","Doublet","NK cell","Proliferating","Doublet","Mast cell","MP")

lineage.data <- data.table(celltype = cluster.annotation, lineage = lineage.annotation)
meta.data <- merge(cell.data, lineage.data, by = "celltype")
meta.data <- data.frame(meta.data, row.names = meta.data$barcode)
meta.data$barcode <- NULL
meta.data$celltype <- NULL
liver.harmony <- AddMetaData(liver.harmony, meta.data, col.name = "lineage")

### save R session ###
save.image(file = "/home/s1987963/ds_group/Niklas/healthy_liver/harmonize_samples/HVGs_studies/healthy_liver_harmony.RData")

### save data ###
saveRDS(liver.harmony, "/home/s1987963/ds_group/Niklas/healthy_organs/healthy_liver_annotated.rds")
