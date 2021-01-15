# A Human Macrophage Atlas in Homeostasis and Fibrosis
GitHub repository accompanying my MSc dissertation with Prakash Ramachandran (@Ramachandranlab), Centre for Inflammation Research @ University Of Edinburgh, Summer 2020.

## Overview 

Mononuclear phagocytes (MPs) play a critical role in the pathogenesis of fibrosis, a conserved final pathway of most acute and chronic inflammatory conditions and a leading cause of death worldwide. 
MPs are an exceptionally heterogeneous and plastic population exhibiting distinct and even antagonistic functional phenotypes during different stages of tissue inflammation, destruction, repair and fibrosis. 
However, the conserved molecular definitions of fibrogenic MP subpop- ulations across organs remain ill defined, hindering their therapeutic targeting and the development of urgently needed anti-fibrotic treatments. 
To dissect the cellular heterogeneity of MPs and compare fibrogenic MP subpopulations across organs, here we perform integrated analyses of more than 450,000 single human cells from publicly available scRNA-seq datasets for the first time. 
Using a consistent analysis workflow, we generate integrated single-cell atlases of the human liver, lung, kidney and blood and robustly define MP subpopulations in homeostasis and fibrosis. 
We observe changes in both MP cell composition and transcriptional profiles in lung and potentially kidney in the context of fibrosis, whilst changes in MP cell composition predominate in the liver.
Specifically, our analyses revealed subsets of tissue-resident macrophages in liver and lung, which demonstrate modified composition in fibrosis: 
the proportion MRC1− alveolar macrophages and TIMD4− Kupffer cells decreases, while the proportion of their counterparts, MRC1+ alveolar macrophages and TIMD4+ Kupffer cells, increases in fibrosis. 
We further identified a conserved TREM2+ CD9+ pro-fibrotic macrophage population which is present in liver, lung and kidney in small numbers during homeostasis and expands in fibrosis. 
We further demonstrate that expression of the conserved transcriptional signature of the TREM2+ CD9+ pro-fibrotic macrophage population is specific to the tissue, suggesting that tissue-specific rather than systemic cues are required to drive its expansion in fibrosis.
Our work represents the most comprehensive comparative analyses of MP populations in human organ fibrosis to date, with important implications for the development of macrophage-specific therapies.

## Requirements

The majority of the scripts in this repository, were run on the Bioinformatics Course Servers administered by [Dr. Simon Tomlinsson](mailto:simon.tomlinson@ed.ac.uk)
bioinfmsc3.mvm.ed.ac.uk and bioinfmsc4.med.ed.ac.uk using R (v3.6.0) and the following packages:

- Seurat (v3.1.5)
- harmony (v1.0)
- edgeR (v3.28.1)

Doublet removal was performed on my local machine using R (v4.0.0) and the R interface to Python reticulate (v1.18) to run the Python (v3.6) package scrublet (v1.0). 
The results were added to the Seurat Object, which were then uploaded to the Bioinformatics Course Servers for fruther processing and analysis.

## Structure

### Processing the data from the original publications
Quality Control (QC) and doublet removal were performed at the level of individual datasets. The respective scripts can be found in the folders *firstauthor*_*organ*.

### Generating integrated atlases of single cells from healthy human organs
Scripts to generate the integrated healthy atlases for each organ can be found in the respective healthy_*organ* folder.

### Generating integrated atlases of single cells from fibrotic human organs
Scripts to generate the integrated fibrotic atlases for each organ can be found in the respective fibrotic_*organ* folder.

### Generating integrated atlases of single cells from healthy and fibrotic human organs
Scripts to generate the integrated healthy and fibrotic atlases for each organ can be found in the respective all_*organ* folder.

#### Differential abundance (DA) analysis 
Scripts can be found in the folder for the respective organ and have the suffix *_DA.R.

#### Differential gene expression (DE) analysis
Scripts can be found in the folder for the respective organ and have the suffix *_DE.R.

### Integrated analysis of mononuclear phagocytes across organs
Scripts to generate and analyse the the MP subset for each organ can be found in the respective all_*organ*_MP folder.

### Integrated analysis of mesenchymal cells across organs
Scripts to generate and analyse the the Mesenchyme subset for each organ can be found in the respective all_*organ*_mesenchyme folder.

### I do have a question, what should I do now?
Don’t hesitate to send me an email (niklas.lang@helmholtz-muenchen.de), I'll do my very best you answer them.

