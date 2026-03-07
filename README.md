# MS-26-015924-for-peer-review
ALS Brain IMC Study – Analysis Code
Authors: Momoka Hikosaka, and Gen Ohtsuki, Ph.D.
Date: March 7, 2026
This repository contains the analysis code used in the submitted manuscript (NCOMMS-26-015924). The code is provided for review purposes and reproduces the computational analyses performed in the study.
________________________________________
Repository Contents
The repository includes the following analysis modules:
1.	IMC Preprocessing (MATLAB)

o	IMC signal filtering

o	Nuclear detection correction

o	Field-of-view (FOV) correction
2.	IMC Proteomics Clustering (R)
3.	ALS Distance Analysis (MATLAB)
4.	Public Data Analysis (R)
________________________________________
Module Descriptions
1. IMC Preprocessing (MATLAB)
We provide the preprocessing code used to correct IMC signal detection errors, nuclear detection errors, and field-of-view (FOV) boundaries. After applying these preprocessing steps, a flag called Sel01 is added to each cell. Cells marked with Sel01 = 1 indicate cells retained for downstream analyses.
o	IMC signal filtering
To remove putative signal detection errors, manual upper-threshold filtering is applied. Cells with expression values above a manually selected threshold are replaced with NaN in the output table.
o	Nuclear detection correction (Matlab)
For refinement of single-cell identification with Cell Profiler v4.2.1 and histoCAT, we developed a custom MATLAB-based computational algorithm to identify multiple detections originating from the same cell, thereby maximizing the accuracy of cell-identification.
o	Field-of-view (FOV) correction (Matlab)
We also implemented a field-of-view (FOV) correction code that refines the region used for analysis in IMC spatial data, particularly for spinal cord tissue.
________________________________________
2. IMC Proteomics Clustering (R)
We provide the code used for clustering IMC-derived single-cell proteomic data.
IMC_Seurat_clustering.rmd
This script performs single-cell protein clustering of Imaging Mass Cytometry (IMC) data using the Seurat framework.
The workflow includes data loading, preprocessing, integration, dimensionality reduction, clustering, and visualization of cell populations.
IMC_subclustering.rmd
This script performs subclustering analysis of major cell populations identified from Imaging Mass Cytometry (IMC) single-cell protein data.
IMC_DEPs_heatmap.rmd
This script computes differential marker expression between ALS and control groups within each major cell type from the integrated IMC Seurat object and visualizes the results as a heatmap.
IMC_marker_expression_profiles.rmd
This script generates dot plots showing marker expression across subcelltypes for selected major cell classes from the integrated IMC Seurat object.
IMC_cell_transition.rmd
This script performs cell-state transition analysis of astrocyte subcelltypes identified from integrated IMC single-cell protein data, using Monocle 3 (https://cole-trapnell-lab.github.io/monocle3/).
 
The parameters in these scripts are currently optimized for the Precentral gyrus dataset.
For parameter settings used for other brain regions, please refer to IMC_analysis_parameters.txt.
 
________________________________________
3. ALS Distance Analysis (MATLAB)
This module contains MATLAB scripts used to analyze spatial relationships among distinct cell clusters.
The algorithm:
•	Calculates pairwise distances between cells within defined clusters
•	Evaluates whether spatial differences between clusters are statistically significant
•	Summarizes the results to determine whether significant spatial correlations are present
Running the Distance Analysis
Download all files and place them in a single directory.
Run the scripts in the following order:
ALS_Distance_Analysis_1_ALS_F_FrontalCrtx_240511_2_260306_hikosaka.m
...
ALS_Distance_Analysis_6_ALS_M_SpinalCord_240424_1_260306_hikosaka.m

After completing the above scripts, run the final script:
ALS_Distance_Analysis_Total_analysis_260306_hikosaka.m
________________________________________
4. Public Data Analysis (R)
We provide R scripts used to analyze publicly available single-cell RNA-seq
data from ALS and control samples (GEO: GSE174332; Pineda et al., Cell 2024).
scRNAseq_Seurat.rmd
This script performs Seurat-based integration and clustering analysis of publicly available single-cell RNA-seq data from ALS and control samples.
 
scRNAseq_expression.rmd
This script loads a Seurat object and performs gene expression visualization and statistical analysis.
________________________________________
Notes
This repository is provided for peer review and reproducibility purposes. Additional documentation may be added after manuscript acceptance.
