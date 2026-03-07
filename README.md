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

The IMC proteomics data have been deposited in Zenodo (an open-access repository developed under the European OpenAIRE program and operated by CERN) under the accession number DOI: 10.5281/zenodo.18823151. Please locate the folder “IMC_preprocessing_raw_csv” and use the data contained within.

o	IMC signal filtering

To remove putative signal detection errors, manual upper-threshold filtering is applied. Cells with expression values above a manually selected threshold are replaced with NaN in the output table.

o	Nuclear detection correction (Matlab)

For refinement of single-cell identification with Cell Profiler v4.2.1 and histoCAT, we developed a custom MATLAB-based computational algorithm to identify multiple detections originating from the same cell, thereby maximizing the accuracy of cell-identification.

o	Field-of-view (FOV) correction (Matlab)

We also implemented a field-of-view (FOV) correction code that refines the region used for analysis in IMC spatial data, particularly for spinal cord tissue.

________________________________________
2. IMC Clustering (R)
We provide the code used for proteomics clustering IMC-derived single-cell proteomic data, "preprocessed_csv” files in the associated folder.

The IMC proteomics data have been deposited in Zenodo (an open-access repository developed under the European OpenAIRE program and operated by CERN) under the accession number DOI: 10.5281/zenodo.18823151. Please locate the folder “IMC_clustering_rds” and use the data contained within.

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
We provide R scripts used to analyze publicly available single-cell RNA-seq data from ALS and control samples (GEO: GSE174332; Pineda et al., Cell 2024) (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174332).

GSM5292194  201019_ALS_101_snRNA-B1

GSM5292195  201019_ALS_102_snRNA-B2

GSM5292196  201019_ALS_103_snRNA-B3

GSM5292197  201019_ALS_104_snRNA-B4

GSM5292198  201019_ALS_106_snRNA-B5

GSM5292199  201019_ALS_108_snRNA-B6

GSM5292143  191112_ALS_110_snRNA-B9

GSM5292144  191112_ALS_111_snRNA-B10

GSM5292145  191112_ALS_112_snRNA-D1

GSM5292148  191112_ALS_116_snRNA-C2

GSM5292149  191112_ALS_117_snRNA-C3

GSM5292150  191112_ALS_118_snRNA-C4

GSM5292151  191112_ALS_120_snRNA-A1

GSM5292152  191112_ALS_122_snRNA-D2

GSM5292153  191112_ALS_124_snRNA-A3

GSM5292154  191112_ALS_126_snRNA-A6

GSM5292155  191112_ALS_127_snRNA-A4


GSM5292174  191114_PN_301_snRNA-E7

GSM5292193  200721_PN_302_snRNA-B4

GSM5292175  191114_PN_303_snRNA-E8

GSM5292176  191114_PN_304_snRNA-E9

GSM5292177  191114_PN_306_snRNA-E10

GSM5292178  191114_PN_307_snRNA-E11

GSM5292180  191114_PN_309_snRNA-F1

GSM5292201  201019_PN_311_snRNA-B8

GSM5292181  191114_PN_317_snRNA-F2

GSM5292182  191114_PN_318_snRNA-F3

GSM5292183  191114_PN_319_snRNA-F4

GSM5292184  191114_PN_322_snRNA-F5

GSM5292185  191114_PN_323_snRNA-F6

GSM5292186  191114_PN_324_snRNA-F7

GSM5292187  191114_PN_325_snRNA-F8

GSM5292188  191114_PN_328_snRNA-F9


scRNAseq_Seurat.rmd

This script performs Seurat-based integration and clustering analysis of publicly available single-cell RNA-seq data from ALS and control samples.
 
scRNAseq_expression.rmd

This script loads a Seurat object and performs gene expression visualization and statistical analysis.

________________________________________
Notes
This repository is provided for peer review and reproducibility purposes. Additional documentation may be added after manuscript acceptance.
