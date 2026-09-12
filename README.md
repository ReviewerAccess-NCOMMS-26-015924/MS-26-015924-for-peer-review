## 2. IMC Proteomics Spillover Correction (MATLAB)

We provide an in-house MATLAB program for spillover-signal compensation.

The program takes single-cell data in CSV format extracted from MCD files and corrects signal spillover using the signal-overlap matrices provided by Standard BioTools, Inc.:

* `Signal Overlap Matrix_ALS_Hyperion.xlsx` for the Hyperion dataset
* `Signal Overlap Matrix_ALS_HypXTi.xlsx` for the Hyperion XTi dataset

The single-cell CSV files to be corrected must be named with the suffix　`_SOc.csv`. 
The program reads each `_SOc.csv` file, applies spillover compensation, and overwrites the same file in place with the corrected values.

### Hyperion spillover-correction code

`SOcorrection_all_csv_ALS_Hyperion.m`
`Spillover_compensation_ALS_Hyperion_func.m`

This script corrects single-cell data acquired using the Hyperion platform. It should be used with:

* `Signal Overlap Matrix_ALS_Hyperion.xlsx`
* `IMC_csv/Hyperion_PrecentralGyrus_01_ALS_SOc.csv` (demo dataset)

### Hyperion XTi spillover-correction code

`SOcorrection_all_csv_ALS_HyperionXTi.m`
`Spillover_compensation_ALS_HyperionXTi_func.m`

This script corrects single-cell data acquired using the Hyperion XTi platform. It should be used with:

* `Signal Overlap Matrix_ALS_HypXTi.xlsx`
* `IMC_csv/HyperionXTi_PrecentralGyrus_13_ALS_1_SOc.csv` (demo dataset)

The corresponding raw MCD data and spillover-corrected single-cell data are available in the following folders in the [[Zenodo](https://doi.org/10.5281/zenodo.22708058)]:

* `IMC-Hyperion_mcd`
* `IMC-Hyperion_csv_SOc`
* `IMC-HyperionXTi_mcd`
* `IMC-HyperionXTi_csv_SOc`
