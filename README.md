## 2. IMC Proteomics Spillover Correction (MATLAB)

We provide an in-house MATLAB program for spillover-signal compensation.

The program takes single-cell data in CSV format extracted from MCD files and corrects signal spillover using the signal-overlap matrices provided by Standard BioTools, Inc.:

* `Signal Overlap Matrix_ALS_Hyperion.xlsx` for the Hyperion dataset
* `Signal Overlap Matrix_ALS_HypXTi.xlsx` for the Hyperion XTi dataset

The program outputs spillover-corrected values as CSV files with the suffix `_SOc.csv`.

### Hyperion spillover-correction code

`[spillover correction code for Hyperion].m`

This script corrects single-cell data acquired using the Hyperion platform. It should be used with:

* `Signal Overlap Matrix_ALS_Hyperion.xlsx`
* `Hyperion_PrecentralGyrus_01_ALS.csv` (demo dataset)

### Hyperion XTi spillover-correction code

`[spillover correction code for Hyperion XTi].m`

This script corrects single-cell data acquired using the Hyperion XTi platform. It should be used with:

* `Signal Overlap Matrix_ALS_HypXTi.xlsx`
* `HyperionXTi_PrecentralGyrus_13_ALS_1.csv` (demo dataset)

The corresponding raw MCD data and spillover-corrected single-cell data are available in the following folders in the [Zenodo record](https://doi.org/10.5281/zenodo.22708058):

* `IMC-Hyperion_mcd`
* `IMC-Hyperion_csv_SOc`
* `IMC-HyperionXTi_mcd`
* `IMC-HyperionXTi_csv_SOc`
