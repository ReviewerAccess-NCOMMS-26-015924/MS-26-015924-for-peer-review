## 1. IMC Preprocessing (MATLAB)

We provide the preprocessing code used to correct IMC signal-detection errors, nuclear-detection errors, and field-of-view (FOV) boundaries.

After these preprocessing steps are applied, a flag named `Sel01` is added to each cell. Cells marked with `Sel01 = 1` are retained for downstream analyses.

### Demo data

One representative sample is included as a demonstration dataset:

* `SpinalCord_02_ALS.csv` — single-cell signal table
* `Ir_image_SpinalCord_02_ALS.tiff` — corresponding Ir-channel image

### Full dataset

The complete IMC proteomics dataset has been deposited in Zenodo, an open-access repository developed under the European OpenAIRE program and operated by CERN.

* **Zenodo record:** 10.5281/zenodo.22708058
* **DOI:** https://doi.org/10.5281/zenodo.22708058

### IMC signal filtering

To remove putative signal-detection errors, manual upper-threshold filtering is applied. Expression values above a manually selected threshold are replaced with `NaN` in the output table.

### Nuclear-detection correction

To refine single-cell identification performed using CellProfiler v4.2.1 and histoCAT, we developed a custom MATLAB-based computational algorithm to identify multiple detections originating from the same cell, thereby improving the accuracy of cell identification.

### Field-of-view correction

We also implemented a field-of-view correction algorithm to refine the regions included in the analysis of IMC spatial data, particularly for spinal cord tissue.
