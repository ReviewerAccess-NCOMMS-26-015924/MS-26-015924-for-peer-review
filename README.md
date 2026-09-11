## 5. Public Data Analysis (R)

We provide R scripts used to analyze publicly available single-nucleus RNA-seq data from ALS and control samples obtained from GEO accession [GSE174332](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174332) (Pineda et al., *Cell*, 2024).

### Samples

#### ALS samples

| GEO accession | Sample ID                  |
| ------------- | -------------------------- |
| GSM5292194    | `201019_ALS_101_snRNA-B1`  |
| GSM5292195    | `201019_ALS_102_snRNA-B2`  |
| GSM5292196    | `201019_ALS_103_snRNA-B3`  |
| GSM5292197    | `201019_ALS_104_snRNA-B4`  |
| GSM5292198    | `201019_ALS_106_snRNA-B5`  |
| GSM5292199    | `201019_ALS_108_snRNA-B6`  |
| GSM5292143    | `191112_ALS_110_snRNA-B9`  |
| GSM5292144    | `191112_ALS_111_snRNA-B10` |
| GSM5292145    | `191112_ALS_112_snRNA-D1`  |
| GSM5292148    | `191112_ALS_116_snRNA-C2`  |
| GSM5292149    | `191112_ALS_117_snRNA-C3`  |
| GSM5292150    | `191112_ALS_118_snRNA-C4`  |
| GSM5292151    | `191112_ALS_120_snRNA-A1`  |
| GSM5292152    | `191112_ALS_122_snRNA-D2`  |
| GSM5292153    | `191112_ALS_124_snRNA-A3`  |
| GSM5292154    | `191112_ALS_126_snRNA-A6`  |
| GSM5292155    | `191112_ALS_127_snRNA-A4`  |

#### Control samples

| GEO accession | Sample ID                 |
| ------------- | ------------------------- |
| GSM5292174    | `191114_PN_301_snRNA-E7`  |
| GSM5292193    | `200721_PN_302_snRNA-B4`  |
| GSM5292175    | `191114_PN_303_snRNA-E8`  |
| GSM5292176    | `191114_PN_304_snRNA-E9`  |
| GSM5292177    | `191114_PN_306_snRNA-E10` |
| GSM5292178    | `191114_PN_307_snRNA-E11` |
| GSM5292180    | `191114_PN_309_snRNA-F1`  |
| GSM5292201    | `201019_PN_311_snRNA-B8`  |
| GSM5292181    | `191114_PN_317_snRNA-F2`  |
| GSM5292182    | `191114_PN_318_snRNA-F3`  |
| GSM5292183    | `191114_PN_319_snRNA-F4`  |
| GSM5292184    | `191114_PN_322_snRNA-F5`  |
| GSM5292185    | `191114_PN_323_snRNA-F6`  |
| GSM5292186    | `191114_PN_324_snRNA-F7`  |
| GSM5292187    | `191114_PN_325_snRNA-F8`  |
| GSM5292188    | `191114_PN_328_snRNA-F9`  |

### `scRNAseq_Seurat.Rmd`

This script performs Seurat-based integration and clustering of publicly available single-nucleus RNA-seq data from ALS and control samples.

### `scRNAseq_expression.Rmd`

This script loads the integrated Seurat object and performs gene-expression visualization and statistical analysis.
