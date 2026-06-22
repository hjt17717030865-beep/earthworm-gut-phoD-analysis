# Code Availability Package

This folder contains the custom R code required to reproduce the main and supplementary figures and statistical tables for the Nature Communications submission.

## Input data

The input workbook is `Supplementary Data.xlsx`, supplied separately as supplementary data. To run the code, either:

1. place `Supplementary Data.xlsx` in this folder, next to `run_all.R`; or
2. set the environment variable `SUPPLEMENTARY_DATA_XLSX` to the full workbook path.

## How to run

```r
source("run_all.R")
```

or from a terminal:

```bash
Rscript run_all.R
```

Outputs are written to `outputs/` by default. To use another output folder, set the environment variable `OUTPUT_DIR`.

## Software

The code was tested with R 4.4.2. Required CRAN packages include `readxl`, `openxlsx`, `dplyr`, `tidyr`, `tibble`, `ggplot2`, `cowplot`, `patchwork`, `vegan`, `car`, `broom`, `emmeans`, `multcompView`, `officer`, `flextable`, `svglite`, `stringr`, `scales`, `ggpubr`, and `ggrepel`.

## Scope

This package intentionally includes only the analysis scripts needed to regenerate the paper figures and supplementary tables. Internal audit utilities, DOCX image-comparison scripts, forensic QC scripts, intermediate outputs, and archived manuscript versions are not included.

See `manifest.csv` for the mapping between scripts and figure/table outputs.
