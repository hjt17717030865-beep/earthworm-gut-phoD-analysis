# Reproduce the manuscript and supplementary figures/tables.
# Place 'Supplementary Data.xlsx' in this directory, or set SUPPLEMENTARY_DATA_XLSX.

CODE_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
DATA_XLSX <- Sys.getenv("SUPPLEMENTARY_DATA_XLSX",
                        unset = file.path(CODE_ROOT, "Supplementary Data.xlsx"))
OUTPUT_DIR <- Sys.getenv("OUTPUT_DIR", unset = file.path(CODE_ROOT, "outputs"))
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(123)

if (!file.exists(DATA_XLSX)) {
  stop("Input workbook not found. Put 'Supplementary Data.xlsx' next to run_all.R, or set SUPPLEMENTARY_DATA_XLSX.")
}

manifest <- read.csv(file.path(CODE_ROOT, "manifest.csv"), stringsAsFactors = FALSE)
for (script in manifest$script) {
  message("Running: ", script)
  source(file.path(CODE_ROOT, script), local = new.env(parent = globalenv()))
}

message("Done. Outputs written to: ", OUTPUT_DIR)
