# ---- Portable path configuration ----
if (!exists("CODE_ROOT", inherits = FALSE)) {
  this_file <- tryCatch(normalizePath(dirname(sys.frame(1)$ofile), winslash = "/", mustWork = FALSE),
                        error = function(e) getwd())
  CODE_ROOT <- normalizePath(file.path(this_file, ".."), winslash = "/", mustWork = FALSE)
}
DATA_XLSX <- Sys.getenv("SUPPLEMENTARY_DATA_XLSX",
                        unset = file.path(CODE_ROOT, "Supplementary Data.xlsx"))
OUTPUT_DIR <- Sys.getenv("OUTPUT_DIR", unset = file.path(CODE_ROOT, "outputs"))
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(123)
if (!file.exists(DATA_XLSX)) {
  stop("Input workbook not found. Put 'Supplementary Data.xlsx' next to run_all.R, or set SUPPLEMENTARY_DATA_XLSX.")
}
if (requireNamespace("dplyr", quietly = TRUE)) {
  select <- dplyr::select
}

# ============================================================
# Supplementary Table S5
# Amplicon sequencing summary and compositional metrics
# with phoD DADA2 read tracking
# ============================================================

pkg_needed <- c("readxl", "dplyr", "tidyr", "stringr", "tibble", "openxlsx", "readr")
pkg_to_install <- pkg_needed[!sapply(pkg_needed, requireNamespace, quietly = TRUE)]
if (length(pkg_to_install) > 0) install.packages(pkg_to_install)

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(openxlsx)
  library(readr)
})

# -------------------------
# 1) Paths
# -------------------------
out_dir <- OUTPUT_DIR
data_xlsx <- DATA_XLSX
main_xlsx  <- data_xlsx
main_sheet <- "Master_all_data"
out_xlsx   <- file.path(out_dir, "Supplementary_Table_S5_amplicon_summary_compositional_metrics.xlsx")

# 2) Helper functions
# -------------------------
normalise_plotid <- function(x) {
  x %>%
    as.character() %>%
    trimws() %>%
    str_replace("\\.fastq(\\.gz)?$", "") %>%
    str_replace("_R1.*$", "") %>%
    str_replace("_R2.*$", "") %>%
    str_replace("^Soil__", "") %>%
    str_replace("^Gut__", "") %>%
    str_replace("^([0-9]+y)EB", "\\1") %>%
    str_replace("^([0-9]+y)E", "\\1") %>%
    str_replace("^3y", "5y")
}

parse_duration <- function(plotid) {
  out <- str_extract(plotid, "^\\d+y")
  out <- ifelse(out == "3y", "5y", out)
  factor(out, levels = c("5y", "8y", "10y"))
}

parse_regime <- function(plotid) {
  out <- case_when(
    str_detect(plotid, "NPKOM") ~ "NPKOM",
    str_detect(plotid, "NPK") ~ "NPK",
    str_detect(plotid, "CK") ~ "CK",
    TRUE ~ NA_character_
  )
  factor(out, levels = c("CK", "NPK", "NPKOM"))
}

parse_rep <- function(plotid) {
  suppressWarnings(as.integer(str_extract(plotid, "(?<=-)\\d+$")))
}

make_16S_qc <- function(mat, compartment) {
  tibble(
    SampleID = rownames(mat),
    PlotID = normalise_plotid(rownames(mat)),
    Duration = parse_duration(PlotID),
    Regime = parse_regime(PlotID),
    Replicate = parse_rep(PlotID),
    Compartment = compartment,
    Target = "16S",
    Raw_reads = NA_real_,
    Filtered_reads = rowSums(mat, na.rm = TRUE),
    Nonchim_reads = rowSums(mat, na.rm = TRUE),
    ASV_richness = rowSums(mat > 0, na.rm = TRUE),
    Notes = "16S raw reads not available here; non-chimeric reads calculated from seqtab row sums."
  )
}

find_col <- function(df, patterns) {
  nm <- names(df)
  hit <- nm[str_detect(tolower(nm), paste(patterns, collapse = "|"))]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

read_phoD_tracking <- function(xlsx_path, sheet_name, compartment) {
  tr <- readxl::read_excel(xlsx_path, sheet = sheet_name)
  
  sample_col  <- find_col(tr, c("sample", "sampleid", "sample_id", "filename", "file"))
  input_col   <- find_col(tr, c("input", "raw", "in_reads", "input_reads"))
  filt_col    <- find_col(tr, c("filtered", "filter"))
  nonchim_col <- find_col(tr, c("nonchim", "non_chim", "nonchimeric", "nochim"))
  
  if (is.na(sample_col)) stop("Cannot identify sample column in sheet: ", sheet_name)
  if (is.na(input_col))  stop("Cannot identify raw/input reads column in sheet: ", sheet_name)
  
  tr %>%
    transmute(
      SampleID       = as.character(.data[[sample_col]]),
      PlotID         = normalise_plotid(SampleID),
      Duration       = parse_duration(PlotID),
      Regime         = parse_regime(PlotID),
      Replicate      = parse_rep(PlotID),
      Compartment    = compartment,
      Target         = "phoD",
      Raw_reads      = as.numeric(.data[[input_col]]),
      Filtered_reads = if (!is.na(filt_col))    as.numeric(.data[[filt_col]])    else NA_real_,
      Nonchim_reads  = if (!is.na(nonchim_col)) as.numeric(.data[[nonchim_col]]) else NA_real_,
      ASV_richness   = NA_real_,
      Notes          = "phoD read counts imported from DADA2 reads tracking table."
    )
}

make_phoD_richness <- function(xlsx_path, sheet_name, compartment) {
  raw <- readxl::read_excel(xlsx_path, sheet = sheet_name) %>%
    as.data.frame()
  asv_cols <- setdiff(colnames(raw), meta_cols)
  mat <- as.matrix(sapply(raw[, asv_cols, drop = FALSE], as.numeric))
  tibble(
    PlotID = trimws(as.character(raw$PlotID)),
    Compartment = compartment,
    Target = "phoD",
    ASV_richness_calc = rowSums(mat > 0, na.rm = TRUE)
  )
}

pad_matrix <- function(mat, all_cols) {
  missing_cols <- setdiff(all_cols, colnames(mat))
  if (length(missing_cols) > 0) {
    zero_block <- matrix(
      0,
      nrow = nrow(mat),
      ncol = length(missing_cols),
      dimnames = list(rownames(mat), missing_cols)
    )
    mat <- cbind(mat, zero_block)
  }
  mat[, all_cols, drop = FALSE]
}

clr_transform <- function(mat) {
  mat <- as.matrix(mat) + 1
  log_mat <- log(mat)
  sweep(log_mat, 1, rowMeans(log_mat), "-")
}

calc_paired_aitchison <- function(soil_mat, gut_mat) {
  common_pairs <- sort(intersect(rownames(soil_mat), rownames(gut_mat)))
  
  all_asv_names <- union(colnames(soil_mat), colnames(gut_mat))
  soil_pad <- pad_matrix(soil_mat[common_pairs, , drop = FALSE], all_asv_names)
  gut_pad  <- pad_matrix(gut_mat[common_pairs, , drop = FALSE], all_asv_names)
  
  soil_clr <- clr_transform(soil_pad)
  gut_clr  <- clr_transform(gut_pad)
  
  tibble(
    PlotID = normalise_plotid(common_pairs),
    D_16S_Aitchison = sqrt(rowSums((soil_clr - gut_clr)^2))
  )
}

# -------------------------
# 3) Main metrics
main_metrics <- readxl::read_excel(main_xlsx, sheet = main_sheet) %>%
  transmute(
    PlotID     = trimws(as.character(PlotID)),
    Duration   = factor(as.character(Duration), levels = c("5y", "8y", "10y")),
    Regime     = factor(as.character(Regime),   levels = c("CK", "NPK", "NPKOM")),
    gut16S_PC1 = as.numeric(gut16S_PC1),
    D_phoD     = as.numeric(D_phoD),
    D_phoD_New = as.numeric(D_phoD_BC)
  )
# -------------------------
# -------------------------

meta_cols <- c("PlotID", "Duration", "Regime", "Compartment")

soil_raw <- readxl::read_excel(data_xlsx, sheet = "16S_soil_ASV_table") %>%
  as.data.frame()
gut_raw  <- readxl::read_excel(data_xlsx, sheet = "16S_gut_ASV_table") %>%
  as.data.frame()

to_seqtab <- function(df) {
  asv_cols <- setdiff(colnames(df), meta_cols)
  mat <- as.matrix(sapply(df[, asv_cols, drop = FALSE], as.numeric))
  rownames(mat) <- trimws(as.character(df$PlotID))
  mat
}

soil16S <- to_seqtab(soil_raw)
gut16S  <- to_seqtab(gut_raw)

qc_16S <- bind_rows(
  make_16S_qc(soil16S, "Soil"),
  make_16S_qc(gut16S,  "Gut")
)

D_16S <- calc_paired_aitchison(soil16S, gut16S)

# -------------------------
# 5) phoD tracking
qc_phoD <- bind_rows(
  read_phoD_tracking(data_xlsx, "phoD_reads_tracking_soil", "Soil"),
  read_phoD_tracking(data_xlsx, "phoD_reads_tracking_gut",  "Gut")
) %>%
  select(-ASV_richness) %>%
  left_join(
    bind_rows(
      make_phoD_richness(data_xlsx, "phoD_soil_ASV_table", "Soil"),
      make_phoD_richness(data_xlsx, "phoD_gut_ASV_table", "Gut")
    ),
    by = c("PlotID", "Compartment", "Target")
  ) %>%
  rename(ASV_richness = ASV_richness_calc)
# -------------------------
# 6) Sample-level QC table
# -------------------------
sample_qc <- bind_rows(qc_16S, qc_phoD) %>%
  left_join(
    main_metrics %>% dplyr::select(PlotID, gut16S_PC1),
    by = "PlotID"
  ) %>%
  mutate(
    gut16S_PC1 = if_else(
      Compartment == "Gut" & Target == "16S",
      gut16S_PC1,
      NA_real_
    )
  ) %>%
  dplyr::select(
    SampleID,
    PlotID,
    Duration,
    Regime,
    Replicate,
    Compartment,
    Target,
    Raw_reads,
    Filtered_reads,
    Nonchim_reads,
    ASV_richness,
    gut16S_PC1,
    Notes
  ) %>%
  arrange(
    PlotID,
    factor(Target, levels = c("16S", "phoD")),
    factor(Compartment, levels = c("Soil", "Gut"))
  )

# -------------------------
# 7) Pair-level metrics table
# -------------------------
pair_metrics <- main_metrics %>%
  left_join(D_16S, by = "PlotID") %>%
  mutate(
    PairID = PlotID,
    Interpretation = "Lower D_phoD indicates stronger soil-gut phoD community coupling."
  ) %>%
  dplyr::select(
    PlotID,
    PairID,
    Duration,
    Regime,
    gut16S_PC1,
    D_16S_Aitchison,
    D_phoD,
    D_phoD_New,
    Interpretation
  ) %>%
  arrange(Duration, Regime, PlotID)

# -------------------------
# 8) Export Excel
# -------------------------
wb <- createWorkbook()

addWorksheet(wb, "S5_sample_level_QC")
addWorksheet(wb, "S5_pair_level_metrics")
addWorksheet(wb, "README")

writeData(wb, "S5_sample_level_QC", sample_qc)
writeData(wb, "S5_pair_level_metrics", pair_metrics)

readme <- tibble(
  Field = c(
    "Supplementary Table S5",
    "Raw_reads",
    "Filtered_reads",
    "Nonchim_reads",
    "ASV_richness",
    "gut16S_PC1",
    "D_16S_Aitchison",
    "D_phoD",
    "D_phoD_New"
  ),
  Description = c(
    "Amplicon sequencing summary and compositional metrics.",
    "For phoD, imported from DADA2 reads tracking tables. For 16S, left as NA unless 16S tracking tables are provided.",
    "Reads after filtering when available from tracking tables; for 16S, calculated from seqtab row sums.",
    "Non-chimeric reads when available; for 16S, calculated from seqtab row sums.",
    "Number of non-zero ASV/features in the seqtab, calculated from 16S and phoD ASV tables.",
    "Gut microbial state derived from Aitchison PCA of 16S communities.",
    "Paired soil-gut 16S Aitchison distance.",
    "Aitchison-based paired soil-gut phoD community dissimilarity.",
    "Bray-Curtis-based robustness version of paired soil-gut phoD community dissimilarity."
  )
)

writeData(wb, "README", readme)

header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  border = "Bottom"
)

for (sh in names(wb)) {
  addStyle(
    wb,
    sheet = sh,
    style = header_style,
    rows = 1,
    cols = 1:50,
    gridExpand = TRUE
  )
  setColWidths(wb, sheet = sh, cols = 1:50, widths = "auto")
  freezePane(wb, sheet = sh, firstRow = TRUE)
}

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

cat("\nSaved Supplementary Table S5 to:\n")
cat(out_xlsx, "\n")

cat("\nCheck missing values:\n")
cat("16S Raw_reads missing:", sum(is.na(sample_qc$Raw_reads[sample_qc$Target == "16S"])), "\n")
cat("phoD Raw_reads missing:", sum(is.na(sample_qc$Raw_reads[sample_qc$Target == "phoD"])), "\n")
cat("D_phoD missing:", sum(is.na(pair_metrics$D_phoD)), "\n")
cat("D_phoD_New missing:", sum(is.na(pair_metrics$D_phoD_New)), "\n")
cat("D_16S_Aitchison missing:", sum(is.na(pair_metrics$D_16S_Aitchison)), "\n")
