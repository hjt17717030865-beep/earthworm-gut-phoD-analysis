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
# Supplementary Table 10 | Robustness analyses of the
# phoD-independent P-index
#
# Purpose:
#   Generate NC-ready robustness tables for the phoD-independent
#   P-index, corresponding to Fig. S7.
#
# Analyses:
#   AI ~ P_index across:
#     1) Main 6-gene set
#     2) Without gcd
#     3) Without pqqC
#   and across three scaling methods:
#     1) Z-score
#     2) log10(x + 1) + Z-score
#     3) Min-max normalisation
#
# Critical rule:
#   phoD is explicitly excluded from all P-index versions to avoid
#   circularity with phoD amplification.
#
# Output:
#   Supplementary_Table_10_Robustness_phoD_independent_Pindex.xlsx
#   Supplementary_Table_10_Robustness_phoD_independent_Pindex.docx
#
# Input:
#
# Required P-gene table structure:
#   - first column named gene
#   - remaining columns are sample IDs
#
# Required main table columns:
#   PlotID, AI
# ============================================================

# -------------------------
# 0) Packages
# -------------------------
pkg_needed <- c(
  "readxl", "dplyr", "tidyr", "stringr", "tibble",
  "purrr", "broom", "openxlsx", "officer", "flextable"
)

pkg_to_install <- pkg_needed[!sapply(pkg_needed, requireNamespace, quietly = TRUE)]
if (length(pkg_to_install) > 0) {
  install.packages(pkg_to_install)
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(broom)
  library(openxlsx)
  library(officer)
  library(flextable)
})
select <- dplyr::select
filter <- dplyr::filter

# -------------------------
# 1) Paths
# -------------------------
# -------------------------
# 1) Paths
# -------------------------
out_dir <- OUTPUT_DIR
file_suppl <- DATA_XLSX

pgene_sheet <- "P_gene_qPCR_matrix"
main_sheet  <- "Master_all_data"

out_xlsx <- file.path(
  out_dir,
  "Supplementary_Table_10_Robustness_phoD_independent_Pindex.xlsx"
)

out_docx <- file.path(
  out_dir,
  "Supplementary_Table_10_Robustness_phoD_independent_Pindex.docx"
)

if (!file.exists(file_suppl)) {
  stop("Supplementary Data file not found: ", file_suppl)
}

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# -------------------------
# 2) Helper functions
# -------------------------
fmt_p <- function(p) {
  ifelse(
    is.na(p), NA_character_,
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), NA_character_, sprintf(paste0("%.", digits, "f"), x))
}

safe_z <- function(x) {
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))

  if (sum(ok) >= 2 && sd(x[ok], na.rm = TRUE) > 0) {
    out[ok] <- as.numeric(scale(x[ok]))
  }

  out
}

safe_log1p_z <- function(x) {
  safe_z(log10(x + 1))
}

safe_minmax <- function(x) {
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))

  if (sum(ok) >= 2) {
    rng <- range(x[ok], na.rm = TRUE)
    if ((rng[2] - rng[1]) > 0) {
      out[ok] <- (x[ok] - rng[1]) / (rng[2] - rng[1])
    }
  }

  out
}

clean_duration <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\s+", "", x)
  x <- gsub("^3y$", "5y", x)
  x <- gsub("^3$", "5y", x)
  x <- gsub("^5$", "5y", x)
  x <- gsub("^8$", "8y", x)
  x <- gsub("^10$", "10y", x)
  factor(x, levels = c("5y", "8y", "10y"))
}

clean_gene_name <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("\\s+", "", x)
  x <- gsub("-", "", x)
  x <- ifelse(x == "ppk3", "ppk", x)
  x
}

gene_set_label <- function(x) {
  dplyr::case_when(
    x == "main_6"  ~ "Main 6-gene set",
    x == "no_gcd"  ~ "Without gcd",
    x == "no_pqqc" ~ "Without pqqC",
    TRUE ~ x
  )
}

scaling_label <- function(x) {
  dplyr::case_when(
    x == "z"       ~ "Z-score",
    x == "log1p_z" ~ "log10(x+1) + Z-score",
    x == "minmax"  ~ "Min-max",
    TRUE ~ x
  )
}

fit_ai_pindex <- function(df) {
  fit <- lm(AI ~ P_index, data = df)
  sm <- summary(fit)

  broom::tidy(fit, conf.int = TRUE) %>%
    filter(term == "P_index") %>%
    transmute(
      N = nobs(fit),
      Estimate = estimate,
      SE = std.error,
      CI_low = conf.low,
      CI_high = conf.high,
      t = statistic,
      P_value = p.value,
      P = fmt_p(p.value),
      R2 = sm$r.squared,
      Adj_R2 = sm$adj.r.squared,
      Df_residual = df.residual(fit)
    )
}

# -------------------------
# 3) Read P-gene table and convert to long format
# -------------------------
dat_raw <- readxl::read_excel(
  file_suppl,
  sheet = pgene_sheet
)

if (!"gene" %in% names(dat_raw)) {
  stop("Input P-gene table must contain a column named 'gene'.")
}

dat_long <- dat_raw %>%
  pivot_longer(
    cols = -gene,
    names_to = "SampleID",
    values_to = "copies"
  ) %>%
  mutate(
    gene = clean_gene_name(gene),
    SampleID = as.character(SampleID),
    copies = suppressWarnings(as.numeric(copies)),
    Duration = str_extract(SampleID, "^\\d+y"),
    Duration = clean_duration(Duration),

    Habitat = case_when(
      str_detect(SampleID, "^\\d+yEB") ~ "gut",
      str_detect(SampleID, "^\\d+yE")  ~ "gut",
      TRUE ~ "soil"
    ),

    Taxon = case_when(
      str_detect(SampleID, "^\\d+yEB") ~ "EB",
      str_detect(SampleID, "^\\d+yE")  ~ "E",
      TRUE ~ "soil"
    ),

    Regime = case_when(
      str_detect(SampleID, "NPKOM") ~ "NPKOM",
      str_detect(SampleID, "NPK")   ~ "NPK",
      str_detect(SampleID, "CK")    ~ "CK",
      TRUE ~ NA_character_
    ),

    Rep = suppressWarnings(as.numeric(str_extract(SampleID, "(?<=-)\\d+$"))),

    PlotID = SampleID %>%
      str_replace("^([0-9]+y)EB", "\\1") %>%
      str_replace("^([0-9]+y)E", "\\1") %>%
      str_replace("^3y", "5y")
  ) %>%
  dplyr::select(SampleID, PlotID, Duration, Habitat, Taxon, Regime, Rep, gene, copies)
# -------------------------
# 4) Define phoD-independent P-index versions
# -------------------------
gene_sets <- list(
  main_6  = c("phnk", "phox", "gcd", "pqqc", "ppk", "ppx"),
  no_gcd  = c("phnk", "phox", "pqqc", "ppk", "ppx"),
  no_pqqc = c("phnk", "phox", "gcd", "ppk", "ppx")
)

scalings <- c("z", "log1p_z", "minmax")

# Critical circularity check
all_genes_requested <- unique(unlist(gene_sets))
if ("phod" %in% all_genes_requested) {
  stop("phoD must not be included in any P-index gene set.")
}

calc_pindex_general <- function(df_long,
                                genes_use,
                                gene_set_name,
                                scaling = c("z", "log1p_z", "minmax")) {

  scaling <- match.arg(scaling)

  genes_use <- clean_gene_name(genes_use)

  if ("phod" %in% genes_use) {
    stop("phoD is included in gene set ", gene_set_name, ". Remove phoD to avoid circularity.")
  }

  df_use <- df_long %>%
    filter(gene %in% genes_use)

  if (nrow(df_use) == 0) {
    stop("No rows found for gene set: ", gene_set_name)
  }

  missing_genes <- setdiff(genes_use, unique(df_use$gene))

  if (length(missing_genes) > 0) {
    warning(
      "Missing genes for ", gene_set_name, ": ",
      paste(missing_genes, collapse = ", ")
    )
  }

  df_scaled <- df_use %>%
    group_by(gene) %>%
    mutate(
      scaled_value = case_when(
        scaling == "z"       ~ safe_z(copies),
        scaling == "log1p_z" ~ safe_log1p_z(copies),
        scaling == "minmax"  ~ safe_minmax(copies)
      )
    ) %>%
    ungroup()

  df_sample <- df_scaled %>%
    filter(
      Habitat == "gut",
      str_detect(SampleID, "^\\d+yE"),
      !str_detect(SampleID, "^\\d+yEB")
    ) %>%
    group_by(SampleID, PlotID, Duration, Habitat, Taxon, Regime, Rep) %>%
    summarise(
      P_index = mean(scaled_value, na.rm = TRUE),
      n_selected_genes = length(genes_use),
      n_detected_genes = n_distinct(gene),
      n_valid_gene_values = sum(is.finite(scaled_value)),
      usable_for_Pindex = n_valid_gene_values > 0,
      .groups = "drop"
    ) %>%
    mutate(
      gene_set_name = gene_set_name,
      scaling = scaling,
      included_genes = paste(genes_use, collapse = ", "),
      missing_genes = ifelse(
        length(missing_genes) == 0,
        "None",
        paste(missing_genes, collapse = ", ")
      ),
      phoD_status = "Excluded"
    ) %>%
    arrange(Duration, Regime, Rep, SampleID)

  df_sample
}

pindex_robust <- bind_rows(
  lapply(names(gene_sets), function(gs) {
    bind_rows(
      lapply(scalings, function(sc) {
        calc_pindex_general(
          df_long = dat_long,
          genes_use = gene_sets[[gs]],
          gene_set_name = gs,
          scaling = sc
        )
      })
    )
  })
)

# -------------------------
# 5) Read main table and merge AI
# -------------------------
dat_main <- readxl::read_excel(
  file_suppl,
  sheet = main_sheet
) %>%
  mutate(
    PlotID = as.character(PlotID),
    PlotID = gsub("^3y", "5y", PlotID),
    AI = as.numeric(AI)
  ) %>%
  select(PlotID, AI)


if (anyDuplicated(dat_main$PlotID) > 0) {
  stop("Duplicated PlotID found in main data table.")
}

pindex_ai <- pindex_robust %>%
  left_join(dat_main, by = "PlotID")

if (sum(is.na(pindex_ai$AI)) > 0) {
  missing_ids <- pindex_ai %>%
    filter(is.na(AI)) %>%
    distinct(PlotID) %>%
    pull(PlotID)

  stop(
    "Some AI values are missing after merging. Check PlotID matching: ",
    paste(missing_ids, collapse = ", ")
  )
}

# -------------------------
# 6) Fit AI ~ P_index for each robustness version
# -------------------------
s10_model_results <- pindex_ai %>%
  group_by(gene_set_name, scaling) %>%
  nest() %>%
  mutate(
    model_table = map(data, fit_ai_pindex)
  ) %>%
  select(-data) %>%
  unnest(model_table) %>%
  ungroup() %>%
  mutate(
    Gene_set = gene_set_label(gene_set_name),
    Scaling = scaling_label(scaling),
    Response = "AI",
    Predictor = "phoD-independent P-index",
    Formula = "AI ~ P_index"
  ) %>%
  arrange(
    factor(gene_set_name, levels = c("main_6", "no_gcd", "no_pqqc")),
    factor(scaling, levels = c("z", "log1p_z", "minmax"))
  )

# -------------------------
# 7) P-index version summary
# -------------------------
s10_version_summary <- pindex_robust %>%
  group_by(gene_set_name, scaling, included_genes, missing_genes, phoD_status) %>%
  summarise(
    N_samples = n_distinct(SampleID),
    N_pairs = n_distinct(PlotID),
    Mean_valid_gene_values = mean(n_valid_gene_values, na.rm = TRUE),
    Min_valid_gene_values = min(n_valid_gene_values, na.rm = TRUE),
    Max_valid_gene_values = max(n_valid_gene_values, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Gene_set = gene_set_label(gene_set_name),
    Scaling = scaling_label(scaling),
    Excluded_gene_in_sensitivity = case_when(
      gene_set_name == "main_6"  ~ "None",
      gene_set_name == "no_gcd"  ~ "gcd",
      gene_set_name == "no_pqqc" ~ "pqqC",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(
    factor(gene_set_name, levels = c("main_6", "no_gcd", "no_pqqc")),
    factor(scaling, levels = c("z", "log1p_z", "minmax"))
  )

# -------------------------
# 8) Word-ready compact tables
# -------------------------
s10a_versions_Word <- s10_version_summary %>%
  transmute(
    Gene_set,
    Scaling,
    Included_genes = included_genes,
    Sensitivity_exclusion = Excluded_gene_in_sensitivity,
    phoD_status,
    N_pairs,
    N_samples,
    Valid_gene_values = paste0(
      fmt_num(Mean_valid_gene_values, 1),
      " [",
      fmt_num(Min_valid_gene_values, 0),
      "-",
      fmt_num(Max_valid_gene_values, 0),
      "]"
    )
  )

s10b_models_Word <- s10_model_results %>%
  transmute(
    Gene_set,
    Scaling,
    Response,
    Predictor,
    Estimate = round(Estimate, 3),
    SE = round(SE, 3),
    CI = paste0("[", fmt_num(CI_low, 3), ", ", fmt_num(CI_high, 3), "]"),
    t = round(t, 3),
    P,
    R2 = round(R2, 3),
    N
  )

# -------------------------
# 9) Caption / README
# -------------------------
caption <- tibble(
  Caption = c(
    "Supplementary Table 10 | Robustness analyses of the phoD-independent P-index.",
    "Panel a summarises phoD-independent P-index versions constructed using alternative gene subsets and scaling methods.",
    "Panel b reports associations between phoD amplification (AI) and each P-index version.",
    "The main P-index used six P-acquisition genes excluding phoD: phnK, phoX, gcd, pqqC, ppk and ppx.",
    "Sensitivity analyses excluded gcd or pqqC to evaluate the influence of P-solubilisation genes.",
    "All P-index versions explicitly excluded phoD to avoid circularity with phoD amplification.",
    "Models are interpreted as descriptive and associational."
  )
)

readme <- tibble(
  Field = c(
    "Table title",
    "Panel a",
    "Panel b",
    "Main gene set",
    "Sensitivity gene sets",
    "Scaling methods",
    "Critical circularity control",
    "Input P-gene file",
    "Input main table"
  ),
  Description = c(
    "Supplementary Table 10 | Robustness analyses of the phoD-independent P-index",
    "Definitions and sample coverage of P-index robustness versions.",
    "Linear models testing AI ~ P_index across gene subsets and scaling methods.",
    "phnK, phoX, gcd, pqqC, ppk and ppx.",
    "Without gcd; without pqqC.",
    "Z-score; log10(x+1) + Z-score; min-max normalisation.",
    "phoD is excluded from all P-index versions.",
    paste0(file_suppl, " | sheet: ", pgene_sheet),
    paste0(file_suppl, " | sheet: ", main_sheet)
  )
)

# -------------------------
# 10) Write Excel workbook
# -------------------------
wb <- createWorkbook()

addWorksheet(wb, "README")
writeData(wb, "README", readme)

addWorksheet(wb, "Caption")
writeData(wb, "Caption", caption)

addWorksheet(wb, "S10a_versions_Word")
writeData(wb, "S10a_versions_Word", s10a_versions_Word)

addWorksheet(wb, "S10b_models_Word")
writeData(wb, "S10b_models_Word", s10b_models_Word)

addWorksheet(wb, "S10_full_model_results")
writeData(wb, "S10_full_model_results", s10_model_results)

addWorksheet(wb, "S10_full_version_summary")
writeData(wb, "S10_full_version_summary", s10_version_summary)

addWorksheet(wb, "S10_sample_level_Pindex")
writeData(wb, "S10_sample_level_Pindex", pindex_ai)

for (sh in names(wb)) {
  freezePane(wb, sh, firstRow = TRUE)
  setColWidths(wb, sh, cols = 1:80, widths = "auto")
}

num_style <- createStyle(numFmt = "0.000")

sheet_objects <- list(
  S10a_versions_Word = s10a_versions_Word,
  S10b_models_Word = s10b_models_Word,
  S10_full_model_results = s10_model_results,
  S10_full_version_summary = s10_version_summary,
  S10_sample_level_Pindex = pindex_ai
)

for (sh in names(sheet_objects)) {
  dat_sh <- sheet_objects[[sh]]
  numeric_cols <- which(vapply(dat_sh, is.numeric, logical(1)))
  if (length(numeric_cols) > 0 && nrow(dat_sh) > 0) {
    addStyle(
      wb, sh, style = num_style,
      rows = 2:(nrow(dat_sh) + 1),
      cols = numeric_cols,
      gridExpand = TRUE,
      stack = TRUE
    )
  }
}

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

# -------------------------
# 11) Write Word-ready document
# -------------------------
make_ft <- function(df, font_size = 7.5) {
  flextable(df) %>%
    theme_booktabs() %>%
    fontsize(size = font_size, part = "all") %>%
    fontsize(size = font_size + 0.5, part = "header") %>%
    align(align = "center", part = "all") %>%
    align(
      j = intersect(c("Gene_set", "Scaling", "Included_genes", "Sensitivity_exclusion",
                      "phoD_status", "Response", "Predictor", "CI"), names(df)),
      align = "left",
      part = "body"
    ) %>%
    autofit()
}

doc <- read_docx()

doc <- body_add_par(
  doc,
  "Supplementary Table 10 | Robustness analyses of the phoD-independent P-index.",
  style = "heading 1"
)

doc <- body_add_par(
  doc,
  paste(
    "Panel a summarises phoD-independent P-index versions constructed using alternative gene subsets and scaling methods.",
    "Panel b reports associations between phoD amplification (AI) and each P-index version.",
    "The main P-index used six P-acquisition genes excluding phoD: phnK, phoX, gcd, pqqC, ppk and ppx.",
    "Sensitivity analyses excluded gcd or pqqC to evaluate the influence of P-solubilisation genes.",
    "All P-index versions explicitly excluded phoD to avoid circularity with phoD amplification.",
    "Models are interpreted as descriptive and associational."
  ),
  style = "Normal"
)

doc <- body_add_par(doc, "a, phoD-independent P-index versions.", style = "heading 2")
doc <- body_add_flextable(doc, make_ft(s10a_versions_Word, font_size = 7.0))

doc <- body_add_par(doc, "b, Associations between AI and phoD-independent P-index versions.", style = "heading 2")
doc <- body_add_flextable(doc, make_ft(s10b_models_Word, font_size = 7.5))

print(doc, target = out_docx)

cat("\nSupplementary Table 10 generated successfully:\n")
cat(out_xlsx, "\n")
cat(out_docx, "\n\n")

cat("Word-ready sheets:\n")
cat("- S10a_versions_Word\n")
cat("- S10b_models_Word\n\n")

cat("Model rows:", nrow(s10b_models_Word), "\n")
