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
# Supplementary Table 11 / Supplementary Data 1 |
# Genome-resolved P-acquisition repertoire in recovered MAGs
#
# This complete version:
#   1) uses the same input files as your Figure S9/S11 MAG repertoire script;
#   2) merges the sample-information mapping table to create readable MAG_IDs;
#   3) merges soil MAG taxonomy from GTDB-Tk ar53/bac120 summaries;
#   4) outputs Word-ready Supplementary Table 11 and machine-readable
#      Supplementary Data 1.
#
# Outputs:
#   Supplementary_Table_11_MAG_P_acquisition_repertoire.docx
#   Supplementary_Data_1_MAG_P_acquisition_repertoire.xlsx
#
# Interpretation:
#   Metagenome/MAG results are functional context only, not causal evidence.
# ============================================================

# -------------------------
# 0) Packages
# -------------------------
pkg_needed <- c(
  "tidyverse", "readr", "readxl", "openxlsx", "officer", "flextable"
)

pkg_to_install <- pkg_needed[!sapply(pkg_needed, requireNamespace, quietly = TRUE)]
if (length(pkg_to_install) > 0) install.packages(pkg_to_install)

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(readxl)
  library(openxlsx)
  library(officer)
  library(flextable)
})
select <- dplyr::select
# -------------------------
# -------------------------
# 1) Paths
# -------------------------
out_dir <- OUTPUT_DIR
supp_xlsx <- DATA_XLSX
out_xlsx <- file.path(out_dir, "Supplementary_Data_1_MAG_P_acquisition_repertoire.xlsx")
out_docx <- file.path(out_dir, "Supplementary_Table_11_MAG_P_acquisition_repertoire.docx")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!file.exists(supp_xlsx)) stop("Supplementary Data workbook not found: ", supp_xlsx)

required_sheets <- c(
  "P_MAG_summary_tax_compartment",
  "P_MAG_gene_presence",
  "MAG_sample_map",
  "GTDB_taxonomy_summary"
)

available_sheets <- readxl::excel_sheets(supp_xlsx)
missing_sheets <- setdiff(required_sheets, available_sheets)

if (length(missing_sheets) > 0) {
  stop("Supplementary Data.xlsx missing sheets: ", paste(missing_sheets, collapse = ", "))
}

# -------------------------
# 2) Helper functions
# -------------------------
fmt_num <- function(x, digits = 1) {
  ifelse(is.na(x), NA_character_, sprintf(paste0("%.", digits, "f"), x))
}

pick_optional_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

safe_percent <- function(n, denom) {
  if (denom == 0 || is.na(denom)) return(NA_real_)
  100 * n / denom
}

clean_mag_name <- function(x) {
  x <- as.character(x)
  x <- basename(x)
  x <- str_remove(x, "\\.fa$")
  x <- str_remove(x, "\\.fna$")
  x <- str_remove(x, "\\.fasta$")
  x <- str_remove(x, "\\.faa$")
  x <- str_remove(x, "\\.tsv$")
  x <- str_remove(x, "\\.txt$")
  x <- trimws(x)
  x
}

clean_duration <- function(x) {
  x <- as.character(x)
  case_when(
    str_detect(x, "10") ~ "10y",
    str_detect(x, "8")  ~ "8y",
    str_detect(x, "5|3") ~ "5y",
    TRUE ~ x
  )
}

presence_label <- function(df, cols, labels = NULL) {
  if (is.null(labels)) labels <- setNames(cols, cols)
  
  m <- df[, cols, drop = FALSE]
  
  vapply(seq_len(nrow(m)), function(i) {
    x <- suppressWarnings(as.integer(unlist(m[i, cols, drop = TRUE])))
    x[is.na(x)] <- 0L
    present <- cols[x == 1L]
    
    if (length(present) == 0) {
      "None"
    } else {
      paste(unname(labels[present]), collapse = "; ")
    }
  }, character(1))
}

standardise_gtdb_taxonomy <- function(x) {
  x <- as.character(x)
  x <- ifelse(is.na(x) | x == "" | x == "NA", NA_character_, x)
  x
}

# -------------------------
# 3) Read input tables
# -------------------------
mag_detail_raw <- readxl::read_excel(
  supp_xlsx,
  sheet = "P_MAG_summary_tax_compartment"
)

mag_gene_raw <- readxl::read_excel(
  supp_xlsx,
  sheet = "P_MAG_gene_presence"
)

need_meta <- c("MAG", "compartment", "KO_list")
miss_meta <- setdiff(need_meta, names(mag_detail_raw))
if (length(miss_meta) > 0) {
  stop("P_MAG_summary_tax_compartment.tsv missing columns: ", paste(miss_meta, collapse = ", "))
}

gene_sets <- c("phoD", "phoX", "pqqC", "ppk", "ppx", "phnK")
need_gene <- c("MAG", "compartment", gene_sets)
miss_gene <- setdiff(need_gene, names(mag_gene_raw))
if (length(miss_gene) > 0) {
  stop("P_MAG_gene_presence.tsv missing columns: ", paste(miss_gene, collapse = ", "))
}

# -------------------------
# 4) Read sample-information mapping table
# -------------------------
sample_map_raw <- readxl::read_excel(
  supp_xlsx,
  sheet = "MAG_sample_map"
)
if (ncol(sample_map_raw) < 9) {
  stop("Sample map must contain at least 9 columns.")
}

names(sample_map_raw)[1:9] <- c(
  "MAG_file",
  "Sample_code",
  "Sample_object",
  "Duration_raw",
  "Regime",
  "Replicate_info",
  "Total_replicates",
  "Merged",
  "Merge_source"
)

sample_map <- sample_map_raw %>%
  mutate(
    MAG = clean_mag_name(MAG_file),
    Sample_code = as.character(Sample_code),
    
    Compartment_map = case_when(
      str_detect(str_to_lower(as.character(Sample_object)), "^soil$") ~ "Soil",
      str_detect(str_to_lower(as.character(Sample_object)), "^gut$|foetida|earthworm") ~ "Gut_Efoetida",
      TRUE ~ as.character(Sample_object)
    ),
    
    Duration_map = clean_duration(Duration_raw),
    Regime_map = as.character(Regime),
    Replicate_info = as.character(Replicate_info),
    Merged = as.character(Merged),
    Merge_source = as.character(Merge_source),
    
    Bin_ID = str_extract(MAG, "bin\\.?[0-9]+$"),
    Bin_ID = str_replace(Bin_ID, "\\.", ""),
    Bin_ID = ifelse(is.na(Bin_ID), MAG, Bin_ID),
    
    MAG_ID = paste(
      Compartment_map,
      Duration_map,
      Regime_map,
      Sample_code,
      Bin_ID,
      sep = "_"
    )
  ) %>%
  dplyr::select(
    MAG,
    MAG_ID,
    Sample_code,
    Compartment_map,
    Duration_map,
    Regime_map,
    Replicate_info,
    Total_replicates,
    Merged,
    Merge_source
  ) %>%
  distinct(MAG, .keep_all = TRUE)

# -------------------------
# 5) Read GTDB taxonomy from Supplementary Data.xlsx
# -------------------------
gtdb_tax_raw <- readxl::read_excel(
  supp_xlsx,
  sheet = "GTDB_taxonomy_summary"
)

if (!all(c("user_genome", "classification") %in% names(gtdb_tax_raw))) {
  stop("GTDB_taxonomy_summary must contain user_genome and classification columns.")
}

gtdb_all <- gtdb_tax_raw %>%
  transmute(
    MAG = clean_mag_name(user_genome),
    GTDB_taxonomy = standardise_gtdb_taxonomy(classification)
  ) %>%
  dplyr::distinct(MAG, .keep_all = TRUE)

# -------------------------
# 6) Prepare MAG module and gene tables
# -------------------------
mag_detail <- mag_detail_raw %>%
  mutate(
    MAG = clean_mag_name(MAG),
    compartment = str_to_lower(as.character(compartment)),
    KO_list = replace_na(as.character(KO_list), "")
  ) %>%
  filter(compartment %in% c("soil", "gut"))

mag_gene <- mag_gene_raw %>%
  mutate(
    MAG = clean_mag_name(MAG),
    compartment = str_to_lower(as.character(compartment)),
    across(all_of(gene_sets), ~ suppressWarnings(as.integer(.x)))
  ) %>%
  filter(compartment %in% c("soil", "gut"))

if (!"P_score" %in% names(mag_gene)) {
  mag_gene <- mag_gene %>%
    mutate(P_score = rowSums(across(all_of(gene_sets)), na.rm = TRUE))
} else {
  mag_gene <- mag_gene %>%
    mutate(P_score = suppressWarnings(as.integer(P_score)))
}

# Collapse possible duplicate MAG rows.
mag_gene <- mag_gene %>%
  group_by(MAG, compartment) %>%
  summarise(
    across(all_of(gene_sets), ~ max(.x, na.rm = TRUE)),
    P_score = max(P_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    across(all_of(gene_sets), ~ ifelse(is.infinite(.x), 0L, as.integer(.x))),
    P_score = ifelse(is.infinite(P_score), 0L, as.integer(P_score))
  )

# -------------------------
# 7) Module presence from KO_list
#    These definitions reproduce the Figure S9/S11 script.
# -------------------------
module_cols <- c(
  "P_transport",
  "PolyP_turnover",
  "Organic_P_mineralization",
  "Phosphonate",
  "Inorganic_P_solubilization"
)

module_labels <- c(
  P_transport = "P transport",
  PolyP_turnover = "PolyP turnover",
  Organic_P_mineralization = "Organic P mineralization",
  Phosphonate = "Phosphonate",
  Inorganic_P_solubilization = "Inorganic P solubilization"
)

mag_module <- mag_detail %>%
  transmute(
    MAG,
    compartment,
    
    # Phosphate transport / uptake
    P_transport = as.integer(str_detect(KO_list, "K02036|K02037|K02038|K02040|K02033")),
    
    # Polyphosphate turnover
    PolyP_turnover = as.integer(str_detect(KO_list, "K00937|K01507")),
    
    # Organic P mineralization
    Organic_P_mineralization = as.integer(str_detect(KO_list, "K01113|K01093|K01077")),
    
    # Phosphonate metabolism
    Phosphonate = as.integer(str_detect(KO_list, "K06164|K06165|K05780")),
    
    # Inorganic P solubilization / PQQ-related
    Inorganic_P_solubilization = as.integer(str_detect(KO_list, "K00117|K06130"))
  ) %>%
  group_by(MAG, compartment) %>%
  summarise(
    across(all_of(module_cols), ~ max(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(across(all_of(module_cols), ~ ifelse(is.infinite(.x), 0L, as.integer(.x))))

# -------------------------
# 8) Optional MAG metadata columns from file_meta
# -------------------------
tax_col <- pick_optional_col(
  mag_detail_raw,
  c("taxonomy", "Taxonomy", "GTDB_taxonomy", "gtdb_taxonomy", "classification", "Classification")
)

sample_col <- pick_optional_col(
  mag_detail_raw,
  c("SampleID", "sample_id", "Sample", "sample", "PlotID", "plot_id")
)

mag_meta <- mag_detail_raw %>%
  mutate(
    MAG = clean_mag_name(MAG),
    SampleID = if (!is.na(sample_col)) as.character(.data[[sample_col]]) else NA_character_,
    Taxonomy_raw = if (!is.na(tax_col)) as.character(.data[[tax_col]]) else NA_character_
  ) %>%
  select(MAG, SampleID, Taxonomy_raw) %>%
  distinct(MAG, .keep_all = TRUE) %>%
  left_join(sample_map, by = "MAG") %>%
  left_join(gtdb_all, by = "MAG") %>%
  mutate(
    Taxonomy = coalesce(GTDB_taxonomy, Taxonomy_raw)
  ) %>%
  dplyr::select(
    MAG,
    MAG_ID,
    Sample_code,
    Compartment_map,
    Duration_map,
    Regime_map,
    Replicate_info,
    Total_replicates,
    Merged,
    Merge_source,
    SampleID,
    Taxonomy
  )

# -------------------------
# 9) MAG-level repertoire summary
# -------------------------
gene_module_map <- tibble(
  Gene = gene_sets,
  Functional_module = c(
    "Organic P mineralization",
    "Organic P mineralization",
    "Inorganic P solubilization",
    "PolyP turnover",
    "PolyP turnover",
    "Phosphonate"
  ),
  Functional_interpretation = c(
    "Alkaline phosphatase D",
    "Alkaline phosphatase X",
    "Pyrroloquinoline quinone biosynthesis protein C",
    "Polyphosphate kinase",
    "Exopolyphosphatase",
    "C-P lyase complex component"
  )
)

mag_all <- full_join(
  mag_gene,
  mag_module,
  by = c("MAG", "compartment")
) %>%
  mutate(
    across(all_of(gene_sets), ~ replace_na(as.integer(.x), 0L)),
    across(all_of(module_cols), ~ replace_na(as.integer(.x), 0L)),
    P_score = ifelse(is.na(P_score), rowSums(across(all_of(gene_sets)), na.rm = TRUE), P_score)
  ) %>%
  left_join(mag_meta, by = "MAG") %>%
  mutate(
    P_genes = presence_label(., gene_sets),
    P_modules = presence_label(., module_cols, module_labels),
    P_module_count = rowSums(across(all_of(module_cols)), na.rm = TRUE),
    P_gene_count = rowSums(across(all_of(gene_sets)), na.rm = TRUE),
    
    Compartment_final = coalesce(Compartment_map, str_to_title(compartment)),
    MAG_ID_final = coalesce(MAG_ID, paste(Compartment_final, MAG, sep = "_")),
    Duration_final = Duration_map,
    Regime_final = Regime_map
  ) %>%
  arrange(desc(P_score), desc(P_module_count), MAG)

p_carrying_mags <- mag_all %>%
  filter(P_score > 0 | P_module_count > 0)

n_p_mags <- nrow(p_carrying_mags)

# -------------------------
# 10) Module and gene detection summaries
# -------------------------
module_detection_summary <- p_carrying_mags %>%
  summarise(across(all_of(module_cols), ~ sum(.x == 1, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Module_key", values_to = "N_MAGs") %>%
  mutate(
    Functional_module = unname(module_labels[Module_key]),
    Percent_P_carrying_MAGs = map_dbl(N_MAGs, ~ safe_percent(.x, n_p_mags))
  ) %>%
  select(Functional_module, N_MAGs, Percent_P_carrying_MAGs) %>%
  arrange(desc(N_MAGs), Functional_module)

gene_detection_summary <- p_carrying_mags %>%
  summarise(across(all_of(gene_sets), ~ sum(.x == 1, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Gene", values_to = "N_MAGs") %>%
  left_join(gene_module_map, by = "Gene") %>%
  mutate(
    Percent_P_carrying_MAGs = map_dbl(N_MAGs, ~ safe_percent(.x, n_p_mags))
  ) %>%
  select(Gene, Functional_module, Functional_interpretation, N_MAGs, Percent_P_carrying_MAGs) %>%
  arrange(Functional_module, desc(N_MAGs), Gene)

module_cooccurrence <- p_carrying_mags %>%
  mutate(Module_combination = ifelse(P_modules == "None", "None", P_modules)) %>%
  count(Module_combination, name = "N_MAGs") %>%
  arrange(desc(N_MAGs), Module_combination)

# -------------------------
# 11) Machine-readable matrices
# -------------------------
gene_matrix <- mag_all %>%
  select(
    MAG_ID = MAG_ID_final,
    Original_MAG = MAG,
    Sample_code,
    Compartment = Compartment_final,
    Duration = Duration_final,
    Regime = Regime_final,
    Replicate = Replicate_info,
    Merged,
    Merge_source,
    Taxonomy,
    P_score,
    P_gene_count,
    P_module_count,
    all_of(gene_sets)
  ) %>%
  arrange(
    factor(Compartment, levels = c("Soil", "Gut_Efoetida")),
    factor(Duration, levels = c("5y", "8y", "10y")),
    factor(Regime, levels = c("CK", "NPK", "NPKOM")),
    Sample_code,
    desc(P_score),
    Original_MAG
  )

module_matrix <- mag_all %>%
  select(
    MAG_ID = MAG_ID_final,
    Original_MAG = MAG,
    Sample_code,
    Compartment = Compartment_final,
    Duration = Duration_final,
    Regime = Regime_final,
    Replicate = Replicate_info,
    Merged,
    Merge_source,
    Taxonomy,
    P_score,
    P_gene_count,
    P_module_count,
    all_of(module_cols)
  ) %>%
  arrange(
    factor(Compartment, levels = c("Soil", "Gut_Efoetida")),
    factor(Duration, levels = c("5y", "8y", "10y")),
    factor(Regime, levels = c("CK", "NPK", "NPKOM")),
    Sample_code,
    desc(P_score),
    Original_MAG
  )

long_gene_presence <- gene_matrix %>%
  pivot_longer(
    cols = all_of(gene_sets),
    names_to = "Gene",
    values_to = "Present"
  ) %>%
  left_join(gene_module_map, by = "Gene")

long_module_presence <- module_matrix %>%
  pivot_longer(
    cols = all_of(module_cols),
    names_to = "Module_key",
    values_to = "Present"
  ) %>%
  mutate(Functional_module = unname(module_labels[Module_key])) %>%
  select(-Module_key)

# -------------------------
# 12) Word-ready compact tables
# -------------------------
s11a_MAG_summary_Word <- mag_all %>%
  filter(P_score > 0 | P_module_count > 0) %>%
  mutate(
    Compartment_short = case_when(
      Compartment_final == "Gut_Efoetida" ~ "Gut",
      TRUE ~ Compartment_final
    ),
    Source = paste0(
      Compartment_short,
      "; ",
      Duration_final,
      " ",
      Regime_final,
      "; sample ",
      Sample_code,
      "; rep ",
      Replicate_info,
      ifelse(Merged == "Yes", "; merged", "")
    )
  ) %>%
  arrange(
    factor(Compartment_final, levels = c("Soil", "Gut_Efoetida")),
    factor(Duration_final, levels = c("5y", "8y", "10y")),
    factor(Regime_final, levels = c("CK", "NPK", "NPKOM")),
    Sample_code,
    desc(P_score),
    MAG
  ) %>%
  transmute(
    MAG_ID = MAG_ID_final,
    Original_MAG = MAG,
    Source,
    Taxonomy,
    P_score,
    P_modules,
    P_genes
  )
s11b_module_summary_Word <- module_detection_summary %>%
  transmute(
    Functional_module,
    N_MAGs,
    Percent_MAGs = paste0(fmt_num(Percent_P_carrying_MAGs, 1), "%")
  )

s11c_gene_summary_Word <- gene_detection_summary %>%
  transmute(
    Gene,
    Functional_module,
    Functional_interpretation,
    N_MAGs,
    Percent_MAGs = paste0(fmt_num(Percent_P_carrying_MAGs, 1), "%")
  )

# -------------------------
# 13) Caption and README
# -------------------------
caption <- tibble(
  Caption = c(
    "Supplementary Table 11 / Supplementary Data 1 | Genome-resolved P-acquisition repertoire in recovered MAGs.",
    "Supplementary Table 11 summarises P-acquisition functional modules and genes detected in recovered metagenome-assembled genomes (MAGs), and Supplementary Data 1 provides the full machine-readable MAG-level gene and module matrices.",
    "MAG_ID is a standardised identifier constructed from compartment, duration, fertilisation regime, sample code and bin ID; Original_MAG retains the raw MAG/bin name.",
    "P_score denotes the number of selected P-acquisition genes detected in each MAG.",
    "The metagenome analysis was restricted to the 5-year subset and is used as functional context rather than causal evidence."
  )
)

readme <- tibble(
  Field = c(
    "Table title",
    "Input workbook",
    "Input module/source sheet",
    "Input gene presence sheet",
    "Input sample map sheet",
    "Input GTDB taxonomy sheet",
    "P_score definition",
    "Gene set",
    "Functional modules",
    "Interpretation"
  ),
  Description = c(
    "Supplementary Table 11 / Supplementary Data 1 | Genome-resolved P-acquisition repertoire in recovered MAGs",
    supp_xlsx,
    "P_MAG_summary_tax_compartment",
    "P_MAG_gene_presence",
    "MAG_sample_map",
    "GTDB_taxonomy_summary",
    "Number of selected P-acquisition genes detected in each MAG.",
    paste(gene_sets, collapse = ", "),
    paste(unname(module_labels), collapse = "; "),
    "Functional context only; not used as causal evidence."
  )
)

# -------------------------
# 14) Write Excel: Supplementary Data 1
# -------------------------
wb <- createWorkbook()

addWorksheet(wb, "README")
writeData(wb, "README", readme)

addWorksheet(wb, "Caption")
writeData(wb, "Caption", caption)

addWorksheet(wb, "S11a_MAG_summary_Word")
writeData(wb, "S11a_MAG_summary_Word", s11a_MAG_summary_Word)

addWorksheet(wb, "S11b_module_summary_Word")
writeData(wb, "S11b_module_summary_Word", s11b_module_summary_Word)

addWorksheet(wb, "S11c_gene_summary_Word")
writeData(wb, "S11c_gene_summary_Word", s11c_gene_summary_Word)

addWorksheet(wb, "SD1_MAG_gene_matrix")
writeData(wb, "SD1_MAG_gene_matrix", gene_matrix)

addWorksheet(wb, "SD1_MAG_module_matrix")
writeData(wb, "SD1_MAG_module_matrix", module_matrix)

addWorksheet(wb, "SD1_MAG_gene_long")
writeData(wb, "SD1_MAG_gene_long", long_gene_presence)

addWorksheet(wb, "SD1_MAG_module_long")
writeData(wb, "SD1_MAG_module_long", long_module_presence)

addWorksheet(wb, "SD1_gene_detection_summary")
writeData(wb, "SD1_gene_detection_summary", gene_detection_summary)

addWorksheet(wb, "SD1_module_detection_summary")
writeData(wb, "SD1_module_detection_summary", module_detection_summary)

addWorksheet(wb, "SD1_module_cooccurrence")
writeData(wb, "SD1_module_cooccurrence", module_cooccurrence)

addWorksheet(wb, "SD1_gene_module_map")
writeData(wb, "SD1_gene_module_map", gene_module_map)

addWorksheet(wb, "SD1_sample_map_used")
writeData(wb, "SD1_sample_map_used", sample_map)

addWorksheet(wb, "SD1_GTDB_taxonomy_used")
writeData(wb, "SD1_GTDB_taxonomy_used", gtdb_all)

for (sh in names(wb)) {
  freezePane(wb, sh, firstRow = TRUE)
  setColWidths(wb, sh, cols = 1:100, widths = "auto")
}

num_style <- createStyle(numFmt = "0.000")
sheet_objects <- list(
  S11a_MAG_summary_Word = s11a_MAG_summary_Word,
  S11b_module_summary_Word = s11b_module_summary_Word,
  S11c_gene_summary_Word = s11c_gene_summary_Word,
  SD1_MAG_gene_matrix = gene_matrix,
  SD1_MAG_module_matrix = module_matrix,
  SD1_MAG_gene_long = long_gene_presence,
  SD1_MAG_module_long = long_module_presence,
  SD1_gene_detection_summary = gene_detection_summary,
  SD1_module_detection_summary = module_detection_summary,
  SD1_module_cooccurrence = module_cooccurrence,
  SD1_gene_module_map = gene_module_map,
  SD1_sample_map_used = sample_map,
  SD1_GTDB_taxonomy_used = gtdb_all
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

tryCatch({
  saveWorkbook(wb, out_xlsx, overwrite = TRUE)
  message("Excel written: ", out_xlsx)
}, error = function(e) {
  stop("Excel write failed: ", e$message)
})

# -------------------------
# 15) Write Word: Supplementary Table 11
# -------------------------
make_ft <- function(df, font_size = 6.5) {
  flextable(df) %>%
    theme_booktabs() %>%
    fontsize(size = font_size, part = "all") %>%
    fontsize(size = font_size + 0.5, part = "header") %>%
    align(align = "center", part = "all") %>%
    align(
      j = intersect(
        c(
          "MAG_ID", "Original_MAG", "Source", "Taxonomy",
          "P_modules", "P_genes",
          "Functional_module", "Functional_interpretation", "Gene"
        ),
        names(df)
      ),
      align = "left",
      part = "body"
    ) %>%
    autofit()
}

doc <- read_docx()

doc <- body_add_par(
  doc,
  "Supplementary Table 11 / Supplementary Data 1 | Genome-resolved P-acquisition repertoire in recovered MAGs.",
  style = "heading 1"
)

doc <- body_add_par(
  doc,
  paste(
    "Supplementary Table 11 summarises P-acquisition functional modules and genes detected in recovered",
    "metagenome-assembled genomes (MAGs), and Supplementary Data 1 provides the full machine-readable",
    "MAG-level gene and module matrices. MAG_ID is a standardised identifier constructed from compartment,",
    "duration, fertilisation regime, sample code and bin ID; Original_MAG retains the raw MAG/bin name.",
    "P_score denotes the number of selected P-acquisition genes detected in each MAG.",
    "The metagenome analysis was restricted to the 5-year subset and is used as functional context rather than causal evidence."
  ),
  style = "Normal"
)

doc <- body_add_par(doc, "a, MAG-level P-acquisition repertoire.", style = "heading 2")
doc <- body_add_flextable(doc, make_ft(s11a_MAG_summary_Word, font_size = 5.8))

doc <- body_add_par(doc, "b, P-acquisition module detection across P-carrying MAGs.", style = "heading 2")
doc <- body_add_flextable(doc, make_ft(s11b_module_summary_Word, font_size = 7.2))

doc <- body_add_par(doc, "c, P-acquisition gene detection across P-carrying MAGs.", style = "heading 2")
doc <- body_add_flextable(doc, make_ft(s11c_gene_summary_Word, font_size = 6.8))

print(doc, target = out_docx)

cat("\nSupplementary Table 11 / Supplementary Data 1 generated successfully:\n")
cat(out_docx, "\n")
cat(out_xlsx, "\n\n")

cat("P-carrying MAGs:", n_p_mags, "\n")
cat("Rows in S11a:", nrow(s11a_MAG_summary_Word), "\n")
cat("Rows in S11b:", nrow(s11b_module_summary_Word), "\n")
cat("Rows in S11c:", nrow(s11c_gene_summary_Word), "\n")
cat("MAGs without sample-map match:", sum(is.na(mag_all$MAG_ID)), "\n")
cat("MAGs without taxonomy:", sum(is.na(mag_all$Taxonomy)), "\n")
