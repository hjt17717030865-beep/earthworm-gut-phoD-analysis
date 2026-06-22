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
# Supplementary Table 8 | PERMANOVA and dispersion test results
# Read all ASV matrices from Supplementary Data.xlsx
# ============================================================

pkg_needed <- c(
  "dplyr", "tibble", "tidyr", "readxl", "openxlsx",
  "vegan", "officer", "flextable"
)

pkg_to_install <- pkg_needed[!sapply(pkg_needed, requireNamespace, quietly = TRUE)]
if (length(pkg_to_install) > 0) {
  install.packages(pkg_to_install)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(readxl)
  library(openxlsx)
  library(vegan)
  library(officer)
  library(flextable)
})

# -------------------------
# 1) Paths
# -------------------------
out_dir <- OUTPUT_DIR
main_xlsx <- DATA_XLSX

soil_16s_sheet <- "16S_soil_ASV_table"
gut_16s_sheet  <- "16S_gut_ASV_table"

include_phoD <- FALSE
soil_phod_sheet <- "soil_phoD_ASV_table"
gut_phod_sheet  <- "gut_phoD_ASV_table"

out_xlsx <- file.path(
  out_dir,
  "Supplementary_Table_8_PERMANOVA_dispersion_test_results.xlsx"
)

out_docx <- file.path(
  out_dir,
  "Supplementary_Table_8_PERMANOVA_dispersion_test_results.docx"
)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

if (!file.exists(main_xlsx)) {
  stop("Supplementary Data file not found: ", main_xlsx)
}

need_sheets <- c(soil_16s_sheet, gut_16s_sheet)

if (isTRUE(include_phoD)) {
  need_sheets <- c(need_sheets, soil_phod_sheet, gut_phod_sheet)
}

miss_sheets <- setdiff(need_sheets, readxl::excel_sheets(main_xlsx))

if (length(miss_sheets) > 0) {
  stop(
    "These sheets are missing from Supplementary Data.xlsx: ",
    paste(miss_sheets, collapse = ", ")
  )
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

clean_plot_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\s+", "", x)
  x <- gsub("\\.fastq(\\.gz)?$", "", x)
  x <- gsub("_R1.*$", "", x)
  x <- gsub("_R2.*$", "", x)
  x <- gsub("^Soil__", "", x)
  x <- gsub("^Gut__", "", x)
  x <- gsub("^([0-9]+y)EB", "\\1", x)
  x <- gsub("^([0-9]+y)E", "\\1", x)
  x <- gsub("^3y", "5y", x)
  x
}

get_duration <- function(x) {
  d <- sub("^([0-9]+y).*", "\\1", x)
  d <- ifelse(d == "3y", "5y", d)
  factor(d, levels = c("5y", "8y", "10y"))
}

get_regime <- function(x) {
  r <- dplyr::case_when(
    grepl("NPKOM", x) ~ "NPKOM",
    grepl("NPK", x) ~ "NPK",
    grepl("CK", x) ~ "CK",
    TRUE ~ NA_character_
  )
  factor(r, levels = c("CK", "NPK", "NPKOM"))
}

read_asv_sheet <- function(xlsx_path, sheet_name) {
  x <- readxl::read_excel(xlsx_path, sheet = sheet_name) %>%
    as.data.frame(check.names = FALSE)
  
  id_candidates <- c("SampleID", "sample_id", "PlotID", "plot_id", "ID", "id")
  id_col <- intersect(id_candidates, names(x))
  
  if (length(id_col) > 0) {
    id_col <- id_col[1]
  } else {
    id_col <- names(x)[1]
  }
  
  rn <- clean_plot_id(x[[id_col]])
  
  meta_cols <- intersect(
    c(
      "SampleID", "sample_id", "PlotID", "plot_id", "ID", "id",
      "Duration", "Regime", "Compartment", "Replicate", "Replicate_plot"
    ),
    names(x)
  )
  
  asv_df <- x[, setdiff(names(x), meta_cols), drop = FALSE]
  
  asv_df[] <- lapply(
    asv_df,
    function(z) suppressWarnings(as.numeric(as.character(z)))
  )
  
  keep_cols <- vapply(
    asv_df,
    function(z) any(is.finite(z)),
    logical(1)
  )
  
  asv_df <- asv_df[, keep_cols, drop = FALSE]
  
  mat <- as.matrix(asv_df)
  rownames(mat) <- rn
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  
  mat
}

pad_matrix <- function(mat, all_cols) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  
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

clr_transform <- function(mat, pseudocount = 1) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  
  mat <- mat + pseudocount
  log_mat <- log(mat)
  sweep(log_mat, 1, rowMeans(log_mat), "-")
}

make_paired_clr_object <- function(soil_mat, gut_mat) {
  rownames(soil_mat) <- clean_plot_id(rownames(soil_mat))
  rownames(gut_mat)  <- clean_plot_id(rownames(gut_mat))
  
  common_pairs <- sort(intersect(rownames(soil_mat), rownames(gut_mat)))
  
  if (length(common_pairs) == 0) {
    stop("No shared PlotID / PairID found between soil and gut matrices.")
  }
  
  soil_mat <- soil_mat[common_pairs, , drop = FALSE]
  gut_mat  <- gut_mat[common_pairs,  , drop = FALSE]
  
  all_asv <- union(colnames(soil_mat), colnames(gut_mat))
  
  soil_pad <- pad_matrix(soil_mat, all_asv)
  gut_pad  <- pad_matrix(gut_mat, all_asv)
  
  soil_clr <- clr_transform(soil_pad)
  gut_clr  <- clr_transform(gut_pad)
  
  rownames(soil_clr) <- paste0("Soil__", common_pairs)
  rownames(gut_clr)  <- paste0("Gut__",  common_pairs)
  
  clr_mat <- rbind(soil_clr, gut_clr)
  
  meta <- tibble(
    SampleID = rownames(clr_mat),
    Compartment = ifelse(grepl("^Soil__", SampleID), "Soil", "Gut"),
    PairID = sub("^(Soil|Gut)__", "", SampleID)
  ) %>%
    mutate(
      Compartment = factor(Compartment, levels = c("Soil", "Gut")),
      Duration = get_duration(PairID),
      Regime = get_regime(PairID)
    )
  
  list(
    clr_mat = clr_mat,
    meta = meta,
    common_pairs = common_pairs
  )
}

run_permanova_dispersion <- function(clr_obj,
                                     figure,
                                     panel,
                                     dataset,
                                     response,
                                     distance_name = "Aitchison",
                                     model_formula = "Distance ~ Compartment",
                                     permutations = 999,
                                     seed = 123,
                                     notes = NA_character_) {
  
  clr_mat <- clr_obj$clr_mat
  meta <- clr_obj$meta
  
  dist_mat <- dist(clr_mat, method = "euclidean")
  
  set.seed(seed)
  perm <- vegan::adonis2(
    dist_mat ~ Compartment,
    data = as.data.frame(meta),
    permutations = permutations,
    strata = meta$PairID
  )
  
  bd <- vegan::betadisper(dist_mat, group = meta$Compartment)
  
  set.seed(seed)
  bd_perm <- vegan::permutest(bd, permutations = permutations)
  
  perm_df <- as.data.frame(perm)
  perm_df$Term <- rownames(perm_df)
  
  permanova_out <- perm_df %>%
    as_tibble() %>%
    filter(Term != "Total") %>%
    transmute(
      Figure = figure,
      Panel = panel,
      Dataset = dataset,
      Response = response,
      Test = "PERMANOVA",
      Distance = distance_name,
      Model = model_formula,
      Term = Term,
      Df = Df,
      Sum_of_squares = SumOfSqs,
      R2 = R2,
      F_value = F,
      P_value = `Pr(>F)`,
      P = fmt_p(`Pr(>F)`),
      Permutations = permutations,
      N_samples = nrow(meta),
      N_pairs = length(unique(meta$PairID)),
      Strata = "PairID",
      Notes = notes
    )
  
  disp_tab <- as.data.frame(bd_perm$tab)
  disp_tab$Term <- rownames(disp_tab)
  
  dispersion_out <- disp_tab %>%
    as_tibble() %>%
    filter(Term == "Groups") %>%
    transmute(
      Figure = figure,
      Panel = panel,
      Dataset = dataset,
      Response = response,
      Test = "Dispersion test",
      Distance = distance_name,
      Model = "betadisper(distance, Compartment)",
      Term = "Compartment",
      Df = Df,
      Sum_of_squares = `Sum Sq`,
      R2 = NA_real_,
      F_value = F,
      P_value = `Pr(>F)`,
      P = fmt_p(`Pr(>F)`),
      Permutations = permutations,
      N_samples = nrow(meta),
      N_pairs = length(unique(meta$PairID)),
      Strata = NA_character_,
      Notes = "Permutation test of multivariate dispersion among compartments."
    )
  
  list(
    permanova = permanova_out,
    dispersion = dispersion_out
  )
}

# -------------------------
# 3) 16S paired soil-gut analysis
# -------------------------
soil_16s <- read_asv_sheet(main_xlsx, soil_16s_sheet)
gut_16s  <- read_asv_sheet(main_xlsx, gut_16s_sheet)

obj_16s <- make_paired_clr_object(
  soil_mat = soil_16s,
  gut_mat = gut_16s
)

res_16s <- run_permanova_dispersion(
  clr_obj = obj_16s,
  figure = "Fig. 3",
  panel = "a",
  dataset = "16S",
  response = "Paired soil-gut 16S community composition",
  distance_name = "Aitchison",
  model_formula = "Aitchison distance ~ Compartment",
  permutations = 999,
  seed = 123,
  notes = "PERMANOVA used restricted permutations stratified by PairID."
)

# -------------------------
# 4) Optional phoD paired soil-gut analysis
# -------------------------
permanova_all <- res_16s$permanova
dispersion_all <- res_16s$dispersion

if (isTRUE(include_phoD)) {
  
  soil_phod <- read_asv_sheet(main_xlsx, soil_phod_sheet)
  gut_phod  <- read_asv_sheet(main_xlsx, gut_phod_sheet)
  
  obj_phod <- make_paired_clr_object(
    soil_mat = soil_phod,
    gut_mat = gut_phod
  )
  
  res_phod <- run_permanova_dispersion(
    clr_obj = obj_phod,
    figure = "Supplementary",
    panel = "optional",
    dataset = "phoD",
    response = "Paired soil-gut phoD community composition",
    distance_name = "Aitchison",
    model_formula = "Aitchison distance ~ Compartment",
    permutations = 999,
    seed = 123,
    notes = "Optional analysis; include only if reported in the manuscript or Supplementary Information."
  )
  
  permanova_all <- bind_rows(permanova_all, res_phod$permanova)
  dispersion_all <- bind_rows(dispersion_all, res_phod$dispersion)
}

# -------------------------
# 5) Word-ready compact table
# -------------------------
s8_compact <- bind_rows(permanova_all, dispersion_all) %>%
  filter(Term != "Residual") %>%
  transmute(
    Figure,
    Panel,
    Dataset,
    Test,
    Distance,
    Term,
    Df,
    F = round(F_value, 3),
    R2 = ifelse(is.na(R2), NA, round(R2, 3)),
    P,
    N_samples,
    N_pairs,
    Permutations,
    Strata
  )

# -------------------------
# 6) README / caption
# -------------------------
readme <- tibble(
  Field = c(
    "Table title",
    "Main analysis",
    "Distance",
    "PERMANOVA",
    "Dispersion test",
    "Pairing",
    "Permutation number",
    "Interpretation"
  ),
  Description = c(
    "Supplementary Table 8 | PERMANOVA and dispersion test results",
    "Fig. 3a paired soil-gut 16S community composition.",
    "Aitchison distance, calculated as Euclidean distance after CLR transformation with a pseudocount of 1.",
    "adonis2(distance ~ Compartment), with permutations stratified by PairID.",
    "betadisper followed by permutest to assess homogeneity of multivariate dispersion.",
    "Soil and gut samples are paired by PlotID / PairID.",
    "999 permutations; random seed = 123.",
    "PERMANOVA tests compositional separation between soil and gut communities; the dispersion test checks whether this result is confounded by differences in within-group dispersion."
  )
)

caption <- tibble(
  Caption = c(
    "Supplementary Table 8 | PERMANOVA and dispersion test results.",
    "PERMANOVA and dispersion-test outputs for paired soil-gut amplicon community comparisons.",
    "Aitchison distances were calculated as Euclidean distances after centred log-ratio transformation.",
    "PERMANOVA was performed using 999 permutations stratified by PairID to preserve the paired soil-gut design.",
    "Differences in multivariate dispersion were assessed using betadisper followed by permutation tests."
  )
)

# -------------------------
# 7) Write Excel
# -------------------------
wb <- createWorkbook()

addWorksheet(wb, "README")
writeData(wb, "README", readme)

addWorksheet(wb, "Caption")
writeData(wb, "Caption", caption)

addWorksheet(wb, "S8_compact_Word_ready")
writeData(wb, "S8_compact_Word_ready", s8_compact)

addWorksheet(wb, "S8a_PERMANOVA_full")
writeData(wb, "S8a_PERMANOVA_full", permanova_all)

addWorksheet(wb, "S8b_Dispersion_full")
writeData(wb, "S8b_Dispersion_full", dispersion_all)

for (sh in names(wb)) {
  freezePane(wb, sh, firstRow = TRUE)
  setColWidths(wb, sh, cols = 1:40, widths = "auto")
}

num_style <- createStyle(numFmt = "0.000")

sheet_objects <- list(
  S8_compact_Word_ready = s8_compact,
  S8a_PERMANOVA_full = permanova_all,
  S8b_Dispersion_full = dispersion_all
)

for (sh in names(sheet_objects)) {
  dat_sh <- sheet_objects[[sh]]
  numeric_cols <- which(vapply(dat_sh, is.numeric, logical(1)))
  
  if (length(numeric_cols) > 0 && nrow(dat_sh) > 0) {
    addStyle(
      wb,
      sh,
      style = num_style,
      rows = 2:(nrow(dat_sh) + 1),
      cols = numeric_cols,
      gridExpand = TRUE,
      stack = TRUE
    )
  }
}

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

# -------------------------
# 8) Write Word-ready document
# -------------------------
ft <- flextable(s8_compact)
ft <- autofit(ft)
ft <- fontsize(ft, size = 8, part = "all")
ft <- fontsize(ft, size = 8.5, part = "header")
ft <- align(ft, align = "center", part = "all")
ft <- align(
  ft,
  j = c("Figure", "Panel", "Dataset", "Test", "Distance", "Term", "Strata"),
  align = "left",
  part = "body"
)
ft <- theme_booktabs(ft)

doc <- read_docx()

doc <- body_add_par(
  doc,
  "Supplementary Table 8 | PERMANOVA and dispersion test results.",
  style = "heading 1"
)

doc <- body_add_par(
  doc,
  paste(
    "PERMANOVA and dispersion-test outputs for paired soil-gut amplicon community comparisons.",
    "Aitchison distances were calculated as Euclidean distances after centred log-ratio transformation.",
    "PERMANOVA was performed using 999 permutations stratified by PairID.",
    "Differences in multivariate dispersion were assessed using betadisper followed by permutation tests."
  ),
  style = "Normal"
)

doc <- body_add_par(
  doc,
  "a, PERMANOVA and dispersion-test results.",
  style = "heading 2"
)

doc <- body_add_flextable(doc, ft)

print(doc, target = out_docx)

cat("\nSupplementary Table 8 generated successfully:\n")
cat(out_xlsx, "\n")
cat(out_docx, "\n\n")

cat("Rows in compact table:\n")
print(s8_compact)
