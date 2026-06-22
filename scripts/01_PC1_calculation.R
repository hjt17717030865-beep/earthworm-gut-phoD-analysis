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
# Compute gut16S_PC1 + ChemPC1 and export tables
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(writexl)
  library(tibble)
})

# -------------------------
# 0) Paths
# -------------------------
supp_xlsx <- DATA_XLSX
out_dir <- OUTPUT_DIR
out_xlsx  <- file.path(out_dir, "PC1_gut16S_PC1_AI.xlsx")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# -------------------------
gut_raw <- readxl::read_excel(supp_xlsx, sheet = "16S_gut_ASV_table") %>%
  as.data.frame()

meta_cols <- c("PlotID", "Duration", "Regime", "Compartment")
asv <- as.matrix(gut_raw[, setdiff(colnames(gut_raw), meta_cols)])
storage.mode(asv) <- "numeric"
rownames(asv) <- trimws(as.character(gut_raw$PlotID))
rn <- rownames(asv)

# -------------------------
# -------------------------
master <- readxl::read_excel(supp_xlsx, sheet = "Master_all_data") %>%
  as.data.frame() %>%
  mutate(PlotID = trimws(as.character(PlotID)))

idx <- match(rn, master$PlotID)
if (anyNA(idx)) stop("PlotID in ASV table not found in Master_all_data: ",
                     paste(rn[is.na(idx)], collapse = ", "))

master2 <- master[idx, , drop = FALSE]

# -------------------------
# 3) Factors
# -------------------------
Duration2 <- trimws(as.character(master2$Duration))
Duration2[Duration2 == "3y"] <- "5y"
Duration2 <- factor(Duration2, levels = c("5y", "8y", "10y"))

Regime <- factor(trimws(as.character(master2$Regime)),
                 levels = c("CK", "NPK", "NPKOM"))

# -------------------------
# 4) gut16S_PC1
# -------------------------
asv_mat <- as.matrix(asv) + 1
log_mat <- log(asv_mat)
clr_mat <- log_mat - rowMeans(log_mat)

pca_gut    <- prcomp(clr_mat, center = TRUE, scale. = FALSE)
gut16S_PC1 <- as.numeric(pca_gut$x[, 1])

gut_var_explained <- (pca_gut$sdev^2) / sum(pca_gut$sdev^2)
gut_pca_variance <- tibble(
  PC                 = paste0("PC", seq_along(gut_var_explained)),
  Eigenvalue         = pca_gut$sdev^2,
  Proportion         = gut_var_explained,
  Percent            = 100 * gut_var_explained,
  Cumulative         = cumsum(gut_var_explained),
  Cumulative_percent = 100 * cumsum(gut_var_explained)
)

# -------------------------
# 5) ChemPC1
# -------------------------
chem_vars <- c("pH", "water content(%)", "TOC(%)",
               "NH4+(mg/kg)", "NO3-(mg/kg)", "TN(mg/kg)",
               "Olsen-P(mg/kg)", "TP(mg/kg)")

miss_chem <- setdiff(chem_vars, names(master2))
if (length(miss_chem) > 0)
  stop("Missing columns in Master_all_data: ", paste(miss_chem, collapse = ", "))

chem_mat <- master2 %>%
  select(all_of(chem_vars)) %>%
  mutate(across(everything(), as.numeric))

pca_chem <- prcomp(chem_mat, center = TRUE, scale. = TRUE)
ChemPC1  <- as.numeric(pca_chem$x[, 1])

chem_var_explained <- (pca_chem$sdev^2) / sum(pca_chem$sdev^2)
chem_pca_variance <- tibble(
  PC                 = paste0("PC", seq_along(chem_var_explained)),
  Eigenvalue         = pca_chem$sdev^2,
  Proportion         = chem_var_explained,
  Percent            = 100 * chem_var_explained,
  Cumulative         = cumsum(chem_var_explained),
  Cumulative_percent = 100 * cumsum(chem_var_explained)
)

chem_pc1_load    <- pca_chem$rotation[, 1]
chem_pc1_contrib <- (chem_pc1_load^2) / sum(chem_pc1_load^2)
chem_pc1_loadings <- tibble(
  Variable                  = names(chem_pc1_load),
  Loading_PC1               = as.numeric(chem_pc1_load),
  Abs_loading_PC1           = abs(as.numeric(chem_pc1_load)),
  Contribution_PC1          = as.numeric(chem_pc1_contrib),
  Contribution_PC1_percent  = 100 * as.numeric(chem_pc1_contrib)
) %>% arrange(desc(Contribution_PC1_percent))

# -------------------------
# -------------------------
if (!("AI" %in% names(master2)))
  stop("Master_all_data is missing the AI column")

base_tbl <- tibble(
  PlotID     = rn,
  Duration   = factor(as.character(Duration2), levels = c("5y","8y","10y")),
  Regime     = factor(as.character(Regime),    levels = c("CK","NPK","NPKOM")),
  AI         = as.numeric(master2$AI),
  gut16S_PC1 = gut16S_PC1,
  ChemPC1    = ChemPC1
) %>%
  filter(is.finite(AI), is.finite(gut16S_PC1), is.finite(ChemPC1))

AI_res         <- resid(lm(AI         ~ ChemPC1 + Duration * Regime, data = base_tbl))
gut16S_PC1_res <- resid(lm(gut16S_PC1 ~ ChemPC1 + Duration * Regime, data = base_tbl))

out_tbl <- base_tbl %>%
  mutate(AI_res         = as.numeric(AI_res),
         gut16S_PC1_res = as.numeric(gut16S_PC1_res))

# -------------------------
# -------------------------
cat("\n================ PCA summary ================\n")
cat(sprintf("gut16S PC1 explained variance: %.2f%%\n", gut_pca_variance$Percent[1]))
cat(sprintf("Chem   PC1 explained variance: %.2f%%\n", chem_pca_variance$Percent[1]))
cat("Top contributors to ChemPC1:\n")
print(chem_pc1_loadings)

# -------------------------
# -------------------------
writexl::write_xlsx(
  list(
    fig4_inputs       = out_tbl,
    gut_pca_variance  = gut_pca_variance,
    chem_pca_variance = chem_pca_variance,
    chem_pc1_loadings = chem_pc1_loadings
  ),
  out_xlsx
)

cat("\nSaved to:\n", out_xlsx, "\n")
