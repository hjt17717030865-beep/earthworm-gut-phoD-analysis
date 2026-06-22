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

library(dplyr)
library(vegan)
library(readr)
library(readxl)

xlsx_path <- DATA_XLSX
# ============================================================
# ============================================================

soil_raw <- read_excel(xlsx_path, sheet = "phoD_soil_ASV_table")
gut_raw  <- read_excel(xlsx_path, sheet = "phoD_gut_ASV_table")

meta_cols <- c("PlotID", "Duration", "Regime", "Compartment")

soil_raw$Compartment <- "Soil"
gut_raw$Compartment  <- "Gut"

soil <- soil_raw %>% arrange(PlotID)
gut  <- gut_raw  %>% arrange(PlotID)

stopifnot(nrow(soil) == nrow(gut))
stopifnot(all(soil$PlotID == gut$PlotID))

non_asv <- c("Duration", "Regime", "Compartment")
soil_asv <- soil[, setdiff(colnames(soil), c("PlotID", non_asv)), drop = FALSE]
gut_asv  <- gut[,  setdiff(colnames(gut),  c("PlotID", non_asv)), drop = FALSE]

soil_asv <- as.data.frame(lapply(soil_asv, as.numeric))
gut_asv  <- as.data.frame(lapply(gut_asv,  as.numeric))

all_asv <- union(colnames(soil_asv), colnames(gut_asv))
align_asv <- function(df, all_names) {
  out <- matrix(0, nrow = nrow(df), ncol = length(all_names))
  colnames(out) <- all_names
  common <- intersect(colnames(df), all_names)
  out[, common] <- as.matrix(df[, common])
  as.data.frame(out)
}
soil_u <- align_asv(soil_asv, all_asv)
gut_u  <- align_asv(gut_asv,  all_asv)

# ============================================================
# ============================================================

clr <- function(x) {
  x <- as.matrix(x) + 1
  logx <- log(x)
  gm <- rowMeans(logx)
  sweep(logx, 1, gm)
}
soil_clr <- clr(soil_u)
gut_clr  <- clr(gut_u)

D_plot_Ait <- sqrt(rowSums((soil_clr - gut_clr)^2))

clr_all <- rbind(soil_clr, gut_clr)
pca <- prcomp(clr_all)

# ============================================================
# ============================================================

rel_abund <- function(df) {
  m <- as.matrix(df)
  sweep(m, 1, rowSums(m), "/")
}
soil_rel <- rel_abund(soil_u)
gut_rel  <- rel_abund(gut_u)

D_plot_BC <- sapply(seq_len(nrow(soil_rel)), function(i) {
  m <- rbind(soil_rel[i, ], gut_rel[i, ])
  as.numeric(vegdist(m, method = "bray"))
})

all_rel <- rbind(soil_rel, gut_rel)
set.seed(123)
nmds <- metaMDS(all_rel, distance = "bray", k = 2, trymax = 100, trace = FALSE)

# ============================================================
# ============================================================

meta_soil <- soil[, c("PlotID", "Duration", "Regime")] %>% mutate(Compartment = "Soil")
meta_gut  <- gut[,  c("PlotID", "Duration", "Regime")] %>% mutate(Compartment = "Gut")
meta_all  <- bind_rows(meta_soil, meta_gut)

# ============================================================
# ============================================================

ordination_ait <- meta_all
ordination_ait$PC1   <- pca$x[, 1]
ordination_ait$PC2   <- pca$x[, 2]
ordination_ait$NMDS1 <- nmds$points[, 1]
ordination_ait$NMDS2 <- nmds$points[, 2]

d_plot_df <- tibble(
  PlotID     = soil$PlotID,
  D_plot_Ait = as.numeric(D_plot_Ait),
  D_plot_BC  = as.numeric(D_plot_BC)
)

final_df <- ordination_ait %>%
  left_join(d_plot_df, by = "PlotID") %>%
  mutate(
    Duration    = factor(Duration,    levels = c("5y", "8y", "10y")),
    Regime      = factor(Regime,      levels = c("CK", "NPK", "NPKOM")),
    Compartment = factor(Compartment, levels = c("Gut", "Soil"))
  ) %>%
  arrange(Compartment, Duration, Regime, PlotID)

sapply(final_df, class)

write_csv(final_df, file.path(OUTPUT_DIR, "phoD_ordination_Dplot.csv"))
