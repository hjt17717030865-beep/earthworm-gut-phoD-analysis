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
# Supplementary Table 9 | Robustness analyses of phoD amplification
# and soil-gut phoD community coupling
#
# Purpose:
#   Generate NC-ready robustness tables corresponding to:
#     Fig. S2 | Chronosequence robustness
#     Fig. S4 | Metric robustness of phoD amplification
#     Fig. S6 | Robustness of soil-gut phoD community dissimilarity
#
# Output:
#   Supplementary_Table_9_Robustness_phoD_amplification_and_coupling.xlsx
#   Supplementary_Table_9_Robustness_phoD_amplification_and_coupling.docx
#
# Input:
#
# Required columns:
#   PlotID, Duration, Regime, AI, gut16S_PC1, ChemPC1,
#   D_phoD, D_phoD_BC,
#   logphoD_gut, logphoD_soil, log16S_gut, log16S_soil
#
# Alternative accepted names:
#   gut_logphoD, soil_logphoD, gut_log16S, soil_log16S
#
# Notes:
#   - Duration is harmonised so that 3y is treated as 5y.
#   - All robustness models are interpreted as associational.
#   - Lower D_phoD indicates stronger soil-gut phoD community coupling.
# ============================================================

# -------------------------
# 0) Packages
# -------------------------
pkg_needed <- c(
  "readxl", "dplyr", "tibble", "tidyr", "purrr", "broom",
  "openxlsx", "officer", "flextable"
)

pkg_to_install <- pkg_needed[!sapply(pkg_needed, requireNamespace, quietly = TRUE)]
if (length(pkg_to_install) > 0) {
  install.packages(pkg_to_install)
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(broom)
  library(openxlsx)
  library(officer)
  library(flextable)
})

# -------------------------
# 1) Paths
# -------------------------
out_dir <- OUTPUT_DIR
in_xlsx  <- DATA_XLSX
in_sheet <- "Master_all_data"

out_xlsx <- file.path(
  out_dir,
  "Supplementary_Table_9_Robustness_phoD_amplification_and_coupling.xlsx"
)

out_docx <- file.path(
  out_dir,
  "Supplementary_Table_9_Robustness_phoD_amplification_and_coupling.docx"
)

if (!file.exists(in_xlsx)) {
  stop("Input file not found: ", in_xlsx)
}

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# -------------------------
# 2) Helper functions
# -------------------------
clean_duration <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\s+", "", x)
  x <- gsub("years|Years|year|Year|yr|YR", "y", x)
  x <- gsub("^3$", "5y", x)
  x <- gsub("^5$", "5y", x)
  x <- gsub("^8$", "8y", x)
  x <- gsub("^10$", "10y", x)
  x <- ifelse(x == "3y", "5y", x)
  factor(x, levels = c("5y", "8y", "10y"))
}

clean_regime <- function(x) {
  x <- trimws(as.character(x))
  factor(x, levels = c("CK", "NPK", "NPKOM"))
}

pick_first_existing <- function(df, candidates, label) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) {
    stop(
      "Missing required column for ", label, ". Tried: ",
      paste(candidates, collapse = ", ")
    )
  }
  hit[1]
}

fmt_p <- function(p) {
  ifelse(
    is.na(p), NA_character_,
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

fmt_num <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    NA_character_,
    ifelse(
      abs(x) >= 0.001 | x == 0,
      sprintf(paste0("%.", digits, "f"), x),
      sprintf("%.2e", x)
    )
  )
}

model_glance_safe <- function(fit) {
  gl <- broom::glance(fit)
  
  tibble(
    N = nobs(fit),
    Df_residual = df.residual(fit),
    R2 = unname(gl$r.squared),
    Adj_R2 = unname(gl$adj.r.squared),
    Model_P_value = unname(gl$p.value)
  )
}

tidy_lm_terms <- function(fit, figure, panel, analysis, model_id,
                          response_label, predictor_terms,
                          notes = NA_character_) {
  gl <- model_glance_safe(fit)
  
  broom::tidy(fit, conf.int = TRUE) %>%
    filter(term %in% predictor_terms) %>%
    mutate(
      Figure = figure,
      Panel = panel,
      Analysis = analysis,
      Model_ID = model_id,
      Response = response_label,
      Formula = paste(deparse(formula(fit)), collapse = " "),
      N = gl$N,
      Df_residual = gl$Df_residual,
      R2 = gl$R2,
      Adj_R2 = gl$Adj_R2,
      Model_P_value = gl$Model_P_value,
      Notes = notes
    ) %>%
    transmute(
      Figure, Panel, Analysis, Model_ID, Response, Formula,
      Term = term,
      Estimate = estimate,
      SE = std.error,
      CI_low = conf.low,
      CI_high = conf.high,
      t = statistic,
      P_value = p.value,
      P = fmt_p(p.value),
      N, Df_residual, R2, Adj_R2, Model_P_value, Notes
    )
}

make_word_compact <- function(df) {
  df %>%
    transmute(
      Figure,
      Panel,
      Analysis,
      Response,
      Term,
      Estimate = round(Estimate, 6),
      SE       = round(SE, 6),
      CI = paste0("[", fmt_num(CI_low, 6), ", ", fmt_num(CI_high, 6), "]"),
      t = round(t, 3),
      P,
      R2 = round(R2, 3),
      N,
      Notes
    )
}

# -------------------------
# 3) Read and prepare data
# -------------------------
dat0 <- readxl::read_excel(in_xlsx, sheet = in_sheet)

required_basic <- c(
  "PlotID", "Duration", "Regime",
  "AI", "gut16S_PC1", "ChemPC1",
  "D_phoD", "D_phoD_BC"
)

miss_basic <- setdiff(required_basic, names(dat0))
if (length(miss_basic) > 0) {
  stop("Input table missing columns: ", paste(miss_basic, collapse = ", "))
}

col_gut_phoD  <- pick_first_existing(dat0, c("logphoD_gut",  "gut_logphoD"),  "gut phoD")
col_soil_phoD <- pick_first_existing(dat0, c("logphoD_soil", "soil_logphoD"), "soil phoD")
col_gut_16S   <- pick_first_existing(dat0, c("log16S_gut",   "gut_log16S"),   "gut 16S")
col_soil_16S  <- pick_first_existing(dat0, c("log16S_soil",  "soil_log16S"),  "soil 16S")

dat <- dat0 %>%
  transmute(
    PlotID = trimws(as.character(PlotID)),
    Duration = clean_duration(Duration),
    Regime = clean_regime(Regime),
    AI = as.numeric(AI),
    gut16S_PC1 = as.numeric(gut16S_PC1),
    ChemPC1 = as.numeric(ChemPC1),
    D_phoD = as.numeric(D_phoD),
    D_phoD_BC = as.numeric(D_phoD_BC),
    gut_phoD = as.numeric(.data[[col_gut_phoD]]),
    soil_phoD = as.numeric(.data[[col_soil_phoD]]),
    gut_16S = as.numeric(.data[[col_gut_16S]]),
    soil_16S = as.numeric(.data[[col_soil_16S]])
  ) %>%
  mutate(
    d_phoD = gut_phoD - soil_phoD,
    d_16S = gut_16S - soil_16S
  ) %>%
  filter(!is.na(Duration), !is.na(Regime)) %>%
  droplevels()

if (nrow(dat) == 0) {
  stop("No valid rows after cleaning Duration and Regime.")
}

# ============================================================
# 4) Fig. S2 | Chronosequence robustness
# ============================================================

# -------------------------
# S2a: within-duration contrasts, AI ~ Regime
# -------------------------
s2a_within_duration <- map_dfr(levels(dat$Duration), function(d0) {
  dd <- dat %>%
    filter(Duration == d0, is.finite(AI)) %>%
    droplevels()
  
  fit <- lm(AI ~ Regime, data = dd)
  gl <- model_glance_safe(fit)
  
  broom::tidy(fit, conf.int = TRUE) %>%
    filter(term %in% c("RegimeNPK", "RegimeNPKOM")) %>%
    transmute(
      Figure = "Fig. S2",
      Panel = "a",
      Analysis = "Within-duration contrasts",
      Model_ID = paste0("S2a_AI_Regime_within_", d0),
      Duration_group = d0,
      Response = "AI",
      Formula = "AI ~ Regime, fitted within each duration",
      Term = case_when(
        term == "RegimeNPK" ~ "NPK - CK",
        term == "RegimeNPKOM" ~ "NPKOM - CK",
        TRUE ~ term
      ),
      Estimate = estimate,
      SE = std.error,
      CI_low = conf.low,
      CI_high = conf.high,
      t = statistic,
      P_value = p.value,
      P = fmt_p(p.value),
      N = gl$N,
      Df_residual = gl$Df_residual,
      R2 = gl$R2,
      Adj_R2 = gl$Adj_R2,
      Model_P_value = gl$Model_P_value,
      Notes = "Contrast estimates are relative to CK within each duration."
    )
})
# -------------------------
# S2b: leave-one-duration-out, AI ~ gut16S_PC1 + ChemPC1
# -------------------------
s2b_leave_one_duration_out <- map_dfr(levels(dat$Duration), function(d0) {
  dd <- dat %>%
    filter(
      Duration != d0,
      is.finite(AI),
      is.finite(gut16S_PC1),
      is.finite(ChemPC1)
    ) %>%
    droplevels()
  
  fit <- lm(AI ~ gut16S_PC1 + ChemPC1, data = dd)
  
  tidy_lm_terms(
    fit = fit,
    figure = "Fig. S2",
    panel = "b",
    analysis = "Leave-one-duration-out",
    model_id = paste0("S2b_AI_gut16S_PC1_ChemPC1_without_", d0),
    response_label = "AI",
    predictor_terms = c("gut16S_PC1", "ChemPC1"),
    notes = paste0("Model fitted after excluding duration ", d0, ".")
  ) %>%
    mutate(Excluded_duration = d0, .after = Model_ID)
})

# -------------------------
# S2c: duration-centred model
# -------------------------
dat_centered <- dat %>%
  filter(
    is.finite(AI),
    is.finite(gut16S_PC1),
    is.finite(ChemPC1)
  ) %>%
  group_by(Duration) %>%
  mutate(
    AI_c = AI - mean(AI, na.rm = TRUE),
    gut16S_PC1_c = gut16S_PC1 - mean(gut16S_PC1, na.rm = TRUE),
    ChemPC1_c = ChemPC1 - mean(ChemPC1, na.rm = TRUE)
  ) %>%
  ungroup()

fit_s2c <- lm(AI_c ~ gut16S_PC1_c + ChemPC1_c, data = dat_centered)

s2c_duration_centered <- tidy_lm_terms(
  fit = fit_s2c,
  figure = "Fig. S2",
  panel = "c",
  analysis = "Duration-centred model",
  model_id = "S2c_AIc_gut16S_PC1c_ChemPC1c",
  response_label = "AI_c",
  predictor_terms = c("gut16S_PC1_c", "ChemPC1_c"),
  notes = "Variables were centred within each duration before model fitting."
)

s2_chronosequence <- bind_rows(
  s2a_within_duration,
  s2b_leave_one_duration_out,
  s2c_duration_centered
)

# ============================================================
# 5) Fig. S4 | Metric robustness of phoD amplification
# ============================================================

fit_metric_response <- function(response_var, response_label) {
  dd <- dat %>%
    filter(
      is.finite(.data[[response_var]]),
      is.finite(gut16S_PC1),
      is.finite(ChemPC1)
    )
  
  fit <- lm(
    as.formula(paste(response_var, "~ gut16S_PC1 + ChemPC1")),
    data = dd
  )
  
  tidy_lm_terms(
    fit = fit,
    figure = "Fig. S4",
    panel = "all",
    analysis = "Metric robustness of phoD amplification",
    model_id = paste0("S4_", response_var, "_gut16S_PC1_ChemPC1"),
    response_label = response_label,
    predictor_terms = c("gut16S_PC1", "ChemPC1"),
    notes = "Same predictor structure fitted across alternative response metrics."
  )
}

s4_metric_robustness <- bind_rows(
  fit_metric_response("AI", "AI"),
  fit_metric_response("d_phoD", "Delta log10 phoD"),
  fit_metric_response("d_16S", "Delta log10 16S")
)

# ============================================================
# 6) Fig. S6 | Distance robustness of soil-gut phoD coupling
# ============================================================

dat_s6 <- dat %>%
  transmute(
    PlotID,
    Duration,
    Regime,
    AI,
    D_phoD_Aitchison = D_phoD,
    D_phoD_BrayCurtis = D_phoD_BC
  ) %>%
  filter(
    is.finite(AI),
    is.finite(D_phoD_Aitchison),
    is.finite(D_phoD_BrayCurtis)
  ) %>%
  pivot_longer(
    cols = c(D_phoD_Aitchison, D_phoD_BrayCurtis),
    names_to = "Distance_metric",
    values_to = "D_phoD_value"
  ) %>%
  mutate(
    Distance_metric = case_when(
      Distance_metric == "D_phoD_Aitchison" ~ "Aitchison",
      Distance_metric == "D_phoD_BrayCurtis" ~ "Bray-Curtis",
      TRUE ~ Distance_metric
    ),
    Distance_metric = factor(Distance_metric, levels = c("Aitchison", "Bray-Curtis"))
  )

fit_distance_metric <- function(metric_name) {
  dd <- dat_s6 %>% filter(Distance_metric == metric_name)
  
  fit <- lm(scale(D_phoD_value) ~ AI, data = dd)
  
  tidy_lm_terms(
    fit = fit,
    figure = "Fig. S6",
    panel = "a-b",
    analysis = "Distance-metric robustness of phoD coupling",
    model_id = paste0("S6_scaled_DphoD_", gsub("-", "_", metric_name), "_AI"),
    response_label = paste0("Scaled D_phoD (", metric_name, ")"),
    predictor_terms = c("AI"),
    notes = "Response was scaled within distance metric, matching the coefficient comparison in Fig. S6."
  ) %>%
    mutate(Distance_metric = metric_name, .after = Model_ID)
}

s6_distance_robustness <- bind_rows(
  fit_distance_metric("Aitchison"),
  fit_distance_metric("Bray-Curtis")
)

# ============================================================
# 7) Compact Word-ready tables
# ============================================================

s9a_chronosequence_compact <- s2_chronosequence %>%
  mutate(
    Context = case_when(
      Analysis == "Within-duration contrasts" ~ paste0("Duration ", Duration_group),
      Analysis == "Leave-one-duration-out" ~ paste0("Excluded ", Excluded_duration),
      Analysis == "Duration-centred model" ~ "Within-duration centred",
      TRUE ~ Analysis
    )
  ) %>%
  transmute(
    Figure,
    Panel,
    Analysis,
    Context,
    Response,
    Term,
    Estimate = fmt_num(Estimate, 4),
    SE       = fmt_num(SE, 4),
    CI = paste0("[", fmt_num(CI_low, 4), ", ", fmt_num(CI_high, 4), "]"),
    t = round(t, 3),
    P,
    R2 = round(R2, 3),
    N
  )

s9b_metric_compact <- s4_metric_robustness %>%
  transmute(
    Figure,
    Panel,
    Analysis,
    Response,
    Term,
    Estimate = fmt_num(Estimate, 4),
    SE       = fmt_num(SE, 4),
    CI = paste0("[", fmt_num(CI_low, 4), ", ", fmt_num(CI_high, 4), "]"),
    t = round(t, 3),
    P,
    R2 = round(R2, 3),
    N
  )

s9c_coupling_compact <- s6_distance_robustness %>%
  transmute(
    Figure,
    Panel,
    Analysis,
    Distance_metric,
    Response,
    Term,
    Estimate = fmt_num(Estimate, 4),
    SE       = fmt_num(SE, 4),
    CI = paste0("[", fmt_num(CI_low, 4), ", ", fmt_num(CI_high, 4), "]"),
    t = round(t, 3),
    P,
    R2 = round(R2, 3),
    N
  )

s9_full <- bind_rows(
  s2_chronosequence %>% mutate(Table_panel = "a, Chronosequence robustness"),
  s4_metric_robustness %>% mutate(Table_panel = "b, Metric robustness"),
  s6_distance_robustness %>% mutate(Table_panel = "c, Coupling-distance robustness")
) %>%
  dplyr::select(Table_panel, everything())

# -------------------------
# 8) README / caption
# -------------------------
caption <- tibble(
  Caption = c(
    "Supplementary Table 9 | Robustness analyses of phoD amplification and soil-gut phoD community coupling.",
    "Panel a reports chronosequence-robustness analyses of phoD amplification, including within-duration regime contrasts, leave-one-duration-out models, and a duration-centred model.",
    "Panel b reports metric-robustness analyses using the same predictor structure for AI, Delta log10 phoD, and Delta log10 16S.",
    "Panel c reports robustness of the association between AI and soil-gut phoD community dissimilarity using Aitchison and Bray-Curtis distance metrics.",
    "Duration represents a site-bound chronosequence proxy. All models are interpreted as descriptive and associational. Lower D_phoD values indicate stronger soil-gut phoD community coupling."
  )
)

readme <- tibble(
  Field = c(
    "Table title",
    "Panel a",
    "Panel b",
    "Panel c",
    "Input file",
    "Duration harmonisation",
    "Interpretation"
  ),
  Description = c(
    "Supplementary Table 9 | Robustness analyses of phoD amplification and soil-gut phoD community coupling",
    "Chronosequence robustness: within-duration contrasts, leave-one-duration-out models, and duration-centred model.",
    "Metric robustness: AI, Delta log10 phoD, and Delta log10 16S fitted against gut16S_PC1 and ChemPC1.",
    "Coupling robustness: scaled D_phoD fitted against AI for Aitchison and Bray-Curtis distances.",
    in_xlsx,
    "3y is harmonised to 5y; factor order is 5y, 8y, 10y.",
    "All robustness analyses are associational. Lower D_phoD indicates stronger soil-gut phoD community coupling."
  )
)

# ============================================================
# 9) Write Excel workbook
# ============================================================

wb <- createWorkbook()

addWorksheet(wb, "README")
writeData(wb, "README", readme)

addWorksheet(wb, "Caption")
writeData(wb, "Caption", caption)

addWorksheet(wb, "S9a_chronosequence_Word")
writeData(wb, "S9a_chronosequence_Word", s9a_chronosequence_compact)

addWorksheet(wb, "S9b_metric_Word")
writeData(wb, "S9b_metric_Word", s9b_metric_compact)

addWorksheet(wb, "S9c_coupling_Word")
writeData(wb, "S9c_coupling_Word", s9c_coupling_compact)

addWorksheet(wb, "S9_full_machine_readable")
writeData(wb, "S9_full_machine_readable", s9_full)

addWorksheet(wb, "S9a_chronosequence_full")
writeData(wb, "S9a_chronosequence_full", s2_chronosequence)

addWorksheet(wb, "S9b_metric_full")
writeData(wb, "S9b_metric_full", s4_metric_robustness)

addWorksheet(wb, "S9c_coupling_full")
writeData(wb, "S9c_coupling_full", s6_distance_robustness)

for (sh in names(wb)) {
  freezePane(wb, sh, firstRow = TRUE)
  setColWidths(wb, sh, cols = 1:60, widths = "auto")
}

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

# ============================================================
# 10) Write Word-ready document
# ============================================================

make_ft <- function(df) {
  flextable(df) %>%
    theme_booktabs() %>%
    fontsize(size = 7.5, part = "all") %>%
    fontsize(size = 8, part = "header") %>%
    align(align = "center", part = "all") %>%
    autofit()
}

doc <- read_docx()

doc <- body_add_par(
  doc,
  "Supplementary Table 9 | Robustness analyses of phoD amplification and soil-gut phoD community coupling.",
  style = "heading 1"
)

doc <- body_add_par(
  doc,
  paste(
    "Panel a reports chronosequence-robustness analyses of phoD amplification, including within-duration regime contrasts,",
    "leave-one-duration-out models, and a duration-centred model.",
    "Panel b reports metric-robustness analyses using the same predictor structure for AI, Delta log10 phoD, and Delta log10 16S.",
    "Panel c reports robustness of the association between AI and soil-gut phoD community dissimilarity using Aitchison and Bray-Curtis distance metrics.",
    "Duration represents a site-bound chronosequence proxy.",
    "All models are interpreted as descriptive and associational.",
    "Lower D_phoD values indicate stronger soil-gut phoD community coupling."
  ),
  style = "Normal"
)

doc <- body_add_par(doc, "a, Chronosequence robustness of phoD amplification.", style = "heading 2")
doc <- body_add_flextable(doc, make_ft(s9a_chronosequence_compact))

doc <- body_add_par(doc, "b, Metric robustness of phoD amplification.", style = "heading 2")
doc <- body_add_flextable(doc, make_ft(s9b_metric_compact))

doc <- body_add_par(doc, "c, Distance-metric robustness of soil-gut phoD community coupling.", style = "heading 2")
doc <- body_add_flextable(doc, make_ft(s9c_coupling_compact))

print(doc, target = out_docx)

cat("\nSupplementary Table 9 generated successfully:\n")
cat(out_xlsx, "\n")
cat(out_docx, "\n\n")

cat("Panel sizes:\n")
cat("S9a rows:", nrow(s9a_chronosequence_compact), "\n")
cat("S9b rows:", nrow(s9b_metric_compact), "\n")
cat("S9c rows:", nrow(s9c_coupling_compact), "\n")
