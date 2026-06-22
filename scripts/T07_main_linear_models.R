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
# Supplementary Table 7 | Statistical outputs for main linear models
#
# Purpose:
#   Generate NC-compliant statistical output tables for the main
#   linear models used in Fig. 2, Fig. 3c, Fig. 4b-d, and Fig. 4e.
#
# Output:
#   Supplementary_Table_7_Statistical_outputs_main_linear_models.xlsx
#
# Input:
#
# Required columns:
#   PlotID, Duration, Regime,
#   ChemPC1, Olsen-P(mg/kg), pH,
#   gut16S_PC1, AI, D_phoD
#
# Notes:
#   - Duration is harmonised so that 3y is treated as 5y.
#   - Type-II ANOVA is used for models with Duration × Regime.
#   - Regression coefficients include 95% CI, R2, adjusted R2, AIC and BIC.
#   - SEM component paths reproduce the Fig. 4e structure:
#       gut16S_PC1 ~ ChemPC1
#       AI ~ gut16S_PC1 + ChemPC1
#       D_phoD ~ AI
# ============================================================

# -------------------------
# 0) Packages
# -------------------------
pkg_needed <- c(
  "readxl", "dplyr", "tibble", "purrr", "broom",
  "car", "openxlsx"
)

pkg_to_install <- pkg_needed[!sapply(pkg_needed, requireNamespace, quietly = TRUE)]
if (length(pkg_to_install) > 0) {
  install.packages(pkg_to_install)
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(broom)
  library(car)
  library(openxlsx)
})

# -------------------------
# 1) Paths
# -------------------------
out_dir <- OUTPUT_DIR
in_xlsx  <- DATA_XLSX
in_sheet <- "Master_all_data"

out_xlsx <- file.path(
  out_dir,
  "Supplementary_Table_7_Statistical_outputs_main_linear_models.xlsx"
)

if (!file.exists(in_xlsx)) {
  stop("Input file not found: ", in_xlsx)
}
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
# -------------------------
# 2) Helpers
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

fmt_formula <- function(fit) {
  paste(deparse(formula(fit)), collapse = " ")
}

get_response <- function(fit) {
  as.character(formula(fit)[[2]])
}

safe_confint <- function(fit) {
  out <- tryCatch(
    as.data.frame(confint(fit)),
    error = function(e) NULL
  )
  if (is.null(out)) {
    return(tibble(term = names(coef(fit)), CI_low = NA_real_, CI_high = NA_real_))
  }
  tibble(
    term = rownames(out),
    CI_low = out[[1]],
    CI_high = out[[2]]
  )
}

partial_r2_term <- function(fit, term) {
  sm <- summary(fit)$coefficients
  if (!term %in% rownames(sm)) return(NA_real_)
  tval <- sm[term, "t value"]
  df_res <- df.residual(fit)
  as.numeric((tval^2) / (tval^2 + df_res))
}

standardised_beta <- function(fit, term) {
  sm <- summary(fit)$coefficients
  mf <- model.frame(fit)
  y <- mf[[1]]
  
  if (!term %in% names(mf)) return(NA_real_)
  x <- mf[[term]]
  
  if (!is.numeric(x) || !is.numeric(y)) return(NA_real_)
  if (sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) return(NA_real_)
  
  b <- sm[term, "Estimate"]
  as.numeric(b * sd(x, na.rm = TRUE) / sd(y, na.rm = TRUE))
}

fit_glance <- function(fit) {
  gl <- broom::glance(fit)
  tibble(
    N = nobs(fit),
    Df_residual = df.residual(fit),
    R2 = unname(gl$r.squared),
    Adj_R2 = unname(gl$adj.r.squared),
    Sigma = unname(gl$sigma),
    F_statistic = unname(gl$statistic),
    F_num_df = unname(gl$df),
    F_den_df = unname(gl$df.residual),
    Model_P_value = unname(gl$p.value),
    AIC = AIC(fit),
    BIC = BIC(fit)
  )
}

make_anova_table <- function(fit, figure, panel, model_id, model_group, notes = NA_character_) {
  aov_tab <- as.data.frame(car::Anova(fit, type = 2))
  aov_tab$Term <- rownames(aov_tab)
  aov_tab <- aov_tab[aov_tab$Term != "Residuals", ]
    gl <- fit_glance(fit)
  
  out <- aov_tab %>%
    as_tibble() %>%
    transmute(
      Figure = figure,
      Panel = panel,
      Model_ID = model_id,
      Model_group = model_group,
      Response = get_response(fit),
      Formula = fmt_formula(fit),
      Term = Term,
      Sum_sq = .data$`Sum Sq`,
      Df = .data$Df,
      F_value = .data$`F value`,
      P_value = .data$`Pr(>F)`,
      N = gl$N,
      Df_residual = gl$Df_residual,
      R2 = gl$R2,
      Adj_R2 = gl$Adj_R2,
      AIC = gl$AIC,
      BIC = gl$BIC,
      Notes = notes
    )
  
  out
}

make_coef_table <- function(fit, figure, panel, model_id, model_group, notes = NA_character_) {
  gl <- fit_glance(fit)
  ci <- safe_confint(fit)
  
  broom::tidy(fit) %>%
    left_join(ci, by = "term") %>%
    mutate(
      Partial_R2 = purrr::map_dbl(term, ~ partial_r2_term(fit, .x)),
      Standardised_beta = purrr::map_dbl(term, ~ standardised_beta(fit, .x))
    ) %>%
    transmute(
      Figure = figure,
      Panel = panel,
      Model_ID = model_id,
      Model_group = model_group,
      Response = get_response(fit),
      Formula = fmt_formula(fit),
      Term = term,
      Estimate = estimate,
      Std_error = std.error,
      CI_low = CI_low,
      CI_high = CI_high,
      Statistic = statistic,
      Df_residual = gl$Df_residual,
      P_value = p.value,
      Partial_R2 = Partial_R2,
      Standardised_beta = Standardised_beta,
      N = gl$N,
      R2 = gl$R2,
      Adj_R2 = gl$Adj_R2,
      AIC = gl$AIC,
      BIC = gl$BIC,
      Notes = notes
    )
}

make_model_summary <- function(fit, figure, panel, model_id, model_group, notes = NA_character_) {
  gl <- fit_glance(fit)
  tibble(
    Figure = figure,
    Panel = panel,
    Model_ID = model_id,
    Model_group = model_group,
    Response = get_response(fit),
    Formula = fmt_formula(fit),
    N = gl$N,
    Df_residual = gl$Df_residual,
    R2 = gl$R2,
    Adj_R2 = gl$Adj_R2,
    Sigma = gl$Sigma,
    F_statistic = gl$F_statistic,
    F_num_df = gl$F_num_df,
    F_den_df = gl$F_den_df,
    Model_P_value = gl$Model_P_value,
    AIC = gl$AIC,
    BIC = gl$BIC,
    Notes = notes
  )
}

make_sem_paths <- function(fit, figure, panel, model_id, model_group) {
  make_coef_table(
    fit = fit,
    figure = figure,
    panel = panel,
    model_id = model_id,
    model_group = model_group,
    notes = "SEM component model; standardised beta calculated as b * sd(x) / sd(y)."
  ) %>%
    filter(Term != "(Intercept)") %>%
    transmute(
      Figure,
      Panel,
      SEM_ID = "Fig4e_piecewise_SEM",
      Component_model = Model_ID,
      Response,
      Predictor = Term,
      Formula,
      Estimate,
      Std_error,
      CI_low,
      CI_high,
      Statistic,
      Df_residual,
      P_value,
      Standardised_beta,
      N,
      R2,
      Adj_R2,
      AIC,
      BIC,
      Notes
    )
}

# -------------------------
# 3) Read and prepare data
# -------------------------
dat_raw <- readxl::read_excel(in_xlsx, sheet = in_sheet)

required_cols <- c(
  "PlotID", "Duration", "Regime",
  "ChemPC1", "Olsen-P(mg/kg)", "pH",
  "gut16S_PC1", "AI", "D_phoD"
)

missing_cols <- setdiff(required_cols, names(dat_raw))
if (length(missing_cols) > 0) {
  stop("Input table missing columns: ", paste(missing_cols, collapse = ", "))
}

dat <- dat_raw %>%
  transmute(
    PlotID = trimws(as.character(PlotID)),
    Duration = clean_duration(Duration),
    Regime = clean_regime(Regime),
    ChemPC1 = as.numeric(ChemPC1),
    Olsen_P = as.numeric(.data[["Olsen-P(mg/kg)"]]),
    pH = as.numeric(pH),
    gut16S_PC1 = as.numeric(gut16S_PC1),
    AI = as.numeric(AI),
    D_phoD = as.numeric(D_phoD)
  ) %>%
  filter(!is.na(Duration), !is.na(Regime)) %>%
  droplevels()

if (nrow(dat) == 0) {
  stop("No valid rows after cleaning Duration and Regime.")
}

# -------------------------
# 4) Fit main models
# -------------------------

# Fig. 2
fit_fig2a <- lm(ChemPC1 ~ Duration * Regime, data = dat %>% filter(is.finite(ChemPC1)))
fit_fig2b <- lm(Olsen_P ~ Duration * Regime, data = dat %>% filter(is.finite(Olsen_P)))
fit_fig2c <- lm(pH ~ Duration * Regime, data = dat %>% filter(is.finite(pH)))

# Fig. 3c
fit_fig3c <- lm(
  gut16S_PC1 ~ Duration * Regime,
  data = dat %>% filter(is.finite(gut16S_PC1))
)

# Fig. 4b
fit_fig4b <- lm(
  AI ~ Duration * Regime,
  data = dat %>% filter(is.finite(AI))
)

# Fig. 4c
fit_fig4c <- lm(
  AI ~ gut16S_PC1 + ChemPC1,
  data = dat %>% filter(is.finite(AI), is.finite(gut16S_PC1), is.finite(ChemPC1))
)

# Fig. 4d
fit_fig4d <- lm(
  D_phoD ~ AI,
  data = dat %>% filter(is.finite(D_phoD), is.finite(AI))
)

# Fig. 4e SEM components
dat_sem <- dat %>%
  filter(is.finite(ChemPC1), is.finite(gut16S_PC1), is.finite(AI), is.finite(D_phoD))

fit_sem_1 <- lm(gut16S_PC1 ~ ChemPC1, data = dat_sem)
fit_sem_2 <- lm(AI ~ gut16S_PC1 + ChemPC1, data = dat_sem)
fit_sem_3 <- lm(D_phoD ~ AI, data = dat_sem)

basis_p <- c(
  summary(lm(D_phoD ~ AI + ChemPC1, data = dat_sem))$coefficients["ChemPC1", "Pr(>|t|)"],
  summary(lm(D_phoD ~ AI + gut16S_PC1, data = dat_sem))$coefficients["gut16S_PC1", "Pr(>|t|)"]
)
basis_p <- pmax(pmin(basis_p, 1), .Machine$double.xmin)
sem_fisher <- tibble(
  Fisher.C = -2 * sum(log(basis_p)),
  df = 2 * length(basis_p),
  P.Value = pchisq(-2 * sum(log(basis_p)), df = 2 * length(basis_p), lower.tail = FALSE)
)

# -------------------------
# 5) Build output tables
# -------------------------

typeII_anova <- bind_rows(
  make_anova_table(
    fit_fig2a, "Fig. 2", "a", "Fig2a_ChemPC1_Duration_Regime",
    "Main Type-II ANOVA",
    "Soil-state background model."
  ),
  make_anova_table(
    fit_fig2b, "Fig. 2", "b", "Fig2b_OlsenP_Duration_Regime",
    "Main Type-II ANOVA",
    "Olsen-P response; original input column: Olsen-P(mg/kg)."
  ),
  make_anova_table(
    fit_fig2c, "Fig. 2", "c", "Fig2c_pH_Duration_Regime",
    "Main Type-II ANOVA",
    "Soil pH response."
  ),
  make_anova_table(
    fit_fig3c, "Fig. 3", "c", "Fig3c_gut16S_PC1_Duration_Regime",
    "Main Type-II ANOVA",
    "Gut microbial-state model."
  ),
  make_anova_table(
    fit_fig4b, "Fig. 4", "b", "Fig4b_AI_Duration_Regime",
    "Main Type-II ANOVA",
    "phoD amplification group-pattern model."
  )
)

regression_coefficients <- bind_rows(
  make_coef_table(
    fit_fig4c, "Fig. 4", "c", "Fig4c_AI_gut16S_PC1_ChemPC1",
    "Main regression",
    "Core model: gut microbial state predicts phoD amplification while accounting for soil-state background."
  ),
  make_coef_table(
    fit_fig4d, "Fig. 4", "d", "Fig4d_DphoD_AI",
    "Main regression",
    "Association between phoD amplification and soil-gut phoD community dissimilarity. Lower D_phoD indicates stronger coupling."
  ),
  make_coef_table(
    fit_sem_1, "Fig. 4", "e", "Fig4e_SEM_gut16S_PC1_ChemPC1",
    "SEM component",
    "SEM component: ChemPC1 to gut16S_PC1."
  ),
  make_coef_table(
    fit_sem_2, "Fig. 4", "e", "Fig4e_SEM_AI_gut16S_PC1_ChemPC1",
    "SEM component",
    "SEM component: gut16S_PC1 and ChemPC1 to AI."
  ),
  make_coef_table(
    fit_sem_3, "Fig. 4", "e", "Fig4e_SEM_DphoD_AI",
    "SEM component",
    "SEM component: AI to D_phoD."
  )
)

model_summary <- bind_rows(
  make_model_summary(
    fit_fig2a, "Fig. 2", "a", "Fig2a_ChemPC1_Duration_Regime",
    "Main Type-II ANOVA", "Soil-state background model."
  ),
  make_model_summary(
    fit_fig2b, "Fig. 2", "b", "Fig2b_OlsenP_Duration_Regime",
    "Main Type-II ANOVA", "Olsen-P response."
  ),
  make_model_summary(
    fit_fig2c, "Fig. 2", "c", "Fig2c_pH_Duration_Regime",
    "Main Type-II ANOVA", "Soil pH response."
  ),
  make_model_summary(
    fit_fig3c, "Fig. 3", "c", "Fig3c_gut16S_PC1_Duration_Regime",
    "Main Type-II ANOVA", "Gut microbial-state model."
  ),
  make_model_summary(
    fit_fig4b, "Fig. 4", "b", "Fig4b_AI_Duration_Regime",
    "Main Type-II ANOVA", "phoD amplification group-pattern model."
  ),
  make_model_summary(
    fit_fig4c, "Fig. 4", "c", "Fig4c_AI_gut16S_PC1_ChemPC1",
    "Main regression", "Core AI model."
  ),
  make_model_summary(
    fit_fig4d, "Fig. 4", "d", "Fig4d_DphoD_AI",
    "Main regression", "AI-D_phoD association model."
  ),
  make_model_summary(
    fit_sem_1, "Fig. 4", "e", "Fig4e_SEM_gut16S_PC1_ChemPC1",
    "SEM component", "SEM component model."
  ),
  make_model_summary(
    fit_sem_2, "Fig. 4", "e", "Fig4e_SEM_AI_gut16S_PC1_ChemPC1",
    "SEM component", "SEM component model."
  ),
  make_model_summary(
    fit_sem_3, "Fig. 4", "e", "Fig4e_SEM_DphoD_AI",
    "SEM component", "SEM component model."
  )
)

sem_paths <- bind_rows(
  make_sem_paths(
    fit_sem_1, "Fig. 4", "e", "Fig4e_SEM_gut16S_PC1_ChemPC1",
    "SEM component"
  ),
  make_sem_paths(
    fit_sem_2, "Fig. 4", "e", "Fig4e_SEM_AI_gut16S_PC1_ChemPC1",
    "SEM component"
  ),
  make_sem_paths(
    fit_sem_3, "Fig. 4", "e", "Fig4e_SEM_DphoD_AI",
    "SEM component"
  )
)

sem_fit_summary <- tibble(
  Figure = "Fig. 4",
  Panel = "e",
  SEM_ID = "Fig4e_piecewise_SEM",
  Component_models = paste(
    "gut16S_PC1 ~ ChemPC1",
    "AI ~ gut16S_PC1 + ChemPC1",
    "D_phoD ~ AI",
    sep = " | "
  ),
  Fisher_C = sem_fisher$Fisher.C[1],
  Df = sem_fisher$df[1],
  P_value = sem_fisher$P.Value[1],
  N = nrow(dat_sem),
  Notes = "Fisher's C calculated from the two basis-set tests used by Fig. 4 when piecewiseSEM is unavailable."
)

readme <- tibble(
  Field = c(
    "Table title",
    "Input file",
    "Input sheet",
    "Duration harmonisation",
    "Type-II ANOVA models",
    "Regression models",
    "SEM component models",
    "Important interpretation",
    "Excluded from this table"
  ),
  Description = c(
    "Supplementary Table 7 | Statistical outputs for main linear models",
    in_xlsx,
    in_sheet,
    "3y is harmonised to 5y; factor order is 5y, 8y, 10y.",
    "ChemPC1, Olsen-P, pH, gut16S_PC1, and AI models with Duration × Regime.",
    "AI ~ gut16S_PC1 + ChemPC1 and D_phoD ~ AI.",
    "gut16S_PC1 ~ ChemPC1; AI ~ gut16S_PC1 + ChemPC1; D_phoD ~ AI.",
    "Lower D_phoD indicates stronger soil-gut phoD community coupling.",
    "PERMANOVA, dispersion tests, chronosequence robustness, metric robustness, distance robustness, and P-index robustness should be reported in separate supplementary tables."
  )
)

# -------------------------
# 6) Write Excel workbook
# -------------------------
wb <- createWorkbook()

addWorksheet(wb, "README")
writeData(wb, "README", readme)

addWorksheet(wb, "TypeII_ANOVA")
writeData(wb, "TypeII_ANOVA", typeII_anova)

addWorksheet(wb, "Regression_coefficients")
writeData(wb, "Regression_coefficients", regression_coefficients)

addWorksheet(wb, "Model_summary")
writeData(wb, "Model_summary", model_summary)

addWorksheet(wb, "SEM_paths")
writeData(wb, "SEM_paths", sem_paths)

addWorksheet(wb, "SEM_fit")
writeData(wb, "SEM_fit", sem_fit_summary)

# Freeze first row and set simple widths
for (sh in names(wb)) {
  freezePane(wb, sh, firstRow = TRUE)
  setColWidths(wb, sh, cols = 1:60, widths = "auto")
}

# Numeric style
num_style <- createStyle(numFmt = "0.0000")

sheet_objects <- list(
  TypeII_ANOVA = typeII_anova,
  Regression_coefficients = regression_coefficients,
  Model_summary = model_summary,
  SEM_paths = sem_paths,
  SEM_fit = sem_fit_summary
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

cat("\nSupplementary Table 7 generated successfully:\n")
cat(out_xlsx, "\n\n")

cat("Sheets written:\n")
cat("- README\n")
cat("- TypeII_ANOVA\n")
cat("- Regression_coefficients\n")
cat("- Model_summary\n")
cat("- SEM_paths\n")
cat("- SEM_fit\n")
