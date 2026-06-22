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
# Figure S8 | Model diagnostics
# Row-wise panel labels: a-e = models; columns = diagnostics
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(car)
})

# ------------------------------------------------------------
# ------------------------------------------------------------
# 1) Paths
# ------------------------------------------------------------
out_dir <- OUTPUT_DIR
in_xlsx  <- DATA_XLSX
in_sheet <- "Master_all_data"

if (!file.exists(in_xlsx)) {
  stop("Input file not found: ", in_xlsx)
}

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
# ------------------------------------------------------------
# 2) Read and prepare data
# ------------------------------------------------------------
dat_raw <- readxl::read_excel(in_xlsx, sheet = in_sheet)

need <- c(
  "PlotID", "Duration", "Regime",
  "ChemPC1", "gut16S_PC1", "AI", "D_phoD"
)

miss <- setdiff(need, names(dat_raw))
if (length(miss) > 0) {
  stop("Input table missing columns: ", paste(miss, collapse = ", "))
}

dat <- dat_raw %>%
  transmute(
    PlotID = trimws(as.character(PlotID)),
    Duration = as.character(Duration),
    Regime = as.character(Regime),
    ChemPC1 = as.numeric(ChemPC1),
    gut16S_PC1 = as.numeric(gut16S_PC1),
    AI = as.numeric(AI),
    D_phoD = as.numeric(D_phoD)
  ) %>%
  mutate(
    Duration = trimws(Duration),
    Duration = gsub("\\s+", "", Duration),
    Duration = gsub("years|Years|yr|YR|Year|year", "y", Duration),
    Duration = gsub("^5$", "5y", Duration),
    Duration = gsub("^8$", "8y", Duration),
    Duration = gsub("^10$", "10y", Duration),
    Duration = ifelse(Duration == "3y", "5y", Duration),
    Duration = factor(Duration, levels = c("5y", "8y", "10y")),
    Regime = trimws(Regime),
    Regime = factor(Regime, levels = c("CK", "NPK", "NPKOM"))
  ) %>%
  filter(!is.na(Duration), !is.na(Regime))

# ------------------------------------------------------------
# 3) Fit current main-line models
# ------------------------------------------------------------
fit_soil_state <- lm(
  ChemPC1 ~ Duration * Regime,
  data = dat %>% filter(is.finite(ChemPC1))
)

fit_gut_state <- lm(
  gut16S_PC1 ~ Duration * Regime,
  data = dat %>% filter(is.finite(gut16S_PC1))
)

fit_AI_group <- lm(
  AI ~ Duration * Regime,
  data = dat %>% filter(is.finite(AI))
)

fit_AI_mechanism <- lm(
  AI ~ gut16S_PC1 + ChemPC1,
  data = dat %>% filter(is.finite(AI), is.finite(gut16S_PC1), is.finite(ChemPC1))
)

fit_DphoD <- lm(
  D_phoD ~ AI,
  data = dat %>% filter(is.finite(D_phoD), is.finite(AI))
)

fits <- list(
  "a  ChemPC1 ~ Duration × Regime" = fit_soil_state,
  "b  gut16S_PC1 ~ Duration × Regime" = fit_gut_state,
  "c  AI ~ Duration × Regime" = fit_AI_group,
  "d  AI ~ gut16S_PC1 + ChemPC1" = fit_AI_mechanism,
  "e  D_phoD ~ AI" = fit_DphoD
)

# ------------------------------------------------------------
# 4) Theme
# ------------------------------------------------------------
theme_s8 <- function(base_size = 8.2) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size, colour = "grey20"),
      axis.text  = element_text(size = base_size - 0.4, colour = "grey25"),
      axis.line  = element_line(linewidth = 0.30, colour = "grey25"),
      axis.ticks = element_line(linewidth = 0.30, colour = "grey25"),
      plot.title = element_text(size = base_size, face = "bold", colour = "grey20", hjust = 0),
      plot.margin = margin(1.4, 1.4, 1.4, 1.4, unit = "mm")
    )
}

# ------------------------------------------------------------
# 5) Diagnostics data
# ------------------------------------------------------------
make_diag_df <- function(fit) {
  data.frame(
    fitted = fitted(fit),
    residual = resid(fit),
    stdresid = rstandard(fit),
    sqrt_abs_stdresid = sqrt(abs(rstandard(fit))),
    cooks = cooks.distance(fit),
    obs = seq_len(nobs(fit))
  )
}

# ------------------------------------------------------------
# 6) Plot builders
# ------------------------------------------------------------
p_resid <- function(df, title = "") {
  ggplot(df, aes(x = fitted, y = residual)) +
    geom_hline(yintercept = 0, linewidth = 0.30, colour = "grey60") +
    geom_point(size = 1.15, alpha = 0.85, colour = "grey35") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.35, colour = "grey20") +
    labs(title = title, x = NULL, y = "Residuals") +
    theme_s8()
}

p_qq <- function(df, title = "") {
  qq <- qqnorm(df$stdresid, plot.it = FALSE)
  qq_df <- data.frame(theoretical = qq$x, sample = qq$y)
  
  ggplot(qq_df, aes(x = theoretical, y = sample)) +
    geom_point(size = 1.15, alpha = 0.85, colour = "grey35") +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.30,
                linetype = "dashed", colour = "grey60") +
    labs(title = title, x = NULL, y = "Std residuals") +
    theme_s8()
}

p_scale <- function(df, title = "") {
  ggplot(df, aes(x = fitted, y = sqrt_abs_stdresid)) +
    geom_point(size = 1.15, alpha = 0.85, colour = "grey35") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.35, colour = "grey20") +
    labs(title = title, x = NULL, y = expression(sqrt("|Std residuals|"))) +
    theme_s8()
}

p_cook <- function(df, title = "") {
  n <- nrow(df)
  cook_cut <- 4 / n
  
  ggplot(df, aes(x = obs, y = cooks)) +
    geom_col(width = 0.70, fill = "grey78", colour = "grey35", linewidth = 0.18) +
    geom_hline(yintercept = cook_cut, linewidth = 0.30,
               colour = "grey55", linetype = "dashed") +
    labs(title = title, x = "Obs", y = "Cook's D") +
    theme_s8()
}

# ------------------------------------------------------------
# 7) Build figure by rows
# ------------------------------------------------------------
model_names <- names(fits)

rows <- lapply(seq_along(fits), function(i) {
  df <- make_diag_df(fits[[i]])
  
  p_resid(df, model_names[i]) |
    p_qq(df, if (i == 1) "Normal Q–Q" else "") |
    p_scale(df, if (i == 1) "Scale-location" else "") |
    p_cook(df, if (i == 1) "Cook's distance" else "")
})

fig_s8 <- wrap_plots(rows, ncol = 1) +
  plot_layout(heights = rep(1, length(rows)))

print(fig_s8)

# ------------------------------------------------------------
# 8) Save figure and model summaries
# ------------------------------------------------------------
ggsave(
  filename = file.path(out_dir, "Figure_S8_model_diagnostics.svg"),
  plot = fig_s8,
  width = 210,
  height = 185,
  units = "mm",
  device = svglite::svglite,
  bg = "white"
)

ggsave(
  filename = file.path(out_dir, "Figure_S8_model_diagnostics.png"),
  plot = fig_s8,
  width = 210,
  height = 185,
  units = "mm",
  dpi = 600,
  bg = "white"
)

sink(file.path(out_dir, "Figure_S8_model_diagnostics_model_summaries.txt"))
cat("Figure S8 model diagnostics: fitted model summaries\n")
cat("===================================================\n\n")
for (nm in names(fits)) {
  cat(nm, "\n")
  cat(strrep("-", nchar(nm)), "\n", sep = "")
  print(summary(fits[[nm]]))
  cat("\nType-II ANOVA where applicable:\n")
  suppressWarnings(print(car::Anova(fits[[nm]], type = 2)))
  cat("\n\n")
}
sink()

cat("Saved Figure S8 diagnostics to: ", out_dir, "\n")
