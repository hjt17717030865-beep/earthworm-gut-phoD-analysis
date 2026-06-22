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
# Figure S4 | Metric robustness of phoD amplification
#
# PURPOSE:
#   Compare the same core model across three alternative responses:
#   1) AI
#   2) d_phoD = logphoD_gut - logphoD_soil
#   3) d_16S  = log16S_gut  - log16S_soil
#
# MODEL:
#   response ~ gut16S_PC1 + ChemPC1
#
# INPUT:
#
# REQUIRED BASIC COLUMNS:
#   PlotID, Duration, Regime, AI, gut16S_PC1, ChemPC1
#
# REQUIRED COMPONENT COLUMNS (first existing will be used):
#   gut phoD : logphoD_gut   > gut_logphoD
#   soil phoD: logphoD_soil  > soil_logphoD
#   gut 16S  : log16S_gut    > gut_log16S
#   soil 16S : log16S_soil   > soil_log16S
#
# OUTPUT:
#   FigS4_metric_robustness.svg
#   FigS4_metric_robustness.png
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(broom)
  library(svglite)
  library(tibble)
})

# -------------------------
# -------------------------
# 0) Paths
# -------------------------
in_dir <- CODE_ROOT
out_dir <- OUTPUT_DIR
infile   <- DATA_XLSX
out_svg  <- file.path(out_dir, "FigS4_metric_robustness.svg")
out_png  <- file.path(out_dir, "FigS4_metric_robustness.png")
in_sheet <- "Master_all_data"
# -------------------------
# 1) Helper functions
# -------------------------
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
    is.na(p), "NA",
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

fmt_num <- function(x, digits = 3) {
  sprintf(paste0("%.", digits, "f"), x)
}

theme_nature <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size, colour = "grey20"),
      axis.text  = element_text(size = base_size, colour = "grey25"),
      axis.line  = element_line(linewidth = 0.30, colour = "grey25"),
      axis.ticks = element_line(linewidth = 0.30, colour = "grey25"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, colour = "grey20"),
      legend.position = "none",
      panel.border = element_blank(),
      panel.spacing = unit(0.42, "lines"),
      axis.ticks.length = unit(2, "pt"),
      plot.margin = margin(3, 4, 3, 3, unit = "mm")
    )
}

# -------------------------
# 2) Load data
# -------------------------
dat0 <- read_excel(infile, sheet = in_sheet)

need_basic <- c("PlotID", "Duration", "Regime", "AI", "gut16S_PC1", "ChemPC1")
miss_basic <- setdiff(need_basic, names(dat0))
if (length(miss_basic) > 0) {
  stop("Input table missing: ", paste(miss_basic, collapse = ", "))
}

col_gut_phoD  <- pick_first_existing(dat0, c("logphoD_gut",  "gut_logphoD"),  "gut phoD")
col_soil_phoD <- pick_first_existing(dat0, c("logphoD_soil", "soil_logphoD"), "soil phoD")
col_gut_16S   <- pick_first_existing(dat0, c("log16S_gut",   "gut_log16S"),   "gut 16S")
col_soil_16S  <- pick_first_existing(dat0, c("log16S_soil",  "soil_log16S"),  "soil 16S")

dat <- dat0 %>%
  transmute(
    PlotID      = trimws(as.character(PlotID)),
    Duration    = factor(as.character(Duration), levels = c("5y", "8y", "10y")),
    Regime      = factor(as.character(Regime),   levels = c("CK", "NPK", "NPKOM")),
    AI          = as.numeric(AI),
    gut16S_PC1  = as.numeric(gut16S_PC1),
    ChemPC1     = as.numeric(ChemPC1),
    gut_phoD    = as.numeric(.data[[col_gut_phoD]]),
    soil_phoD   = as.numeric(.data[[col_soil_phoD]]),
    gut_16S     = as.numeric(.data[[col_gut_16S]]),
    soil_16S    = as.numeric(.data[[col_soil_16S]])
  ) %>%
  mutate(
    d_phoD = gut_phoD - soil_phoD,
    d_16S  = gut_16S  - soil_16S
  ) %>%
  filter(
    is.finite(AI),
    is.finite(gut16S_PC1),
    is.finite(ChemPC1),
    is.finite(d_phoD),
    is.finite(d_16S)
  )

# -------------------------
# 3) Fit same model across responses
# -------------------------
fit_one_response <- function(df, response_var) {
  fml <- as.formula(paste(response_var, "~ gut16S_PC1 + ChemPC1"))
  fit <- lm(fml, data = df)
  sm  <- summary(fit)
  
  broom::tidy(fit, conf.int = TRUE) %>%
    filter(term %in% c("gut16S_PC1", "ChemPC1")) %>%
    mutate(
      response = response_var,
      model_r2 = sm$r.squared,
      n = nrow(df)
    )
}

coef_df <- bind_rows(
  fit_one_response(dat, "AI"),
  fit_one_response(dat, "d_phoD"),
  fit_one_response(dat, "d_16S")
) %>%
  mutate(
    response_lab = factor(
      response,
      levels = c("AI", "d_phoD", "d_16S"),
      labels = c("AI", "Δlog10 phoD", "Δlog10 16S")
    ),
    term = factor(
      term,
      levels = c("gut16S_PC1", "ChemPC1"),
      labels = c("gut16S_PC1", "ChemPC1")
    )
  )

# -------------------------
# 4) Panel labels
# -------------------------
lab_df <- coef_df %>%
  group_by(response_lab) %>%
  summarise(
    r2 = unique(model_r2),
    gut_p = p.value[term == "gut16S_PC1"],
    chem_p = p.value[term == "ChemPC1"],
    .groups = "drop"
  )

rng <- range(c(coef_df$conf.low, coef_df$conf.high, 0), na.rm = TRUE)
pad <- diff(rng) * 0.30
if (!is.finite(pad) || pad == 0) pad <- 0.03

label_y <- max(coef_df$conf.high, na.rm = TRUE) + diff(rng) * 0.24

lab_df <- lab_df %>%
  mutate(
    x = 0.58,
    y = label_y,
    label = paste0(
      "R² = ", fmt_num(r2), "\n",
      "gut P = ", fmt_p(gut_p), "\n",
      "Chem P = ", fmt_p(chem_p)
    )
  )

# -------------------------
# 5) Plot
# -------------------------
p_S4 <- ggplot(coef_df, aes(x = term, y = estimate)) +
  geom_hline(yintercept = 0, linewidth = 0.30, colour = "grey60") +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width = 0.08,
    linewidth = 0.36,
    colour = "grey25"
  ) +
  geom_point(
    aes(fill = term),
    shape = 21,
    size = 2.3,
    stroke = 0.30,
    colour = "grey20",
    alpha = 0.92
  ) +
  geom_text(
    data = lab_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.0,
    colour = "grey20",
    lineheight = 1.10
  ) +
  facet_grid(. ~ response_lab) +
  scale_fill_manual(
    values = c("gut16S_PC1" = "grey75", "ChemPC1" = "grey55"),
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Model coefficient"
  ) +
  coord_cartesian(
    ylim = c(rng[1] - pad * 0.35, label_y + pad * 0.12),
    clip = "off"
  ) +
  theme_nature(base_size = 9) +
  theme(
    axis.text.x = element_text(size = 8.5)
  )

print(p_S4)

# -------------------------
# 6) Export
# -------------------------
ggsave(
  filename = out_svg,
  plot = p_S4,
  width = 170, height = 62, units = "mm",
  device = svglite::svglite
)

ggsave(
  filename = out_png,
  plot = p_S4,
  width = 170, height = 62, units = "mm",
  dpi = 600, bg = "white"
)

message("Saved SVG: ", out_svg)
message("Saved PNG: ", out_png)
