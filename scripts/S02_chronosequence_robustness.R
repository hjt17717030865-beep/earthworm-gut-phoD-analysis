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
# Figure S2: Chronosequence robustness
#
# S2A  Within-duration contrasts:
#      AI ~ Regime within each Duration
#      Plot NPK - CK and NPKOM - CK
#
# S2B  Leave-one-duration-out robustness:
#      AI ~ gut16S_PC1 + ChemPC1
#
# S2C  Duration-centred model:
#      AI_c ~ gut16S_PC1_c + ChemPC1_c
#      Plot coefficients only (NOT scatter regression)
#
# INPUT:
#   ???????????????.xlsx (Sheet1)
#
# REQUIRED COLUMNS:
#   PlotID, Duration, Regime, AI, gut16S_PC1, ChemPC1
#
# OUTPUT:
#   FigS2_chronosequence_robustness_v2.svg
#   FigS2_chronosequence_robustness_v2.png
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(car)
  library(svglite)
  library(tibble)
  library(purrr)
  library(tidyr)
  library(broom)
})

# -------------------------
# -------------------------
# 0) Paths
# -------------------------
out_dir <- OUTPUT_DIR
in_xlsx  <- DATA_XLSX

# -------------------------
# 1) Load table
# -------------------------
dat <- read_excel(in_xlsx, sheet = "Master_all_data")
need <- c("PlotID", "Duration", "Regime", "AI", "gut16S_PC1", "ChemPC1")
miss <- setdiff(need, names(dat))
if (length(miss) > 0) stop("Input table missing: ", paste(miss, collapse = ", "))

dat <- dat %>%
  mutate(
    PlotID     = trimws(as.character(PlotID)),
    Duration   = factor(as.character(Duration), levels = c("5y", "8y", "10y")),
    Regime     = factor(as.character(Regime),   levels = c("CK", "NPK", "NPKOM")),
    AI         = as.numeric(AI),
    gut16S_PC1 = as.numeric(gut16S_PC1),
    ChemPC1    = as.numeric(ChemPC1)
  ) %>%
  filter(
    is.finite(AI),
    is.finite(gut16S_PC1),
    is.finite(ChemPC1)
  )

# -------------------------
# 2) Style helpers
# -------------------------
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
      panel.spacing = unit(0.45, "lines"),
      axis.ticks.length = unit(2, "pt"),
      plot.margin = margin(3, 4, 3, 3, unit = "mm")
    )
}

pal_regime <- c(
  "CK"    = "#B8B8B8",
  "NPK"   = "#7FAFB4",
  "NPKOM" = "#D39A92"
)

# ============================================================
# 3) S2A: Within-duration AI ~ Regime
#    show contrasts relative to CK
# ============================================================

within_res <- dat %>%
  group_by(Duration) %>%
  group_modify(~{
    dd <- .x
    
    fit <- lm(AI ~ Regime, data = dd)
    sm  <- summary(fit)
    an  <- anova(fit)
    
    broom::tidy(fit, conf.int = TRUE) %>%
      filter(term %in% c("RegimeNPK", "RegimeNPKOM")) %>%
      mutate(
        contrast = dplyr::case_when(
          term == "RegimeNPK"   ~ "NPK - CK",
          term == "RegimeNPKOM" ~ "NPKOM - CK",
          TRUE ~ NA_character_
        ),
        model_p = an["Regime", "Pr(>F)"],
        r2 = sm$r.squared,
        n = nrow(dd)
      )
  }) %>%
  ungroup() %>%
  mutate(
    Duration = factor(Duration, levels = c("5y", "8y", "10y")),
    contrast = factor(contrast, levels = c("NPK - CK", "NPKOM - CK"))
  )

lab_y_A <- max(within_res$conf.high, na.rm = TRUE) +
  diff(range(c(within_res$conf.low, within_res$conf.high), na.rm = TRUE)) * 0.28

within_lab <- within_res %>%
  distinct(Duration, model_p, r2, n) %>%
  mutate(
    x = 0.58,
    y = lab_y_A,
    label = paste0(
      "AI ~ Regime\n",
      "R2 = ", fmt_num(r2), "\n",
      "P = ", fmt_p(model_p)
    )
  )

rngA <- range(c(within_res$conf.low, within_res$conf.high, 0), na.rm = TRUE)
padA <- diff(rngA) * 0.30
if (!is.finite(padA) || padA == 0) padA <- 0.03

pS2A <- ggplot(within_res, aes(x = contrast, y = estimate)) +
  geom_hline(yintercept = 0, linewidth = 0.30, colour = "grey60") +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width = 0.08,
    linewidth = 0.36,
    colour = "grey25"
  ) +
  geom_point(
    aes(fill = contrast),
    shape = 21,
    size = 2.2,
    stroke = 0.30,
    colour = "grey20",
    alpha = 0.90
  ) +
  geom_label(
    data = within_lab,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.45,
    hjust = 0,
    vjust = 1,
    label.size = NA,
    fill = "white",
    label.padding = unit(0.10, "lines"),
    colour = "grey20"
  ) +
  facet_grid(. ~ Duration) +
  scale_fill_manual(
    values = c("NPK - CK" = pal_regime["NPK"], "NPKOM - CK" = pal_regime["NPKOM"]),
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = expression(Delta * "AI relative to CK")
  ) +
  coord_cartesian(
    ylim = c(rngA[1] - padA * 0.35, lab_y_A + padA * 0.10),
    clip = "off"
  ) +
  theme_nature(base_size = 9) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 8.5)
  )

# ============================================================
# 4) S2B: Leave-one-duration-out
#    AI ~ gut16S_PC1 + ChemPC1
# ============================================================

loo_levels <- levels(dat$Duration)

loo_res <- map_dfr(loo_levels, function(d0) {
  dd <- dat %>% filter(Duration != d0)
  
  fit <- lm(AI ~ gut16S_PC1 + ChemPC1, data = dd)
  sm  <- summary(fit)
  
  broom::tidy(fit, conf.int = TRUE) %>%
    filter(term %in% c("gut16S_PC1", "ChemPC1")) %>%
    mutate(
      left_out = d0,
      n = nrow(dd),
      r2 = sm$r.squared
    )
}) %>%
  mutate(
    left_out = factor(left_out, levels = c("5y", "8y", "10y")),
    term = factor(term, levels = c("gut16S_PC1", "ChemPC1"))
  )

lab_y_B <- max(loo_res$conf.high, na.rm = TRUE) +
  diff(range(c(loo_res$conf.low, loo_res$conf.high), na.rm = TRUE)) * 0.28

loo_lab <- loo_res %>%
  group_by(left_out) %>%
  summarise(
    r2 = unique(r2),
    n = unique(n),
    gut_p = p.value[term == "gut16S_PC1"],
    chem_p = p.value[term == "ChemPC1"],
    .groups = "drop"
  ) %>%
  mutate(
    x = 0.58,
    y = lab_y_B,
    label = paste0(
      "R2 = ", fmt_num(r2), "\n",
      "gut P = ", fmt_p(gut_p), "\n",
      "Chem P = ", fmt_p(chem_p)
    )
  )

rngB <- range(c(loo_res$conf.low, loo_res$conf.high, 0), na.rm = TRUE)
padB <- diff(rngB) * 0.30
if (!is.finite(padB) || padB == 0) padB <- 0.02

pS2B <- ggplot(loo_res, aes(x = term, y = estimate)) +
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
    size = 2.2,
    stroke = 0.30,
    colour = "grey20",
    alpha = 0.90
  ) +
  geom_label(
    data = loo_lab,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.45,
    hjust = 0,
    vjust = 1,
    label.size = NA,
    fill = "white",
    label.padding = unit(0.10, "lines"),
    colour = "grey20"
  ) +
  facet_grid(. ~ left_out, labeller = labeller(left_out = function(x) paste0("Leave ", x, " out"))) +
  scale_fill_manual(
    values = c("gut16S_PC1" = "grey75", "ChemPC1" = "grey55"),
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Model coefficient"
  ) +
  coord_cartesian(
    ylim = c(rngB[1] - padB * 0.35, lab_y_B + padB * 0.10),
    clip = "off"
  ) +
  theme_nature(base_size = 9) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 8.5)
  )

# ============================================================
# 5) S2C: Duration-centred model
#    AI_c ~ gut16S_PC1_c + ChemPC1_c
#    coefficient plot only
# ============================================================

dat_c <- dat %>%
  group_by(Duration) %>%
  mutate(
    AI_c         = AI - mean(AI, na.rm = TRUE),
    gut16S_PC1_c = gut16S_PC1 - mean(gut16S_PC1, na.rm = TRUE),
    ChemPC1_c    = ChemPC1 - mean(ChemPC1, na.rm = TRUE)
  ) %>%
  ungroup()

fit_c <- lm(AI_c ~ gut16S_PC1_c + ChemPC1_c, data = dat_c)
sm_c  <- summary(fit_c)

coef_c <- broom::tidy(fit_c, conf.int = TRUE) %>%
  filter(term %in% c("gut16S_PC1_c", "ChemPC1_c")) %>%
  mutate(
    term = factor(term, levels = c("gut16S_PC1_c", "ChemPC1_c"))
  )

lab_C <- tibble(
  x = 0.45,
  y = max(coef_c$conf.high, na.rm = TRUE) +
    diff(range(c(coef_c$conf.low, coef_c$conf.high), na.rm = TRUE)) * 0.28,
  label = paste0(
    "Duration-centred model\n",
    "AI_c ~ gut16S_PC1_c + ChemPC1_c\n",
    "R2 = ", fmt_num(sm_c$r.squared), "\n",
    "gut P = ", fmt_p(coef_c$p.value[coef_c$term == "gut16S_PC1_c"]), "\n",
    "Chem P = ", fmt_p(coef_c$p.value[coef_c$term == "ChemPC1_c"])
  )
)

rngC <- range(c(coef_c$conf.low, coef_c$conf.high, 0), na.rm = TRUE)
padC <- diff(rngC) * 0.35
if (!is.finite(padC) || padC == 0) padC <- 0.02

pS2C <- ggplot(coef_c, aes(x = term, y = estimate)) +
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
    size = 2.4,
    stroke = 0.30,
    colour = "grey20",
    alpha = 0.92
  ) +
  geom_text(
    data = lab_C,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.4,
    hjust = 0,
    vjust = 1,
    colour = "grey20",
    lineheight = 1.15
  ) +
  scale_fill_manual(
    values = c("gut16S_PC1_c" = "grey75", "ChemPC1_c" = "grey55"),
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Duration-centred coefficient"
  ) +
  coord_cartesian(
    ylim = c(rngC[1] - padC * 0.35, lab_C$y + padC * 0.30),
    clip = "off"
  ) +
  theme_nature(base_size = 9) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 8.5)
  )

# ============================================================
# 6) Combine & export
# ============================================================

figS2 <- pS2A / pS2B / pS2C +
  plot_layout(heights = c(0.95, 0.95, 0.85)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 10, colour = "grey20")
  )

print(figS2)

ggsave(
  filename = file.path(out_dir, "FigS2_chronosequence_robustness.svg"),
  plot = figS2,
  width = 170,
  height = 168,
  units = "mm",
  device = svglite::svglite
)

ggsave(
  filename = file.path(out_dir, "FigS2_chronosequence_robustness.png"),
  plot = figS2,
  width = 170,
  height = 168,
  units = "mm",
  dpi = 600
)

cat("Saved Figure S2 to: ", out_dir, "\n")
