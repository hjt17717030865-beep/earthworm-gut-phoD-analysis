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
# Figure S3 | qPCR technical reliability
#
# INPUT:
#   Step1_qPCR_processed.xlsx
#
# OUTPUT:
#   FigS3_qPCR_technical_reliability.svg
#   FigS3_qPCR_technical_reliability.png
#
# PANELS:
#   S3A | technical replicate consistency for 16S / phoD
#   S3B | replicate variability summary for 16S / phoD
#
# NOTE:
#   16S has 2 technical replicates
#   phoD has 3 technical replicates
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(tibble)
})

# -------------------------
# -------------------------
# 0) I/O paths
# -------------------------
out_dir <- OUTPUT_DIR
out_svg  <- file.path(out_dir, "FigS3_qPCR_technical_reliability.svg")
out_png  <- file.path(out_dir, "FigS3_qPCR_technical_reliability.png")

# -------------------------
# 1) Load data
# -------------------------
df <- readxl::read_excel(
  DATA_XLSX,
  sheet = "qPCR_16S_phoD_replicates"
) %>%
  mutate(
    Duration = factor(Duration, levels = c("5y", "8y", "10y")),
    Regime   = factor(Regime,   levels = c("CK", "NPK", "NPKOM"))
  )
# -------------------------
# 2) Theme
# -------------------------
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
      panel.spacing = unit(0.40, "lines"),
      axis.ticks.length = unit(2, "pt"),
      plot.margin = margin(3, 4, 3, 3, unit = "mm")
    )
}

fmt_r <- function(x) {
  ifelse(is.na(x), "NA", sprintf("%.3f", x))
}

fmt_n <- function(x) {
  ifelse(is.na(x), "NA", as.character(x))
}

# -------------------------
# 3) Helper functions
# -------------------------
safe_log10 <- function(x) {
  ifelse(is.finite(x) & x > 0, log10(x), NA_real_)
}

make_pair_df3 <- function(data, rep1, rep2, rep3, target_label) {
  d <- data %>%
    transmute(
      rep1 = safe_log10(.data[[rep1]]),
      rep2 = safe_log10(.data[[rep2]]),
      rep3 = safe_log10(.data[[rep3]])
    )
  
  bind_rows(
    d %>% transmute(x = rep1, y = rep2),
    d %>% transmute(x = rep1, y = rep3),
    d %>% transmute(x = rep2, y = rep3)
  ) %>%
    filter(is.finite(x), is.finite(y)) %>%
    mutate(Target = target_label)
}

pair_stats <- function(pair_df) {
  pair_df %>%
    group_by(Target) %>%
    summarise(
      r = suppressWarnings(cor(x, y, method = "pearson", use = "complete.obs")),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      x = -Inf,
      y = Inf,
      label = paste0(
        "r = ", fmt_r(r),
        "\n",
        "n = ", fmt_n(n)
      )
    )
}

calc_cv2 <- function(a, b) {
  mat <- cbind(a, b)
  apply(mat, 1, function(v) {
    v <- v[is.finite(v) & v > 0]
    if (length(v) < 2) return(NA_real_)
    m <- mean(v)
    s <- sd(v)
    if (!is.finite(m) || m <= 0) return(NA_real_)
    100 * s / m
  })
}

calc_cv3 <- function(a, b, c) {
  mat <- cbind(a, b, c)
  apply(mat, 1, function(v) {
    v <- v[is.finite(v) & v > 0]
    if (length(v) < 2) return(NA_real_)
    m <- mean(v)
    s <- sd(v)
    if (!is.finite(m) || m <= 0) return(NA_real_)
    100 * s / m
  })
}

# -------------------------
# 4) Panel S3A:
#    technical replicate consistency
# -------------------------

pair_16S <- df %>%
  transmute(
    x = safe_log10(`16S_rep1`),
    y = safe_log10(`16S_rep2`)
  ) %>%
  filter(is.finite(x), is.finite(y)) %>%
  mutate(Target = "16S")

pair_phoD <- make_pair_df3(df, "phoD_rep1", "phoD_rep2", "phoD_rep3", "phoD")

pair_all <- bind_rows(pair_16S, pair_phoD) %>%
  mutate(
    Target = factor(Target, levels = c("16S", "phoD"))
  )

pair_lab <- pair_stats(pair_all) %>%
  mutate(
    Target = factor(Target, levels = c("16S", "phoD"))
  )

p_S3A <- ggplot(pair_all, aes(x = x, y = y)) +
  geom_abline(
    slope = 1, intercept = 0,
    linetype = "dashed", linewidth = 0.35, colour = "grey55"
  ) +
  geom_point(
    size = 1.25, alpha = 0.75, colour = "grey35"
  ) +
  geom_text(
    data = pair_lab,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.15,
    size = 3.0, colour = "grey20", lineheight = 1.10
  ) +
  facet_grid(. ~ Target, scales = "free") +
  labs(
    x = expression(log[10] * "(copy number, technical replicate)"),
    y = expression(log[10] * "(copy number, technical replicate)")
  ) +
  coord_cartesian(clip = "off") +
  theme_nature(base_size = 9)

# -------------------------
# 5) Panel S3B:
#    replicate variability summary
# -------------------------
cv_df <- df %>%
  transmute(
    CV_16S  = calc_cv2(`16S_rep1`, `16S_rep2`),
    CV_phoD = calc_cv3(phoD_rep1, phoD_rep2, phoD_rep3)
  ) %>%
  pivot_longer(
    cols = starts_with("CV_"),
    names_to = "Target",
    values_to = "CV"
  ) %>%
  mutate(
    Target = factor(
      Target,
      levels = c("CV_16S", "CV_phoD"),
      labels = c("16S", "phoD")
    )
  ) %>%
  filter(is.finite(CV))

cv_lab <- cv_df %>%
  group_by(Target) %>%
  summarise(
    med = median(CV, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    x = Target,
    y = Inf,
    label = paste0(
      "median CV = ", sprintf("%.1f", med), "%\n",
      "n = ", n
    )
  )

p_S3B <- ggplot(cv_df, aes(x = Target, y = CV)) +
  geom_boxplot(
    width = 0.55,
    linewidth = 0.32,
    colour = "grey25",
    fill = "grey90",
    outlier.shape = NA
  ) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    size = 0.95, alpha = 0.72, colour = "grey45"
  ) +
  geom_text(
    data = cv_lab,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0.5, vjust = 1.10,
    size = 3.0, colour = "grey20", lineheight = 1.10
  ) +
  labs(
    x = NULL,
    y = "Technical replicate CV (%)"
  ) +
  coord_cartesian(clip = "off") +
  theme_nature(base_size = 9)

# -------------------------
# 6) Assemble
# -------------------------
fig_S3 <- p_S3A / p_S3B +
  plot_layout(heights = c(1.05, 0.95)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 10, colour = "grey20")
  )

print(fig_S3)

# -------------------------
# 7) Export
# -------------------------
ggsave(
  filename = out_svg,
  plot = fig_S3,
  width = 170, height = 125, units = "mm",
  device = svglite::svglite
)

ggsave(
  filename = out_png,
  plot = fig_S3,
  width = 170, height = 125, units = "mm",
  dpi = 600, bg = "white"
)

message("Saved SVG: ", out_svg)
message("Saved PNG: ", out_png)
