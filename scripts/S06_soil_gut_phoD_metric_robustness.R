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

# =========================
# S7 | Distance metric robustness of D_phoD
# =========================

# ---------- 0. packages ----------
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(purrr)
library(patchwork)

# ---------- 1. file path ----------
sum_file   <- DATA_XLSX
sheet_name <- "Master_all_data"
out_dir <- OUTPUT_DIR
# ---------- 2. read summary table ----------
# ---------- 2. read summary table ----------
sum_df <- readxl::read_excel(sum_file, sheet = sheet_name) %>%
  as.data.frame()
# ---------- 3. keep needed columns ----------
# ---------- 3. keep needed columns ----------
plot_df <- sum_df %>%
  dplyr::select(PlotID, Duration, Regime, AI, D_phoD, D_phoD_BC) %>%
  dplyr::rename(D_phoD_Aitchison = D_phoD) %>%
  dplyr::mutate(
    AI               = as.numeric(AI),
    D_phoD_Aitchison = as.numeric(D_phoD_Aitchison),
    D_phoD_BC       = as.numeric(D_phoD_BC),
    Duration = factor(Duration, levels = c("5y", "8y", "10y")),
    Regime   = factor(Regime,   levels = c("CK", "NPK", "NPKOM"))
  ) %>%
  dplyr::filter(!is.na(AI), !is.na(D_phoD_Aitchison), !is.na(D_phoD_BC))

# ---------- 4. long format ----------
dat_long <- tidyr::pivot_longer(
  data = plot_df,
  cols = c("D_phoD_Aitchison", "D_phoD_BC"),
  names_to = "Metric",
  values_to = "D_phoD"
) %>%
  dplyr::mutate(
    Metric = dplyr::case_when(
      Metric == "D_phoD_Aitchison" ~ "Aitchison",
      Metric == "D_phoD_BC"       ~ "Bray-Curtis",
      TRUE ~ Metric
    ),
    Metric = factor(Metric, levels = c("Aitchison", "Bray-Curtis"))
  )

# ---------- 5. fit models ----------
model_tbl <- dat_long %>%
  dplyr::group_by(Metric) %>%
  tidyr::nest() %>%
  dplyr::mutate(
    model  = purrr::map(data, ~ lm(scale(D_phoD) ~ AI, data = .x)),
    tidy   = purrr::map(model, ~ broom::tidy(.x, conf.int = TRUE)),
    glance = purrr::map(model, ~ broom::glance(.x))
  )

coef_df <- model_tbl %>%
  dplyr::select(Metric, tidy) %>%
  tidyr::unnest(tidy) %>%
  dplyr::filter(term == "AI")

stat_df <- model_tbl %>%
  dplyr::select(Metric, glance) %>%
  tidyr::unnest(glance) %>%
  dplyr::select(Metric, r.squared, p.value)

result_df <- coef_df %>%
  dplyr::left_join(stat_df, by = "Metric") %>%
  dplyr::select(Metric, estimate, std.error, conf.low, conf.high, statistic, p.value.x, r.squared) %>%
  dplyr::rename(
    beta    = estimate,
    SE      = std.error,
    CI_low  = conf.low,
    CI_high = conf.high,
    t       = statistic,
    p_value = p.value.x,
    R2      = r.squared
  )

dat_long <- dat_long %>%
  dplyr::group_by(Metric) %>%
  dplyr::mutate(D_phoD_scaled = as.numeric(scale(D_phoD))) %>%
  dplyr::ungroup()

label_df <- result_df

label_pos <- dat_long %>%
  dplyr::group_by(Metric) %>%
  dplyr::summarise(
    x = min(AI, na.rm = TRUE),
    y = max(D_phoD_scaled, na.rm = TRUE),
    .groups = "drop"
  )

label_df <- dplyr::left_join(label_df, label_pos, by = "Metric")

label_df$label <- paste0(
  "italic(R)^2==", sprintf("%.2f", label_df$R2),
  "*','~~italic(P)==", sprintf("%.3f", label_df$p_value)
)

# ---------- 7. theme ----------
theme_s7 <- function() {
  theme_bw(base_size = 10) +
    theme(
      panel.grid      = element_blank(),
      panel.border    = element_rect(colour = "grey25", linewidth = 0.45),
      axis.line       = element_line(colour = "grey25", linewidth = 0.30),
      axis.ticks      = element_line(colour = "grey25", linewidth = 0.30),
      axis.text       = element_text(colour = "black", size = 9),
      axis.title      = element_text(colour = "black", size = 10),
      strip.background = element_rect(fill = "white", colour = "grey25", linewidth = 0.45),
      strip.text      = element_text(face = "bold", size = 9),
      legend.title    = element_text(size = 9),
      legend.text     = element_text(size = 8),
      plot.title      = element_text(face = "bold", size = 10)
    )
}

# ---------- 8. panel A: coefficient plot ----------
pA <- ggplot(result_df, aes(x = Metric, y = beta)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35, colour = "grey50") +
  geom_errorbar(
    aes(ymin = CI_low, ymax = CI_high),
    width = 0.08, linewidth = 0.30, colour = "grey25"
  ) +
  geom_point(
    shape = 21, size = 2.3, stroke = 0.30,
    fill = "white", colour = "grey25"
  ) +
  labs(
    x = NULL,
    y = expression(beta ~ "(AI \u2192 scaled~D"[plot] * ")")
  ) +
  theme_s7()

# ---------- 9. panel B: scatter comparison ----------
pB <- ggplot(dat_long, aes(x = AI, y = D_phoD_scaled)) +
  geom_smooth(
    method = "lm", se = TRUE,
    colour = "grey25", fill = "grey85", linewidth = 0.35
  ) +
  geom_point(
    aes(shape = Duration, fill = Regime),
    size = 2.0, stroke = 0.30, colour = "grey25", alpha = 0.95
  ) +
  facet_wrap(~ Metric, scales = "free_y") +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_fill_manual(values = c(
    "CK"    = "#BFBFBF",
    "NPK"   = "#7FAFB4",
    "NPKOM" = "#D39A92"
  )) +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.1,
    size = 3,
    parse = TRUE
  ) +
  labs(
    x = "Amplification index (AI)",
    y = expression(D[phoD])
  ) +
  theme_s7() +
  theme(legend.position = "right") +
  guides(
    fill  = guide_legend(override.aes = list(shape = 21, size = 3)),
    shape = guide_legend(override.aes = list(fill  = "grey50"))
  )

# ---------- 10. combine ----------
p_s7 <- pA + pB +
  patchwork::plot_layout(widths = c(0.8, 1.8)) +
  patchwork::plot_annotation(tag_levels = "a")

# ---------- 11. save ----------
ggsave(
  filename = file.path(out_dir, "Figure_S6_Dplot_metric_robustness.svg"),
  plot = p_s7,
  width = 8.2, height = 3.6, units = "in"
)

ggsave(
  filename = file.path(out_dir, "Figure_S6_Dplot_metric_robustness.png"),
  plot = p_s7,
  width = 8.2, height = 3.6, units = "in", dpi = 600
)

print(p_s7)
