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
# Fig.5 | Gut filtering reconfigures broader P-acquisition potential
#
# 5A: paired soil vs gut P-index
# 5B: module-level paired ΔP estimates
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(grid)
})

svg_device <- if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite
} else {
  grDevices::svg
}

# -------------------------
# -------------------------
# 1) Paths
# -------------------------
data_xlsx <- DATA_XLSX
out_dir <- OUTPUT_DIR
if (!file.exists(data_xlsx)) stop("File not found: ", data_xlsx)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
# -------------------------
# 2) Style
# -------------------------
pal_regime <- c(
  CK    = "#B8B8B8",
  NPK   = "#2CAEB8",
  NPKOM = "#E8897D"
)

fmt_p_eq <- function(p) {
  ifelse(
    is.na(p),
    "P = NA",
    ifelse(p < 0.001, "P < 0.001", paste0("P = ", sprintf("%.3f", p)))
  )
}

theme_nature <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size, colour = "grey20"),
      axis.text  = element_text(size = base_size - 0.5, colour = "grey25"),
      axis.line  = element_line(linewidth = 0.30, colour = "grey25"),
      axis.ticks = element_line(linewidth = 0.30, colour = "grey25"),
      legend.title = element_blank(),
      legend.text  = element_text(size = base_size - 0.5, colour = "grey25"),
      legend.key   = element_blank(),
      strip.background = element_rect(fill = "grey90", colour = NA, linewidth = 0),
      strip.text = element_text(size = base_size, colour = "grey10", face = "bold"),
      panel.border = element_blank(),
      axis.ticks.length = unit(2, "pt"),
      panel.spacing = unit(1.2, "mm"),
      plot.margin = margin(4, 5, 4, 4, unit = "mm")
    )
}

# -------------------------
# 3) Read data
# -------------------------
p_sample <- read_excel(data_xlsx, sheet = "Master_all_data")
p_gene   <- read_excel(data_xlsx, sheet = "P_gene_analysis_long")

# ============================================================
# 5A | paired soil vs gut P-index
# ============================================================
p_sample2 <- p_sample %>%
  mutate(
    PlotID       = as.character(PlotID),
    Duration     = as.character(Duration),
    Duration     = ifelse(Duration == "3y", "5y", Duration),
    Duration     = factor(Duration, levels = c("5y", "8y", "10y")),
    Regime       = factor(as.character(Regime), levels = c("CK", "NPK", "NPKOM")),
    P_index_soil = as.numeric(P_index_soil),
    P_index_gut  = as.numeric(P_index_gut)
  ) %>%
  filter(is.finite(P_index_soil), is.finite(P_index_gut))

p_wide <- p_sample2 %>%
  select(PlotID, Duration, Regime, P_index_soil, P_index_gut) %>%
  distinct() %>%
  mutate(Delta_P_index = P_index_gut - P_index_soil)

p_long <- p_wide %>%
  pivot_longer(
    cols      = c(P_index_soil, P_index_gut),
    names_to  = "Habitat",
    values_to = "P_index"
  ) %>%
  mutate(
    Habitat = case_when(
      Habitat == "P_index_soil" ~ "Soil",
      Habitat == "P_index_gut"  ~ "Gut"
    ),
    Habitat = factor(Habitat, levels = c("Soil", "Gut"))
  )
lab_5A <- p_wide %>%
  group_by(Duration) %>%
  summarise(
    p = t.test(P_index_gut, P_index_soil, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = fmt_p_eq(p),
    x = 1.50
  )

mean_5A <- p_long %>%
  group_by(Duration, Habitat) %>%
  summarise(
    mean_p = mean(P_index, na.rm = TRUE),
    .groups = "drop"
  )

y_rng_5A <- range(p_long$P_index, na.rm = TRUE)
y_5A_top <- y_rng_5A[2] + diff(y_rng_5A) * 0.19
y_5A_bracket <- y_rng_5A[2] + diff(y_rng_5A) * 0.11
y_5A_tick <- diff(y_rng_5A) * 0.030
y_5A_label <- y_rng_5A[2] + diff(y_rng_5A) * 0.145

lab_5A <- lab_5A %>%
  mutate(
    y = y_5A_label,
    y_bracket = y_5A_bracket,
    y_bracket_low = y_5A_bracket - y_5A_tick
  )

p5A <- ggplot(p_long, aes(x = Habitat, y = P_index)) +
  geom_line(
    aes(group = PlotID),
    colour = "grey78",
    linewidth = 0.24,
    alpha = 0.68
  ) +
  geom_line(
    data = mean_5A,
    aes(y = mean_p, group = Duration),
    colour = "grey12",
    linewidth = 0.68,
    inherit.aes = TRUE
  ) +
  geom_point(
    aes(fill = Regime),
    shape = 21,
    size = 1.75,
    colour = "grey25",
    stroke = 0.22,
    alpha = 0.90
  ) +
  geom_point(
    data = mean_5A,
    aes(y = mean_p),
    shape = 23,
    size = 2.35,
    fill = "white",
    colour = "grey12",
    stroke = 0.42,
    inherit.aes = TRUE
  ) +
  geom_segment(
    data = lab_5A,
    aes(x = 1, xend = 2, y = y_bracket, yend = y_bracket),
    linewidth = 0.30,
    colour = "grey25",
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = lab_5A,
    aes(x = 1, xend = 1, y = y_bracket_low, yend = y_bracket),
    linewidth = 0.30,
    colour = "grey25",
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = lab_5A,
    aes(x = 2, xend = 2, y = y_bracket_low, yend = y_bracket),
    linewidth = 0.30,
    colour = "grey25",
    inherit.aes = FALSE
  ) +
  geom_text(
    data = lab_5A,
    aes(x = x, y = y, label = label),
    hjust = 0.5,
    vjust = 0,
    size = 2.25,
    colour = "grey25",
    inherit.aes = FALSE
  ) +
  facet_grid(. ~ Duration) +
  scale_fill_manual(values = pal_regime, drop = FALSE) +
  scale_x_discrete(expand = expansion(mult = c(0.20, 0.20))) +
  labs(x = NULL, y = "P-index") +
  coord_cartesian(ylim = c(y_rng_5A[1], y_5A_top), clip = "off") +
  theme_nature(base_size = 9) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8.6, colour = "grey10", face = "bold"),
    axis.text.x = element_text(size = 8.4, colour = "grey20"),
    plot.margin = margin(4, 5, 4, 4, unit = "mm")
  )

# ============================================================
# 5B | module-level paired ΔP estimates
# ============================================================

module_map <- tibble(
  gene = c("gcd", "pqqC", "phoX", "phnK", "ppk", "ppx"),
  module = c(
    "Inorganic P\nsolubilization",
    "Inorganic P\nsolubilization",
    "Organic P\nmineralization",
    "Phosphonate\nuse",
    "PolyP\nturnover",
    "PolyP\nturnover"
  )
)

gene_scaled <- p_gene %>%
  mutate(
    PlotID = as.character(PlotID),
    Duration = as.character(Duration),
    Duration = ifelse(Duration == "3y", "5y", Duration),
    Duration = factor(Duration, levels = c("5y", "8y", "10y")),
    Regime = factor(as.character(Regime), levels = c("CK", "NPK", "NPKOM")),
    Habitat = tolower(as.character(Habitat)),
    gene = as.character(gene),
    copies = as.numeric(copies)
  ) %>%
  filter(
    Habitat %in% c("soil", "gut"),
    gene %in% module_map$gene,
    is.finite(copies)
  ) %>%
  left_join(module_map, by = "gene") %>%
  group_by(gene) %>%
  mutate(
    min_pos = suppressWarnings(min(copies[copies > 0], na.rm = TRUE)),
    pseudo = ifelse(is.finite(min_pos), min_pos / 2, 1e-12),
    log_copies = log10(copies + pseudo),
    gene_z = as.numeric(scale(log_copies))
  ) %>%
  ungroup()

module_sample <- gene_scaled %>%
  group_by(PlotID, Duration, Regime, Habitat, module) %>%
  summarise(
    module_index = mean(gene_z, na.rm = TRUE),
    .groups = "drop"
  )

module_delta <- module_sample %>%
  pivot_wider(
    names_from = Habitat,
    values_from = module_index,
    names_prefix = "module_"
  ) %>%
  filter(is.finite(module_soil), is.finite(module_gut)) %>%
  mutate(
    Delta_module = module_gut - module_soil
  )
# -------------------------
# Module effect test
# -------------------------
fit_mod <- lm(Delta_module ~ module, data = module_delta)
anova_res <- anova(fit_mod)

p_module <- anova_res["module", "Pr(>F)"]

lab_module_text <- paste0(
  "Module effect\n", fmt_p_eq(p_module)
)
module_sum <- module_delta %>%
  group_by(module) %>%
  summarise(
    n = n(),
    mean_delta = mean(Delta_module, na.rm = TRUE),
    se = sd(Delta_module, na.rm = TRUE) / sqrt(n),
    ci_low = mean_delta - qt(0.975, df = n - 1) * se,
    ci_high = mean_delta + qt(0.975, df = n - 1) * se,
    p = t.test(Delta_module, mu = 0)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    module = factor(
      module,
      levels = module[order(mean_delta)]
    ),
    p_label = fmt_p_eq(p)
  )

x_rng_5B <- range(c(module_sum$ci_low, module_sum$ci_high), na.rm = TRUE)
x_pad_5B <- diff(x_rng_5B) * 0.12
x_5B_min <- min(0, x_rng_5B[1]) - x_pad_5B * 0.55
x_5B_max <- x_rng_5B[2] + x_pad_5B * 1.65
x_5B_lab_pad <- diff(x_rng_5B) * 0.045

p5B <- ggplot(module_sum, aes(x = mean_delta, y = module)) +
  geom_hline(
    yintercept = seq_along(levels(module_sum$module)),
    colour = "grey93",
    linewidth = 0.26
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.32,
    linetype = 2,
    colour = "grey55"
  ) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    width = 0.18,
    orientation = "y",
    linewidth = 0.45,
    colour = "grey25"
  ) +
  geom_point(
    size = 2.3,
    shape = 21,
    fill = "grey80",
    colour = "grey25",
    stroke = 0.35
  ) +
  geom_text(
    aes(x = ci_high + x_5B_lab_pad, label = p_label),
    hjust = 0,
    size = 2.35,
    colour = "grey25"
  ) +
  labs(
    x = expression(Delta~"module index (gut - soil)"),
    y = NULL
  ) +
  annotate(
    "text",
    x = 0.12,
    y = length(unique(module_sum$module)) + 0.55,
    label = lab_module_text,
    hjust = 0,
    size = 2.35,
    colour = "grey30",
    lineheight = 1.02
  ) +
  coord_cartesian(xlim = c(x_5B_min, x_5B_max), clip = "off") +
  theme_nature(base_size = 9) +
  theme(
    legend.position = "none",
    plot.margin = margin(5, 10, 4, 4, unit = "mm")
  )

# ============================================================
# Combine
# ============================================================

fig5 <- cowplot::plot_grid(
  p5A,
  p5B,
  ncol = 2,
  rel_widths = c(1.12, 1.00),
  align = "hv",
  axis = "tb",
  labels = c("a", "b"),
  label_size = 10,
  label_fontface = "bold",
  label_colour = "grey15"
)

print(fig5)

# ============================================================
# Save
# ============================================================

ggsave(
  filename = file.path(out_dir, "Fig5_P_acquisition_context_polished.svg"),
  plot = fig5,
  width = 183,
  height = 82,
  units = "mm",
  device = svg_device,
  bg = "white"
)

ggsave(
  filename = file.path(out_dir, "Fig5_P_acquisition_context_polished.png"),
  plot = fig5,
  width = 183,
  height = 82,
  units = "mm",
  dpi = 600,
  bg = "white"
)

write.csv(
  p_wide,
  file = file.path(out_dir, "Fig5_paired_Pindex_table.csv"),
  row.names = FALSE
)

write.csv(
  module_delta,
  file = file.path(out_dir, "Fig5_module_delta_table.csv"),
  row.names = FALSE
)

write.csv(
  module_sum,
  file = file.path(out_dir, "Fig5_module_delta_summary.csv"),
  row.names = FALSE
)

cat("Saved Fig.5 to: ", out_dir, "\n")
