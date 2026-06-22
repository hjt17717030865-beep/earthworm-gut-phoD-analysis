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
# Fig.2 | Soil-state background
# Layout:
#   left: ChemPC1
#   right top: Olsen-P
#   right bottom: pH
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(car)
  library(cowplot)
  library(grid)
})

svg_device <- if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite
} else {
  grDevices::svg
}

# 1) paths ----
xlsx_path <- DATA_XLSX
out_dir <- OUTPUT_DIR
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
# 2) read data ----
dat <- read_excel(xlsx_path, sheet = "Master_all_data") %>%
  mutate(
    Duration = as.character(Duration),
    Duration = ifelse(Duration == "3y", "5y", Duration),
    Duration = factor(Duration, levels = c("5y", "8y", "10y")),
    Regime   = factor(Regime, levels = c("CK", "NPK", "NPKOM")),
    ChemPC1  = as.numeric(ChemPC1),
    `Olsen-P(mg/kg)` = as.numeric(`Olsen-P(mg/kg)`),
    pH = as.numeric(pH)
  ) %>%
  filter(
    is.finite(ChemPC1),
    is.finite(`Olsen-P(mg/kg)`),
    is.finite(pH),
    !is.na(Duration),
    !is.na(Regime)
  ) %>%
  droplevels()
# 3) functions ----
get_p <- function(model, term) {
  aov_tab <- as.data.frame(car::Anova(model, type = 2))
  aov_tab$term <- rownames(aov_tab)
  p <- aov_tab$`Pr(>F)`[aov_tab$term == term]
  if (length(p) == 0) return(NA_real_)
  p[1]
}

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}

fmt_p_eq <- function(p) {
  if (is.na(p)) return("P = NA")
  if (p < 0.001) return("P < 0.001")
  paste0("P = ", sprintf("%.3f", p))
}

make_label <- function(model) {
  paste0(
    "Type-II ANOVA P\n",
    "Dur=", fmt_p(get_p(model, "Duration")),
    "  Reg=", fmt_p(get_p(model, "Regime")),
    "\nDur×Reg=", fmt_p(get_p(model, "Duration:Regime"))
  )
}

make_label_nature <- function(model) {
  paste0(
    "ANOVA\n",
    "Dur. ", fmt_p_eq(get_p(model, "Duration")),
    "\nReg. ", fmt_p_eq(get_p(model, "Regime")),
    "\nD x R ", fmt_p_eq(get_p(model, "Duration:Regime"))
  )
}

# 4) models ----
m_chem  <- lm(ChemPC1 ~ Duration * Regime, data = dat)
m_olsen <- lm(`Olsen-P(mg/kg)` ~ Duration * Regime, data = dat)
m_ph    <- lm(pH ~ Duration * Regime, data = dat)

label_chem  <- make_label_nature(m_chem)
label_olsen <- make_label_nature(m_olsen)
label_ph    <- make_label_nature(m_ph)

print(car::Anova(m_chem, type = 2))
print(car::Anova(m_olsen, type = 2))
print(car::Anova(m_ph, type = 2))

# 5) colours ----
pal_fill <- c(
  CK    = "#B8B8B8",
  NPK   = "#2CAEB8",
  NPKOM = "#E8897D"
)

# 6) theme ----
theme_nature_fig2 <- function(base_size = 9.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title   = element_text(colour = "grey18", size = base_size),
      axis.text    = element_text(colour = "grey25", size = base_size - 1),
      axis.text.x  = element_text(
        colour = "grey25",
        size = base_size - 2.2,
        angle = 35,
        hjust = 1
      ),
      axis.line    = element_line(linewidth = 0.35, colour = "grey25"),
      axis.ticks   = element_line(linewidth = 0.35, colour = "grey25"),
      legend.title = element_blank(),
      legend.text  = element_text(colour = "grey25", size = base_size - 0.5),
      legend.key   = element_blank(),
      plot.title   = element_text(colour = "grey15", size = base_size + 0.5, face = "bold"),
      strip.background = element_rect(fill = "grey90", colour = NA, linewidth = 0),
      strip.text = element_text(colour = "grey10", size = base_size, face = "bold"),
      panel.spacing = unit(1.4, "mm"),
      plot.margin  = margin(5, 5, 5, 5, unit = "mm")
    )
}

pos_d <- position_dodge(width = 0.70)

make_boxplot_panel <- function(
  data, yvar, ylab, label_txt, base_size = 9.5,
  top_pad = 0.44, label_pad = 0.37
) {
  plot_data <- data

  y <- data[[yvar]]
  y_rng <- range(y, na.rm = TRUE)
  y_span <- diff(y_rng)
  y_top <- y_rng[2] + y_span * top_pad
  y_label <- y_rng[2] + y_span * label_pad

  label_data <- data.frame(
    Duration = factor("5y", levels = levels(plot_data$Duration)),
    Regime = factor("CK", levels = levels(plot_data$Regime)),
    y = y_label,
    label = label_txt
  )
  
  ggplot(plot_data, aes(x = Regime, y = .data[[yvar]], fill = Regime)) +
    geom_boxplot(
      alpha = 0.42,
      linewidth = 0.35,
      colour = "grey25",
      outlier.shape = NA,
      width = 0.62
    ) +
    geom_point(
      position = position_jitter(width = 0.08, height = 0),
      size = 1.25,
      alpha = 0.85,
      shape = 21,
      colour = "grey25",
      stroke = 0.20
    ) +
    geom_text(
      data = label_data,
      aes(x = Regime, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.12,
      colour = "grey25",
      lineheight = 1.02
    ) +
    facet_grid(. ~ Duration) +
    scale_fill_manual(values = pal_fill) +
    labs(x = NULL, y = ylab) +
    coord_cartesian(ylim = c(y_rng[1], y_top), clip = "off") +
    theme_nature_fig2(base_size = base_size)
}

# 7) panels ----
p_pc1 <- make_boxplot_panel(
  dat,
  yvar = "ChemPC1",
  ylab = "ChemPC1 (PCA score)",
  label_txt = label_chem,
  base_size = 10,
  top_pad = 0.34,
  label_pad = 0.28
) +
  theme(legend.position = "none")

p_olsen <- make_boxplot_panel(
  dat,
  yvar = "Olsen-P(mg/kg)",
  ylab = "Olsen-P (mg kg-1)",
  label_txt = label_olsen,
  base_size = 8.5,
  top_pad = 0.68,
  label_pad = 0.60
) +
  theme(
    legend.position = "none"
  )

p_ph <- make_boxplot_panel(
  dat,
  yvar = "pH",
  ylab = "pH",
  label_txt = label_ph,
  base_size = 8.5,
  top_pad = 0.68,
  label_pad = 0.60
) +
  theme(
    legend.position = "none"
  )

# 9) combine ----
right_panels <- plot_grid(
  p_olsen,
  p_ph,
  ncol = 1,
  rel_heights = c(1, 1),
  labels = c("b", "c"),
  label_size = 10,
  label_fontface = "bold",
  label_colour = "grey15"
)

main_panels <- plot_grid(
  p_pc1,
  right_panels,
  ncol = 2,
  rel_widths = c(1.24, 1),
  labels = c("a", ""),
  label_size = 10,
  label_fontface = "bold",
  label_colour = "grey15"
)

fig2 <- main_panels

print(fig2)

# 10) export ----
ggsave(
  filename = file.path(out_dir, "Fig2_soil_state_PC1_OlsenP_pH_polished.png"),
  plot = fig2,
  width = 7.6,
  height = 4.6,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(out_dir, "Fig2_soil_state_PC1_OlsenP_pH_polished.svg"),
  plot = fig2,
  width = 7.6,
  height = 4.6,
  device = svg_device,
  bg = "white"
)
