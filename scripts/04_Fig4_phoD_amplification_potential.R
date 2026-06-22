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
# Fig.4 | Gut microbial state predicts phoD amplification
#         and soil–gut phoD community coupling
#
# 4A: paired soil vs gut logphoD_rel
# 4B: AI ~ Duration × Regime
# 4C: AI ~ gut16S_PC1 + ChemPC1
# 4D: AI ~ D_phoD
# 4E: piecewise structural equation model
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(car)
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
# 2) Read data
# -------------------------
dat_raw <- read_excel(data_xlsx, sheet = "Master_all_data")

need <- c(
  "PlotID", "Duration", "Regime",
  "logphoD_rel_soil", "logphoD_rel_gut",
  "AI", "gut16S_PC1", "ChemPC1", "D_phoD"
)

miss <- setdiff(need, names(dat_raw))
if (length(miss) > 0) {
  stop("Input table missing: ", paste(miss, collapse = ", "))
}

dat <- dat_raw %>%
  transmute(
    PlotID           = trimws(as.character(PlotID)),
    Duration         = as.character(Duration),
    Regime           = as.character(Regime),
    logphoD_rel_soil = as.numeric(logphoD_rel_soil),
    logphoD_rel_gut  = as.numeric(logphoD_rel_gut),
    AI               = as.numeric(AI),
    gut16S_PC1       = as.numeric(gut16S_PC1),
    ChemPC1          = as.numeric(ChemPC1),
    D_phoD           = as.numeric(D_phoD)
  ) %>%
  mutate(
    Duration = trimws(Duration),
    Duration = gsub("\\s+", "", Duration),
    Duration = gsub("years|Years|yr|YR|Year|year", "y", Duration),
    Duration = gsub("^5$", "5y", Duration),
    Duration = gsub("^8$", "8y", Duration),
    Duration = gsub("^10$", "10y", Duration),
    Duration = factor(Duration, levels = c("5y", "8y", "10y")),
    Regime   = factor(Regime, levels = c("CK", "NPK", "NPKOM"))
  ) %>%
  filter(!is.na(Duration), !is.na(Regime))
# -------------------------
# 3) Style
# -------------------------
pal_regime <- c(
  CK    = "#B8B8B8",
  NPK   = "#2CAEB8",
  NPKOM = "#E8897D"
)

shape_duration <- c(
  "5y" = 24,
  "8y" = 22,
  "10y" = 21
)

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}

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
      strip.background = element_rect(fill = "grey90", colour = NA, linewidth = 0),
      strip.text = element_text(size = base_size, colour = "grey10", face = "bold"),
      legend.title = element_text(size = base_size - 0.2, colour = "grey15", face = "bold"),
      legend.text  = element_text(size = base_size - 0.5, colour = "grey25"),
      legend.key   = element_blank(),
      panel.border = element_blank(),
      axis.ticks.length = unit(2, "pt"),
      panel.spacing = unit(1.2, "mm"),
      plot.margin = margin(4, 5, 4, 4, unit = "mm")
    )
}

# ============================================================
# 4A | paired soil vs gut logphoD_rel
# ============================================================

dat_4A <- dat %>%
  select(PlotID, Duration, Regime, logphoD_rel_soil, logphoD_rel_gut) %>%
  filter(is.finite(logphoD_rel_soil), is.finite(logphoD_rel_gut)) %>%
  pivot_longer(
    cols = c(logphoD_rel_soil, logphoD_rel_gut),
    names_to = "Compartment",
    values_to = "logphoD_rel"
  ) %>%
  mutate(
    Compartment = case_when(
      Compartment == "logphoD_rel_soil" ~ "Soil",
      Compartment == "logphoD_rel_gut"  ~ "Gut",
      TRUE ~ Compartment
    ),
    Compartment = factor(Compartment, levels = c("Soil", "Gut"))
  )

lab_4A_df <- dat %>%
  filter(is.finite(logphoD_rel_soil), is.finite(logphoD_rel_gut)) %>%
  group_by(Duration) %>%
  summarise(
    p = t.test(logphoD_rel_gut, logphoD_rel_soil, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("Paired t-test\n", fmt_p_eq(p)),
    x = 0.70,
    y = -2.035
  )

p4A <- ggplot(dat_4A, aes(x = Compartment, y = logphoD_rel)) +
  geom_line(
    aes(group = PlotID),
    colour = "grey70",
    linewidth = 0.28,
    alpha = 0.75
  ) +
  geom_point(
    aes(fill = Regime),
    shape = 21,
    size = 1.8,
    colour = "grey25",
    stroke = 0.25,
    alpha = 0.92
  ) +
  geom_text(
    data = lab_4A_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1.35,
    size = 2.4,
    colour = "grey25",
    lineheight = 1.02
  ) +
  facet_grid(. ~ Duration) +
  scale_fill_manual(values = pal_regime, drop = FALSE) +
  labs(
    x = NULL,
    y = expression(log[10]~italic(phoD)[rel])
  ) +
  coord_cartesian(ylim = c(-2.95, -2.00), clip = "off") +
  theme_nature(base_size = 9) +
  theme(
    legend.position = "none",
    plot.margin = margin(4, 5, 4, 4, unit = "mm")
  )

# ============================================================
# 4B | AI ~ Duration × Regime
# ============================================================

dat_4B <- dat %>%
  filter(is.finite(AI))

fit_4B <- lm(AI ~ Duration * Regime, data = dat_4B)
tab_4B <- car::Anova(fit_4B, type = 2)

lab_4B <- paste0(
  "ANOVA\n",
  "Dur. ", fmt_p_eq(tab_4B["Duration", "Pr(>F)"]),
  "\nReg. ", fmt_p_eq(tab_4B["Regime", "Pr(>F)"]),
  "\nD x R ", fmt_p_eq(tab_4B["Duration:Regime", "Pr(>F)"])
)

y_4B_top <- max(dat_4B$AI, na.rm = TRUE) +
  diff(range(dat_4B$AI, na.rm = TRUE)) * 0.56
y_4B_label <- max(dat_4B$AI, na.rm = TRUE) +
  diff(range(dat_4B$AI, na.rm = TRUE)) * 0.50

p4B <- ggplot(dat_4B, aes(x = Regime, y = AI, fill = Regime)) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.30,
    colour = "grey65"
  ) +
  geom_boxplot(
    width = 0.62,
    outlier.shape = NA,
    linewidth = 0.32,
    colour = "grey25",
    alpha = 0.55
  ) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    size = 1.25,
    alpha = 0.88,
    shape = 21,
    colour = "grey25",
    stroke = 0.22
  ) +
  geom_text(
    data = data.frame(
      Duration = factor("5y", levels = levels(dat_4B$Duration)),
      Regime = factor("CK", levels = levels(dat_4B$Regime)),
      y = y_4B_label,
      label = lab_4B
    ),
    aes(x = Regime, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1.05,
    size = 2.10,
    colour = "grey25",
    lineheight = 1.02
  ) +
  facet_grid(. ~ Duration) +
  scale_fill_manual(values = pal_regime, drop = FALSE) +
  labs(x = NULL, y = "Amplification index (AI)") +
  coord_cartesian(ylim = c(0, y_4B_top), clip = "off") +
  theme_nature(base_size = 9) +
  theme(
    axis.text.x = element_text(
      size = 7.0,
      angle = 25,
      hjust = 1,
      colour = "grey25"
    ),
    legend.position = "none"
  )

# ============================================================
# 4C | AI ~ gut16S_PC1 + ChemPC1
# ============================================================

dat_4C <- dat %>%
  filter(is.finite(AI), is.finite(gut16S_PC1), is.finite(ChemPC1))

fit_4C <- lm(AI ~ gut16S_PC1 + ChemPC1, data = dat_4C)
coef_4C <- summary(fit_4C)$coefficients

t_gut <- coef_4C["gut16S_PC1", "t value"]
p_gut <- coef_4C["gut16S_PC1", "Pr(>|t|)"]
df_res <- df.residual(fit_4C)
partial_r2 <- (t_gut^2) / (t_gut^2 + df_res)

lab_4C <- paste0(
  "Gut 16S effect\n",
  "partial R2 = ", sprintf("%.3f", partial_r2),
  "\n", fmt_p_eq(p_gut)
)

lab_4C <- paste0(
  "AI ~ gut16S_PC1 + ChemPC1\n",
  "partial R² = ", sprintf("%.3f", partial_r2),
  "\nP = ", fmt_p(p_gut)
)

lab_4C <- paste0(
  "Gut 16S effect\n",
  "partial R2 = ", sprintf("%.3f", partial_r2),
  "\n", fmt_p_eq(p_gut)
)

make_partial_effect_df <- function(fit, data, xvar, n = 200, level = 0.95) {
  x_seq <- seq(
    min(data[[xvar]], na.rm = TRUE),
    max(data[[xvar]], na.rm = TRUE),
    length.out = n
  )
  
  V <- vcov(fit)
  crit <- qt(1 - (1 - level) / 2, df = df.residual(fit))
  
  out <- lapply(x_seq, function(x0) {
    newdata <- data
    newdata[[xvar]] <- x0
    
    X <- model.matrix(delete.response(terms(fit)), data = newdata)
    xbar <- colMeans(X)
    
    fit_mean <- as.numeric(xbar %*% coef(fit))
    se_mean <- sqrt(as.numeric(xbar %*% V %*% xbar))
    
    tibble(
      !!xvar := x0,
      fit = fit_mean,
      lwr = fit_mean - crit * se_mean,
      upr = fit_mean + crit * se_mean
    )
  })
  
  bind_rows(out)
}

pred_4C <- make_partial_effect_df(
  fit = fit_4C,
  data = dat_4C,
  xvar = "gut16S_PC1"
)

x_rng_4C <- range(dat_4C$gut16S_PC1, na.rm = TRUE)
y_rng_4C <- range(c(dat_4C$AI, pred_4C$lwr, pred_4C$upr), na.rm = TRUE)
y_4C_top <- y_rng_4C[2] + diff(y_rng_4C) * 0.30
y_4C_label <- y_rng_4C[2] + diff(y_rng_4C) * 0.26

p4C <- ggplot(dat_4C, aes(x = gut16S_PC1, y = AI)) +
  geom_ribbon(
    data = pred_4C,
    aes(x = gut16S_PC1, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,
    fill = "grey55",
    alpha = 0.20
  ) +
  geom_line(
    data = pred_4C,
    aes(x = gut16S_PC1, y = fit),
    inherit.aes = FALSE,
    colour = "grey20",
    linewidth = 0.42
  ) +
  geom_point(
    aes(fill = Regime, shape = Duration),
    size = 2.0,
    colour = "grey25",
    stroke = 0.25,
    alpha = 0.92
  ) +
  annotate(
    "text",
    x = x_rng_4C[1] + diff(x_rng_4C) * 0.01,
    y = y_4C_label,
    label = lab_4C,
    hjust = 0,
    vjust = 1,
    size = 2.20,
    colour = "grey20",
    lineheight = 1.02
  ) +
  scale_fill_manual(values = pal_regime, drop = FALSE) +
  scale_shape_manual(values = shape_duration, drop = FALSE) +
  labs(x = "Gut 16S PC1", y = "Amplification index (AI)") +
  coord_cartesian(ylim = c(y_rng_4C[1], y_4C_top), clip = "off") +
  theme_nature(base_size = 9) +
  theme(legend.position = "none")

# ============================================================
# 4D | AI ~ D_phoD
# ============================================================

dat_4D <- dat %>%
  filter(is.finite(AI), is.finite(D_phoD))

fit_4D <- lm(D_phoD ~ AI, data = dat_4D)
sum_4D <- summary(fit_4D)

r2_4D <- sum_4D$r.squared
p_4D <- sum_4D$coefficients["AI", "Pr(>|t|)"]

lab_4D <- paste0(
  "R² = ", sprintf("%.3f", r2_4D),
  "\nP = ", fmt_p(p_4D)
)

lab_4D <- paste0(
  "R2 = ", sprintf("%.3f", r2_4D),
  "\n", fmt_p_eq(p_4D)
)

x_rng_4D <- range(dat_4D$AI, na.rm = TRUE)
y_rng_4D <- range(dat_4D$D_phoD, na.rm = TRUE)
y_4D_top <- y_rng_4D[2] + diff(y_rng_4D) * 0.30
y_4D_label <- y_rng_4D[2] + diff(y_rng_4D) * 0.26

p4D <- ggplot(dat_4D, aes(x = AI, y = D_phoD)) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.42,
    colour = "grey20",
    fill = "grey55",
    alpha = 0.18
  ) +
  geom_point(
    aes(fill = Regime, shape = Duration),
    size = 2.0,
    colour = "grey25",
    stroke = 0.25,
    alpha = 0.92
  ) +
  annotate(
    "text",
    x = x_rng_4D[1] + diff(x_rng_4D) * 0.01,
    y = y_4D_label,
    label = lab_4D,
    hjust = 0,
    vjust = 1,
    size = 2.20,
    colour = "grey20",
    lineheight = 1.02
  ) +
  scale_fill_manual(values = pal_regime, drop = FALSE) +
  scale_shape_manual(values = shape_duration, drop = FALSE) +
  guides(
    fill = "none",
    shape = guide_legend(override.aes = list(fill = "white", colour = "grey25", size = 2.4))
  ) +
  labs(
    x = "Amplification index (AI)",
    y = expression("Soil-gut "*italic(phoD)*" dissimilarity")
  ) +
  coord_cartesian(ylim = c(y_rng_4D[1], y_4D_top), clip = "off") +
  theme_nature(base_size = 9) +
  theme(
    legend.position = c(0.78, 0.86),
    legend.justification = c(0, 0.5),
    legend.background = element_blank(),
    legend.key.size = unit(3.2, "mm"),
    legend.title = element_text(size = 8.2, colour = "grey15", face = "bold"),
    legend.text = element_text(size = 7.2, colour = "grey25")
  )

# ============================================================
# 4E | Piecewise structural equation model
# ============================================================

# ============================================================
# 4E | Piecewise structural equation model
# ============================================================

dat_sem <- dat %>%
  filter(
    is.finite(ChemPC1),
    is.finite(gut16S_PC1),
    is.finite(AI),
    is.finite(D_phoD)
  )

m1 <- lm(gut16S_PC1 ~ ChemPC1, data = dat_sem)
m2 <- lm(AI ~ gut16S_PC1 + ChemPC1, data = dat_sem)
m3 <- lm(D_phoD ~ AI, data = dat_sem)

if (requireNamespace("piecewiseSEM", quietly = TRUE)) {
  sem_fit <- piecewiseSEM::psem(m1, m2, m3)
  fc <- piecewiseSEM::fisherC(sem_fit)
} else {
  basis_p <- c(
    summary(lm(D_phoD ~ AI + ChemPC1, data = dat_sem))$coefficients["ChemPC1", "Pr(>|t|)"],
    summary(lm(D_phoD ~ AI + gut16S_PC1, data = dat_sem))$coefficients["gut16S_PC1", "Pr(>|t|)"]
  )
  basis_p <- pmax(pmin(basis_p, 1), .Machine$double.xmin)
  fisher_c <- -2 * sum(log(basis_p))
  fisher_df <- 2 * length(basis_p)

  fc <- tibble(
    Fisher.C = fisher_c,
    df = fisher_df,
    P.Value = pchisq(fisher_c, df = fisher_df, lower.tail = FALSE)
  )
}

fit_lab <- paste0(
  "Fisher's C = ", sprintf("%.2f", fc$Fisher.C[1]),
  ", df = ", fc$df[1],
  ", ", fmt_p_eq(fc$P.Value[1])
)

get_beta_p <- function(fit, term) {
  sm <- summary(fit)$coefficients
  mf <- model.frame(fit)
  y <- mf[[1]]
  x <- mf[[term]]
  
  b <- sm[term, "Estimate"]
  p <- sm[term, "Pr(>|t|)"]
  beta <- as.numeric(b * sd(x, na.rm = TRUE) / sd(y, na.rm = TRUE))
  
  list(beta = beta, p = p)
}

bp_CG <- get_beta_p(m1, "ChemPC1")
bp_GA <- get_beta_p(m2, "gut16S_PC1")
bp_CA <- get_beta_p(m2, "ChemPC1")
bp_AD <- get_beta_p(m3, "AI")

edge_lab <- function(x) {
  paste0("β=", sprintf("%.3f", x$beta), "\nP=", fmt_p(x$p))
}

edge_lab <- function(x) {
  paste0("beta = ", sprintf("%.3f", x$beta), "\n", fmt_p_eq(x$p))
}

nodes <- tibble(
  name = c("Soil chemistry\nChemPC1", "Gut 16S\nPC1", "Amplification\nindex", "phoD\ncoupling"),
  x = c(0, 1.15, 2.30, 3.40),
  y = c(0, 0, 0, 0)
)

edges <- tibble(
  from = c("ChemPC1", "gut16S_PC1", "AI", "ChemPC1"),
  to   = c("gut16S_PC1", "AI", "D_phoD", "AI"),
  x    = c(0, 1.15, 2.30, 0),
  y    = c(0, 0, 0, 0),
  xend = c(1.15, 2.30, 3.40, 2.30),
  yend = c(0, 0, 0, 0),
  lab  = c(edge_lab(bp_CG), edge_lab(bp_GA), edge_lab(bp_AD), edge_lab(bp_CA)),
  curved = c(FALSE, FALSE, FALSE, TRUE),
  beta = c(bp_CG$beta, bp_GA$beta, bp_AD$beta, bp_CA$beta),
  p = c(bp_CG$p, bp_GA$p, bp_AD$p, bp_CA$p)
) %>%
  mutate(
    col = case_when(
      p < 0.05 & beta > 0 ~ "#4F6D8A",
      p < 0.05 & beta < 0 ~ "#C06C6C",
      TRUE ~ "grey70"
    ),
    lw = 0.35 + pmin(abs(beta), 1) * 0.45,
    alpha = ifelse(p < 0.05, 1, 0.55)
  )

p4E <- ggplot() +
  geom_segment(
    data = edges %>% filter(!curved),
    aes(x = x, y = y, xend = xend, yend = yend,
        colour = col, linewidth = lw, alpha = alpha),
    lineend = "round",
    arrow = arrow(length = unit(2.1, "mm"), type = "closed"),
    show.legend = FALSE
  ) +
  geom_curve(
    data = edges %>% filter(curved),
    aes(x = x, y = y - 0.05, xend = xend, yend = yend - 0.05,
        colour = col, linewidth = lw, alpha = alpha),
    curvature = +0.36,
    lineend = "round",
    arrow = arrow(length = unit(2.1, "mm"), type = "closed"),
    show.legend = FALSE
  ) +
  geom_label(
    data = nodes,
    aes(x = x, y = y, label = name),
    size = 3.0,
    fill = "white",
    colour = "grey20",
    label.size = 0.25,
    label.r = unit(0, "lines"),
    label.padding = unit(0.16, "lines")
  ) +
  geom_text(
    data = edges %>% filter(!curved),
    aes(x = (x + xend) / 2, y = 0.24, label = lab),
    size = 2.35,
    colour = "grey20",
    lineheight = 0.95
  ) +
  geom_text(
    data = edges %>% filter(curved),
    aes(x = (x + xend) / 2, y = -0.36, label = lab),
    size = 2.35,
    colour = "grey20",
    lineheight = 0.95
  ) +
  annotate(
    "text",
    x = -0.20,
    y = 0.60,
    label = "Piecewise SEM",
    hjust = 0,
    size = 3.0,
    colour = "grey20",
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 3.55,
    y = 0.60,
    label = fit_lab,
    hjust = 1,
    size = 2.5,
    colour = "grey35"
  ) +
  scale_colour_identity() +
  scale_linewidth_identity() +
  scale_alpha_identity() +
  coord_cartesian(
    xlim = c(-0.25, 3.82),
    ylim = c(-0.55, 0.62),
    clip = "off"
  )+
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(2, 5, 2, 5, unit = "mm")
  )

# ============================================================
# Combine: align panels by columns + bottom SEM
# ============================================================

left_col <- cowplot::plot_grid(
  p4A,
  p4C,
  ncol = 1,
  align = "v",
  axis = "lr",
  labels = c("a", "c"),
  label_size = 10,
  label_fontface = "bold",
  label_colour = "grey15"
)

right_col <- cowplot::plot_grid(
  p4B,
  p4D,
  ncol = 1,
  align = "v",
  axis = "lr",
  labels = c("b", "d"),
  label_size = 10,
  label_fontface = "bold",
  label_colour = "grey15"
)

main_grid <- cowplot::plot_grid(
  left_col,
  right_col,
  ncol = 2,
  align = "h",
  axis = "tb",
  rel_widths = c(1, 1)
)

sem_row <- cowplot::plot_grid(
  p4E,
  ncol = 1,
  labels = "e",
  label_size = 10,
  label_fontface = "bold",
  label_colour = "grey15"
)

fig4 <- cowplot::plot_grid(
  main_grid,
  sem_row,
  ncol = 1,
  rel_heights = c(2, 0.52)
)

print(fig4)

# ============================================================
# Save
# ============================================================

ggsave(
  filename = file.path(out_dir, "Fig4_gut_state_phoD_amplification_coupling_DphoD_SEM_polished.svg"),
  plot = fig4,
  width = 185,
  height = 175,
  units = "mm",
  device = svg_device,
  bg = "white"
)

ggsave(
  filename = file.path(out_dir, "Fig4_gut_state_phoD_amplification_coupling_DphoD_SEM_polished.png"),
  plot = fig4,
  width = 185,
  height = 175,
  units = "mm",
  dpi = 600,
  bg = "white"
)

cat("Saved Fig.4 to: ", out_dir, "\n")
