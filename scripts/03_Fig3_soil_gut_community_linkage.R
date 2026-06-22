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
# Fig.3 | Gut filtering
# 3A: soil vs gut 16S ordination (paired Aitchison PCA)
# 3B: paired soil–gut 16S Aitchison distance
# 3C: gut16S_PC1 ~ Duration × Regime
#
# Revised version:
# - keeps main Fig.3 structure
# - absorbs the cleaner S5 style
# - removes the need for a separate S5 figure
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(vegan)
  library(cowplot)
  library(readxl)
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
out_dir <- OUTPUT_DIR
supp_xlsx <- DATA_XLSX
if (!file.exists(supp_xlsx)) stop("Supplementary Data file not found: ", supp_xlsx)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
# -------------------------
# 2) Style
# -------------------------
pal_regime <- c(
  CK    = "#B8B8B8",
  NPK   = "#2CAEB8",
  NPKOM = "#E8897D"
)

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

theme_nature <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size, colour = "grey20"),
      axis.text  = element_text(size = base_size - 0.5, colour = "grey25"),
      axis.line  = element_line(linewidth = 0.30, colour = "grey25"),
      axis.ticks = element_line(linewidth = 0.30, colour = "grey25"),
      axis.ticks.length = unit(2, "pt"),
      legend.title = element_text(size = base_size - 0.2, colour = "grey15", face = "bold"),
      legend.text  = element_text(size = base_size - 0.5, colour = "grey25"),
      legend.key   = element_blank(),
      legend.spacing.y = unit(1, "pt"),
      strip.background = element_rect(fill = "grey90", colour = NA, linewidth = 0),
      strip.text = element_text(size = base_size, colour = "grey10", face = "bold"),
      panel.border = element_blank(),
      panel.spacing = unit(1.2, "mm"),
      plot.margin = margin(4, 5, 4, 4, unit = "mm")
    )
}

# -------------------------
# 3) Read 16S ASV tables
# -------------------------
soil_raw <- readxl::read_excel(supp_xlsx, sheet = "16S_soil_ASV_table") %>%
  as.data.frame()

gut_raw <- readxl::read_excel(supp_xlsx, sheet = "16S_gut_ASV_table") %>%
  as.data.frame()

meta_cols <- c("PlotID", "Duration", "Regime", "Compartment")

soil_asv_cols <- setdiff(colnames(soil_raw), meta_cols)
gut_asv_cols  <- setdiff(colnames(gut_raw),  meta_cols)

soil <- as.matrix(sapply(soil_raw[, soil_asv_cols, drop = FALSE], as.numeric))
gut  <- as.matrix(sapply(gut_raw[,  gut_asv_cols,  drop = FALSE], as.numeric))

rownames(soil) <- trimws(as.character(soil_raw$PlotID))
rownames(gut)  <- trimws(as.character(gut_raw$PlotID))

common_pairs <- intersect(rownames(soil), rownames(gut))
common_pairs <- sort(common_pairs)

if (length(common_pairs) == 0) {
  stop("No paired sample names shared between soil and gut.")
}

soil <- soil[common_pairs, , drop = FALSE]
gut  <- gut[common_pairs,  , drop = FALSE]

# -------------------------
# 4) Same ASV space
# -------------------------
all_asv_names <- union(colnames(soil), colnames(gut))

pad_matrix <- function(mat, all_cols) {
  missing_cols <- setdiff(all_cols, colnames(mat))
  if (length(missing_cols) > 0) {
    zero_block <- matrix(
      0,
      nrow = nrow(mat),
      ncol = length(missing_cols),
      dimnames = list(rownames(mat), missing_cols)
    )
    mat <- cbind(mat, zero_block)
  }
  mat[, all_cols, drop = FALSE]
}

soil_pad <- pad_matrix(soil, all_asv_names)
gut_pad  <- pad_matrix(gut,  all_asv_names)

# -------------------------
# 5) CLR transform
# -------------------------
clr_transform <- function(mat) {
  mat <- as.matrix(mat) + 1
  log_mat <- log(mat)
  sweep(log_mat, 1, rowMeans(log_mat), "-")
}

soil_clr <- clr_transform(soil_pad)
gut_clr  <- clr_transform(gut_pad)

rownames(soil_clr) <- paste0("Soil__", rownames(soil_clr))
rownames(gut_clr)  <- paste0("Gut__",  rownames(gut_clr))

clr_mat <- rbind(soil_clr, gut_clr)

# -------------------------
# 6) Metadata
# -------------------------
meta <- tibble(
  SampleID = rownames(clr_mat),
  Type     = ifelse(grepl("^Soil__", SampleID), "Soil", "Gut"),
  PairID   = sub("^(Soil|Gut)__", "", SampleID)
) %>%
  mutate(
    Type = factor(Type, levels = c("Soil", "Gut")),
    Duration = sub("^([0-9]+y).*", "\\1", PairID),
    Duration = ifelse(Duration == "3y", "5y", Duration),
    Duration = factor(Duration, levels = c("5y", "8y", "10y")),
    Regime = sub("^[0-9]+y(.*)-[0-9]+$", "\\1", PairID),
    Regime = factor(Regime, levels = c("CK", "NPK", "NPKOM"))
  )

# -------------------------
# 7) Aitchison PCA
# -------------------------
pca_all <- prcomp(clr_mat, center = TRUE, scale. = FALSE)
var_exp <- 100 * (pca_all$sdev^2 / sum(pca_all$sdev^2))

ord_df <- meta %>%
  mutate(
    PC1 = pca_all$x[, 1],
    PC2 = pca_all$x[, 2]
  ) %>%
  arrange(PairID, Type)

# -------------------------
# 8) PERMANOVA + dispersion
# -------------------------
dist_mat <- dist(clr_mat, method = "euclidean")

set.seed(123)
perm <- adonis2(
  dist_mat ~ Type,
  data = as.data.frame(meta),
  permutations = 999,
  strata = meta$PairID
)

perm_F  <- perm$F[1]
perm_R2 <- perm$R2[1]
perm_P  <- perm$`Pr(>F)`[1]

bd <- betadisper(dist_mat, group = meta$Type)
set.seed(123)
bd_perm <- permutest(bd, permutations = 999)

disp_F <- bd_perm$tab[1, "F"]
disp_P <- bd_perm$tab[1, "Pr(>F)"]

stat_lab <- paste0(
  "PERMANOVA: R² = ", sprintf("%.3f", perm_R2),
  ", F = ", sprintf("%.2f", perm_F),
  ", P = ", fmt_p(perm_P), "\n",
  "Dispersion test: F = ", sprintf("%.2f", disp_F),
  ", P = ", fmt_p(disp_P)
)

# -------------------------
# -------------------------
# -------------------------
pair_dist_df <- readxl::read_excel(supp_xlsx, sheet = "Master_all_data") %>%
  transmute(
    PairID             = trimws(as.character(PlotID)),
    Duration           = factor(as.character(Duration), levels = c("5y","8y","10y")),
    Regime             = factor(as.character(Regime),   levels = c("CK","NPK","NPKOM")),
    Aitchison_distance = as.numeric(D_16S_Aitchison)
  )
# -------------------------
# 10) Fig.3A paired ordination
# -------------------------
stat_lab <- paste0(
  "PERMANOVA: R2 = ", sprintf("%.3f", perm_R2),
  ", F = ", sprintf("%.2f", perm_F),
  ", ", fmt_p_eq(perm_P), "\n",
  "Dispersion test: F = ", sprintf("%.2f", disp_F),
  ", ", fmt_p_eq(disp_P)
)

segments_df <- ord_df %>%
  dplyr::select(PairID, Type, PC1, PC2)%>%
  pivot_wider(names_from = Type, values_from = c(PC1, PC2)) %>%
  filter(
    is.finite(PC1_Soil), is.finite(PC2_Soil),
    is.finite(PC1_Gut),  is.finite(PC2_Gut)
  )

x_rng <- range(ord_df$PC1, na.rm = TRUE)
y_rng <- range(ord_df$PC2, na.rm = TRUE)

p3A <- ggplot() +
  geom_segment(
    data = segments_df,
    aes(x = PC1_Soil, y = PC2_Soil, xend = PC1_Gut, yend = PC2_Gut),
    linewidth = 0.30,
    colour = "grey72",
    alpha = 0.80
  ) +
  geom_point(
    data = ord_df,
    aes(x = PC1, y = PC2, fill = Regime, shape = Type),
    size = 2.35,
    colour = "grey25",
    stroke = 0.28,
    alpha = 0.96
  ) +
  annotate(
    "text",
    x = x_rng[1] + 0.03 * diff(x_rng),
    y = y_rng[2] + 0.28 * diff(y_rng),
    label = stat_lab,
    hjust = 0,
    vjust = 1,
    size = 2.18,
    colour = "grey20",
    lineheight = 0.98
  ) +
  scale_fill_manual(values = pal_regime, drop = FALSE) +
  scale_shape_manual(values = c(Soil = 21, Gut = 24), drop = FALSE) +
  labs(
    x = paste0("16S PC1 (", sprintf("%.1f", var_exp[1]), "%)"),
    y = paste0("16S PC2 (", sprintf("%.1f", var_exp[2]), "%)")
  ) +
  coord_cartesian(
    ylim = c(y_rng[1], y_rng[2] + 0.36 * diff(y_rng)),
    clip = "off"
  ) +
  guides(
    fill = guide_legend(
      title = "Treatment",
      order = 1,
      override.aes = list(shape = 21, size = 2.6, colour = "grey25")
    ),
    shape = guide_legend(
      title = "Compartment",
      order = 2,
      override.aes = list(fill = "white", size = 2.6, colour = "grey25")
    )
  ) +
  theme_nature(base_size = 9) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.justification = "center"
  )

# -------------------------
# 11) Fig.3B overall paired soil–gut 16S distance
# -------------------------
dist_median <- median(pair_dist_df$Aitchison_distance, na.rm = TRUE)
dist_n <- nrow(pair_dist_df)

lab_dist <- paste0(
  "n = ", dist_n,
  "\nmedian = ", sprintf("%.1f", dist_median)
)

y_dist_rng <- range(pair_dist_df$Aitchison_distance, na.rm = TRUE)

p3B <- ggplot(pair_dist_df, aes(x = "", y = Aitchison_distance)) +
  geom_violin(
    width = 0.86,
    fill = "grey85",
    colour = NA,
    alpha = 0.90,
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    linewidth = 0.34,
    colour = "grey25",
    fill = "white"
  ) +
  geom_point(
    position = position_jitter(width = 0.07, height = 0),
    size = 1.35,
    alpha = 0.88,
    colour = "grey25"
  ) +
  annotate(
    "text",
    x = 1.34,
    y = y_dist_rng[2],
    label = lab_dist,
    hjust = 0,
    vjust = 1,
    size = 2.45,
    colour = "grey25",
    lineheight = 1.05
  ) +
  labs(
    x = NULL,
    y = "Paired soil–gut\n16S Aitchison distance"
  ) +
  coord_cartesian(clip = "off") +
  theme_nature(base_size = 8.8) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(4, 10, 4, 4, unit = "mm")
  )

# -------------------------
# 12) Fig.3C gut16S_PC1
# -------------------------
dat_pc1 <- readxl::read_excel(supp_xlsx, sheet = "Master_all_data")

need_pc1 <- c("PlotID", "Duration", "Regime", "gut16S_PC1")
miss_pc1 <- setdiff(need_pc1, names(dat_pc1))
if (length(miss_pc1) > 0) {
  stop("Main table missing columns: ", paste(miss_pc1, collapse = ", "))
}

dat_pc1 <- dat_pc1 %>%
  mutate(
    PlotID = trimws(as.character(PlotID)),
    Duration = as.character(Duration),
    Duration = ifelse(Duration == "3y", "5y", Duration),
    Duration = factor(Duration, levels = c("5y", "8y", "10y")),
    Regime = factor(as.character(Regime), levels = c("CK", "NPK", "NPKOM")),
    gut16S_PC1 = as.numeric(gut16S_PC1)
  ) %>%
  filter(is.finite(gut16S_PC1), !is.na(Duration), !is.na(Regime))

fit_pc1 <- lm(gut16S_PC1 ~ Duration * Regime, data = dat_pc1)
tab_pc1 <- car::Anova(fit_pc1, type = 2)

lab_pc1 <- paste0(
  "Type-II ANOVA P\n",
  "Dur=", fmt_p(tab_pc1["Duration", "Pr(>F)"]),
  "  Reg=", fmt_p(tab_pc1["Regime", "Pr(>F)"]),
  "\nDur×Reg=", fmt_p(tab_pc1["Duration:Regime", "Pr(>F)"])
)

y_pc1_rng <- range(dat_pc1$gut16S_PC1, na.rm = TRUE)
y_pc1_top <- y_pc1_rng[2] + diff(y_pc1_rng) * 0.58
y_pc1_label <- y_pc1_rng[2] + diff(y_pc1_rng) * 0.50

lab_pc1 <- paste0(
  "Type-II ANOVA\n",
  "Duration: ", fmt_p_eq(tab_pc1["Duration", "Pr(>F)"]),
  "\nRegime: ", fmt_p_eq(tab_pc1["Regime", "Pr(>F)"]),
  "\nInteraction: ", fmt_p_eq(tab_pc1["Duration:Regime", "Pr(>F)"])
)

lab_pc1_df <- tibble(
  Duration = factor("10y", levels = c("5y", "8y", "10y")),
  Regime = factor("CK", levels = c("CK", "NPK", "NPKOM")),
  y = y_pc1_label,
  label = lab_pc1
)

p3C <- ggplot(dat_pc1, aes(x = Regime, y = gut16S_PC1, fill = Regime)) +
  geom_boxplot(
    width = 0.58,
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
    data = lab_pc1_df,
    aes(x = Regime, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 2.05,
    colour = "grey25",
    lineheight = 1.02
  ) +
  facet_grid(. ~ Duration) +
  scale_x_discrete(labels = c(CK = "CK", NPK = "NPK", NPKOM = "NPKOM")) +
  scale_fill_manual(values = pal_regime, drop = FALSE) +
  labs(x = NULL, y = "Gut 16S PC1") +
  coord_cartesian(ylim = c(y_pc1_rng[1], y_pc1_top), clip = "off") +
  theme_nature(base_size = 8.8) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(4, 5, 4, 6, unit = "mm")
  )

# -------------------------
# 13) Combine
# -------------------------
bottom_row <- cowplot::plot_grid(
  p3B,
  p3C,
  cowplot::ggdraw(),
  ncol = 3,
  rel_widths = c(0.82, 1.12, 0.28),
  labels = c("b", "c", ""),
  label_size = 10,
  label_fontface = "bold",
  label_colour = "grey15"
)

fig3 <- cowplot::plot_grid(
  p3A,
  bottom_row,
  ncol = 1,
  rel_heights = c(1.10, 1.00),
  labels = c("a", ""),
  label_size = 10,
  label_fontface = "bold",
  label_colour = "grey15"
)

print(fig3)

# -------------------------
# 14) Save
# -------------------------
ggsave(
  filename = file.path(out_dir, "Fig3_gut_filtering_revised_polished.svg"),
  plot = fig3,
  width = 170,
  height = 140,
  units = "mm",
  device = svg_device,
  bg = "white"
)

ggsave(
  filename = file.path(out_dir, "Fig3_gut_filtering_revised_polished.png"),
  plot = fig3,
  width = 170,
  height = 140,
  units = "mm",
  dpi = 600,
  bg = "white"
)

write.csv(
  pair_dist_df,
  file = file.path(out_dir, "Fig3_paired_16S_Aitchison_distance.csv"),
  row.names = FALSE
)

sink(file.path(out_dir, "Fig3_16S_filtering_stats.txt"))
cat("Fig.3 statistics\n")
cat("================\n\n")
cat("PERMANOVA: Aitchison distance, strata = PairID\n")
print(perm)
cat("\nDispersion test\n")
print(bd_perm)
cat("\nPaired soil-gut 16S Aitchison distance\n")
cat("n =", dist_n, "\n")
cat("median =", sprintf("%.3f", dist_median), "\n")
cat("\ngut16S_PC1 model: gut16S_PC1 ~ Duration * Regime\n")
print(tab_pc1)
sink()

cat("Saved revised Fig.3 to: ", out_dir, "\n")
