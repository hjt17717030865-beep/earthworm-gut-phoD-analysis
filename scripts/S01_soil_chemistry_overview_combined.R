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
# Supplementary Figure S1: soil chemistry overview
# Panel a: soil chemistry distributions by fertilization regime and duration.
# Panel b: Spearman correlation matrix among soil chemistry variables.
# Panel c: ChemPC1 loadings from scaled soil chemistry PCA.
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(car)
  library(multcompView)
  library(svglite)
  library(grid)
  library(tibble)
  library(writexl)
  library(patchwork)
})

f_sum <- DATA_XLSX
out_dir <- OUTPUT_DIR
dat_sum <- read_excel(f_sum, sheet = "Master_all_data")

chem_vars <- c(
  "pH",
  "water content(%)",
  "TOC(%)",
  "NH4+(mg/kg)",
  "NO3-(mg/kg)",
  "TN(mg/kg)",
  "Olsen-P(mg/kg)",
  "TP(mg/kg)"
)

chem_labels <- c(
  "pH",
  "Water content (%)",
  "TOC (%)",
  "NH4+ (mg kg-1)",
  "NO3- (mg kg-1)",
  "TN (mg kg-1)",
  "Olsen-P (mg kg-1)",
  "TP (mg kg-1)"
)
names(chem_labels) <- chem_vars

miss_chem <- setdiff(c("Duration", "Regime", chem_vars), names(dat_sum))
if (length(miss_chem) > 0) {
  stop("Missing columns in Master_all_data: ", paste(miss_chem, collapse = ", "))
}

compress_letters <- function(letter_vec, alphabet = letters) {
  all_chars <- sort(unique(unlist(strsplit(paste0(letter_vec, collapse = ""), ""))))
  map <- setNames(alphabet[seq_along(all_chars)], all_chars)

  out <- vapply(letter_vec, function(s) {
    chars <- strsplit(s, "")[[1]]
    new_chars <- unname(map[chars])
    paste0(sort(new_chars), collapse = "")
  }, character(1))

  names(out) <- names(letter_vec)
  out
}

relabel_letters_by_mean_desc <- function(letter_vec, mean_vec, alphabet = letters) {
  stopifnot(all(names(letter_vec) %in% names(mean_vec)))

  groups <- sort(unique(unlist(strsplit(paste0(letter_vec, collapse = ""), ""))))

  group_score <- sapply(groups, function(g) {
    regs <- names(letter_vec)[grepl(g, letter_vec, fixed = TRUE)]
    max(mean_vec[regs], na.rm = TRUE)
  })

  groups_sorted <- names(sort(group_score, decreasing = TRUE))
  new_labels <- alphabet[seq_along(groups_sorted)]
  names(new_labels) <- groups_sorted

  new_letter_vec <- sapply(letter_vec, function(s) {
    parts <- strsplit(s, "")[[1]]
    parts_new <- unname(new_labels[parts])
    paste0(sort(parts_new), collapse = "")
  })

  new_letter_vec <- compress_letters(new_letter_vec, alphabet = alphabet)
  names(new_letter_vec) <- names(letter_vec)
  new_letter_vec
}

get_tukey_letters <- function(df, response, group_var, alphabet = letters) {
  mean_by_group <- tapply(df[[response]], df[[group_var]], mean, na.rm = TRUE)

  fit_w <- aov(as.formula(paste(response, "~", group_var)), data = df)
  tk <- TukeyHSD(fit_w, group_var)[[1]]

  letters_raw <- multcompView::multcompLetters(tk[, "p adj"])$Letters
  letters_sorted <- relabel_letters_by_mean_desc(
    letter_vec = letters_raw,
    mean_vec = mean_by_group,
    alphabet = alphabet
  )

  out <- tibble(
    group = names(letters_sorted),
    tukey_letter = unname(letters_sorted)
  )
  names(out)[1] <- group_var
  out
}

pal_regime <- c(
  "CK" = "#BFBFBF",
  "NPK" = "#7FAFB4",
  "NPKOM" = "#D39A92"
)

theme_s1 <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size, colour = "grey20"),
      axis.text = element_text(size = base_size, colour = "grey25"),
      axis.line = element_line(linewidth = 0.28, colour = "grey25"),
      axis.ticks = element_line(linewidth = 0.28, colour = "grey25"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, colour = "grey20"),
      strip.placement = "outside",
      legend.title = element_text(size = base_size - 1, colour = "grey20"),
      legend.text = element_text(size = base_size - 1, colour = "grey20"),
      panel.spacing = unit(0.35, "lines"),
      plot.margin = margin(3, 3, 3, 3, unit = "mm")
    )
}

# Panel a
dat_chem <- dat_sum %>%
  transmute(
    Duration = as.character(Duration),
    Regime = as.character(Regime),
    across(all_of(chem_vars), as.numeric)
  ) %>%
  mutate(
    Duration = trimws(Duration),
    Duration = gsub("\\s+", "", Duration),
    Duration = gsub("years|Years|yr|YR|Year|year", "y", Duration),
    Duration = gsub("^3y$", "5y", Duration),
    Duration = gsub("^5$", "5y", Duration),
    Duration = gsub("^8$", "8y", Duration),
    Duration = gsub("^10$", "10y", Duration),
    Duration = factor(Duration, levels = c("5y", "8y", "10y")),
    Regime = factor(Regime, levels = c("CK", "NPK", "NPKOM"))
  ) %>%
  pivot_longer(
    cols = all_of(chem_vars),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  filter(!is.na(Duration), !is.na(Regime), !is.na(Value)) %>%
  mutate(
    Variable = factor(Variable, levels = chem_vars)
  )

letters_lower <- dat_chem %>%
  group_by(Variable, Duration) %>%
  group_modify(~{
    out <- get_tukey_letters(.x, "Value", "Regime", letters)
    names(out)[names(out) == "tukey_letter"] <- "letter_lower"
    out
  }) %>%
  ungroup() %>%
  mutate(Regime = factor(Regime, levels = c("CK", "NPK", "NPKOM")))

letters_upper <- dat_chem %>%
  group_by(Variable, Regime) %>%
  group_modify(~{
    out <- get_tukey_letters(.x, "Value", "Duration", LETTERS)
    names(out)[names(out) == "tukey_letter"] <- "letter_upper"
    out
  }) %>%
  ungroup() %>%
  mutate(Duration = factor(Duration, levels = c("5y", "8y", "10y")))

ypos_df <- dat_chem %>%
  group_by(Variable, Duration, Regime) %>%
  summarise(
    ymax = max(Value, na.rm = TRUE),
    ymin = min(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Variable) %>%
  mutate(
    var_range = max(ymax, na.rm = TRUE) - min(ymin, na.rm = TRUE),
    pad = ifelse(is.finite(var_range) & var_range > 0, 0.08 * var_range, 0.05),
    y_lab = ymax + pad
  ) %>%
  ungroup()

letters_plot <- ypos_df %>%
  left_join(letters_lower, by = c("Variable", "Duration", "Regime")) %>%
  left_join(letters_upper, by = c("Variable", "Duration", "Regime")) %>%
  mutate(label = paste0(letter_upper, "\u2009", letter_lower))

p_S1a <- ggplot(dat_chem, aes(x = Regime, y = Value)) +
  geom_boxplot(
    aes(fill = Regime),
    width = 0.60,
    outlier.shape = NA,
    linewidth = 0.26,
    colour = "grey30",
    alpha = 0.18
  ) +
  geom_point(
    aes(colour = Regime),
    position = position_jitter(width = 0.09, height = 0),
    size = 0.72,
    alpha = 0.78
  ) +
  geom_text(
    data = letters_plot,
    aes(x = Regime, y = y_lab, label = label),
    inherit.aes = FALSE,
    size = 2.35,
    vjust = 0,
    colour = "grey30"
  ) +
  facet_grid(Variable ~ Duration, scales = "free_y") +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.16))) +
  scale_colour_manual(values = pal_regime, drop = FALSE) +
  scale_fill_manual(values = pal_regime, drop = FALSE) +
  theme_s1(base_size = 8.5) +
  theme(legend.position = "none")

# Panel b
chem_mat <- dat_sum %>%
  select(all_of(chem_vars)) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame()

n_var <- length(chem_vars)
cor_mat <- matrix(NA_real_, n_var, n_var)
p_mat <- matrix(NA_real_, n_var, n_var)
rownames(cor_mat) <- colnames(cor_mat) <- chem_vars
rownames(p_mat) <- colnames(p_mat) <- chem_vars

for (i in seq_len(n_var)) {
  for (j in seq_len(n_var)) {
    x <- chem_mat[[i]]
    y <- chem_mat[[j]]
    ok <- is.finite(x) & is.finite(y)

    if (sum(ok) >= 3) {
      ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
      cor_mat[i, j] <- unname(ct$estimate)
      p_mat[i, j] <- ct$p.value
    }
  }
}

cor_df <- expand.grid(
  Var_y = chem_vars,
  Var_x = chem_vars,
  stringsAsFactors = FALSE
) %>%
  mutate(
    i = match(Var_y, chem_vars),
    j = match(Var_x, chem_vars),
    r = mapply(function(a, b) cor_mat[a, b], i, j),
    p = mapply(function(a, b) p_mat[a, b], i, j)
  ) %>%
  filter(i > j) %>%
  mutate(
    star = case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ ""
    ),
    label = paste0(sprintf("%.2f", r), star),
    Var_x_label = chem_labels[Var_x],
    Var_y_label = chem_labels[Var_y]
  )

x_levels <- chem_labels[chem_vars[-length(chem_vars)]]
y_levels <- rev(chem_labels[chem_vars[-1]])

cor_df <- cor_df %>%
  mutate(
    Var_x_label = factor(Var_x_label, levels = x_levels),
    Var_y_label = factor(Var_y_label, levels = y_levels)
  )

top_label_df <- tibble(
  Var_x_label = factor(x_levels, levels = x_levels),
  x_num = seq_along(x_levels),
  y_num = seq(from = length(y_levels) + 0.55, by = -1, length.out = length(x_levels)),
  label = x_levels
)

p_S1b <- ggplot(cor_df, aes(x = Var_x_label, y = Var_y_label, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = label), size = 2.55, colour = "grey20", fontface = "bold") +
  geom_text(
    data = top_label_df,
    aes(x = x_num, y = y_num, label = label),
    inherit.aes = FALSE,
    angle = 45,
    hjust = 0,
    vjust = 0,
    size = 2.55,
    colour = "grey20"
  ) +
  scale_fill_gradient2(
    low = "#D39A92",
    mid = "white",
    high = "#7FAFB4",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = seq(-1, 1, by = 0.5),
    name = "Spearman r"
  ) +
  coord_fixed(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_s1(base_size = 8.2) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour = "grey20", size = 7.8),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(13, "mm"),
    plot.margin = margin(10, 4, 4, 4, unit = "mm")
  )

# Panel c
pca_chem <- prcomp(chem_mat, center = TRUE, scale. = TRUE)
chem_var_explained <- (pca_chem$sdev^2) / sum(pca_chem$sdev^2)
chem_pc1_load <- pca_chem$rotation[, 1]
chem_pc1_contrib <- (chem_pc1_load^2) / sum(chem_pc1_load^2)

chem_pc1_loadings <- tibble(
  Variable = names(chem_pc1_load),
  Loading_PC1 = as.numeric(chem_pc1_load),
  Abs_loading_PC1 = abs(as.numeric(chem_pc1_load)),
  Contribution_PC1 = as.numeric(chem_pc1_contrib),
  Contribution_PC1_percent = 100 * as.numeric(chem_pc1_contrib)
) %>%
  arrange(Loading_PC1)

chem_pca_variance <- tibble(
  PC = paste0("PC", seq_along(chem_var_explained)),
  Eigenvalue = pca_chem$sdev^2,
  Proportion = chem_var_explained,
  Percent = 100 * chem_var_explained,
  Cumulative = cumsum(chem_var_explained),
  Cumulative_percent = 100 * cumsum(chem_var_explained)
)

p_S1c <- ggplot(
  chem_pc1_loadings,
  aes(x = reorder(Variable, Loading_PC1), y = Loading_PC1)
) +
  geom_col(
    fill = "#7FAFB4",
    colour = "grey30",
    linewidth = 0.26,
    width = 0.62,
    alpha = 0.85
  ) +
  coord_flip() +
  labs(x = NULL, y = "ChemPC1 loading") +
  theme_s1(base_size = 8.5) +
  theme(legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.04)))

fig_S1 <- p_S1a + p_S1b + p_S1c +
  plot_layout(
    design = "
AAA
BBC
",
    heights = c(1.9, 1),
    widths = c(1.15, 1.15, 0.9)
  ) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold", size = 11, colour = "black"),
    plot.tag.position = c(0.01, 0.99),
    plot.margin = margin(2, 2, 2, 2, unit = "mm")
  )

ggsave(
  filename = file.path(out_dir, "FigS1_soil_chemistry_overview_combined.svg"),
  plot = fig_S1,
  width = 180,
  height = 285,
  units = "mm",
  device = svglite::svglite,
  bg = "white"
)

ggsave(
  filename = file.path(out_dir, "FigS1_soil_chemistry_overview_combined.png"),
  plot = fig_S1,
  width = 180,
  height = 285,
  units = "mm",
  dpi = 600,
  bg = "white"
)

write_xlsx(
  list(
    chem_pca_variance = chem_pca_variance,
    chem_pc1_loadings = chem_pc1_loadings
  ),
  file.path(out_dir, "ChemPC1_loadings.xlsx")
)

cat("Saved combined Supplementary Figure S1 PNG + SVG to:", out_dir, "\n")
