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
# Figure S8 | Robustness of the AI–P-index association
#
# One-script reproducible version
#
# Output:
#   Figure_S7_Pindex_robustness.svg
#   Figure_S7_Pindex_robustness.png
#   Figure_S7_Pindex_robustness_stats.csv
#   Figure_S7_reproducible_objects.rds
# ============================================================

# -------------------------
# 0) Packages
# -------------------------
pkg_needed <- c(
  "readxl", "dplyr", "tidyr", "stringr", "tibble",
  "ggplot2", "svglite", "readr"
)

pkg_to_install <- pkg_needed[!sapply(pkg_needed, requireNamespace, quietly = TRUE)]
if (length(pkg_to_install) > 0) install.packages(pkg_to_install)

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(svglite)
  library(readr)
})

# -------------------------
# 1) Paths
# -------------------------
file_suppl <- DATA_XLSX
out_dir <- OUTPUT_DIR
if (!file.exists(file_suppl)) {
  stop("Supplementary Data file not found: ", file_suppl)
}

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# -------------------------
# 2) Read P-gene wide table
# -------------------------
dat_raw <- read_excel(
  file_suppl,
  sheet = "P_gene_qPCR_matrix"
)
if (!"gene" %in% names(dat_raw)) {
  stop("Input P-gene table must contain a column named 'gene'.")
}

# -------------------------
# 3) Wide to long table
# -------------------------
dat_long <- dat_raw %>%
  pivot_longer(
    cols = -gene,
    names_to = "SampleID",
    values_to = "copies"
  ) %>%
  mutate(
    gene = as.character(gene),
    SampleID = as.character(SampleID),
    copies = suppressWarnings(as.numeric(copies)),
    
    Duration = str_extract(SampleID, "^\\d+y"),
    
    Habitat = dplyr::case_when(
      str_detect(SampleID, "^\\d+yEB") ~ "gut",
      str_detect(SampleID, "^\\d+yE")  ~ "gut",
      TRUE ~ "soil"
    ),
    
    Taxon = dplyr::case_when(
      str_detect(SampleID, "^\\d+yEB") ~ "EB",
      str_detect(SampleID, "^\\d+yE")  ~ "E",
      TRUE ~ "soil"
    ),
    
    Regime = dplyr::case_when(
      str_detect(SampleID, "NPKOM") ~ "NPKOM",
      str_detect(SampleID, "NPK")   ~ "NPK",
      str_detect(SampleID, "CK")    ~ "CK",
      TRUE ~ NA_character_
    ),
    
    Rep = suppressWarnings(as.numeric(str_extract(SampleID, "(?<=-)\\d+$"))),
    
    PlotID = SampleID %>%
      str_replace("^([0-9]+y)EB", "\\1") %>%
      str_replace("^([0-9]+y)E", "\\1"),
    
    gene = str_to_lower(gene),
    gene = ifelse(gene == "ppk3", "ppk", gene)
  ) %>%
  dplyr::select(SampleID, PlotID, Duration, Habitat, Taxon, Regime, Rep, gene, copies)

# -------------------------
# 4) Scaling helpers
# -------------------------
safe_z <- function(x) {
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (sum(ok) >= 2 && sd(x[ok], na.rm = TRUE) > 0) {
    out[ok] <- as.numeric(scale(x[ok]))
  }
  out
}

safe_log1p_z <- function(x) {
  safe_z(log10(x + 1))
}

safe_minmax <- function(x) {
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (sum(ok) >= 2) {
    rng <- range(x[ok], na.rm = TRUE)
    if ((rng[2] - rng[1]) > 0) {
      out[ok] <- (x[ok] - rng[1]) / (rng[2] - rng[1])
    }
  }
  out
}

# -------------------------
# 5) General P-index calculator
#
# Important:
#   This follows your original P-index logic:
#   gene-wise scaling is calculated across all selected samples first,
#   then E-gut samples are used for AI association.
# -------------------------
calc_pindex_general <- function(df_long,
                                genes_use,
                                gene_set_name,
                                scaling = c("z", "log1p_z", "minmax")) {
  
  scaling <- match.arg(scaling)
  
  df_use <- df_long %>%
    filter(gene %in% tolower(genes_use))
  
  if (nrow(df_use) == 0) {
    stop("No rows found for gene set: ", gene_set_name)
  }
  
  missing_genes <- setdiff(tolower(genes_use), unique(df_use$gene))
  if (length(missing_genes) > 0) {
    warning(
      "Missing genes for ", gene_set_name, ": ",
      paste(missing_genes, collapse = ", ")
    )
  }
  
  # Gene-wise scaling across all selected samples
  df_scaled <- df_use %>%
    group_by(gene) %>%
    mutate(
      scaled_value = dplyr::case_when(
        scaling == "z"       ~ safe_z(copies),
        scaling == "log1p_z" ~ safe_log1p_z(copies),
        scaling == "minmax"  ~ safe_minmax(copies)
      )
    ) %>%
    ungroup()
  
  # Use E-gut samples only for AI association
  df_sample <- df_scaled %>%
    filter(
      Habitat == "gut",
      str_detect(SampleID, "^\\d+yE"),
      !str_detect(SampleID, "^\\d+yEB")
    ) %>%
    group_by(SampleID, PlotID, Duration, Habitat, Taxon, Regime, Rep) %>%
    summarise(
      P_index = mean(scaled_value, na.rm = TRUE),
      n_selected_genes = n_distinct(gene),
      n_valid_genes = sum(is.finite(scaled_value)),
      usable_for_Pindex = n_valid_genes > 0,
      .groups = "drop"
    ) %>%
    mutate(
      gene_set_name = gene_set_name,
      scaling = scaling
    ) %>%
    arrange(Duration, Regime, Rep, SampleID)
  
  df_sample
}

# -------------------------
# 6) Define robustness versions
# -------------------------
gene_sets <- list(
  main_6  = c("phnk", "phox", "gcd", "pqqc", "ppk", "ppx"),
  no_gcd  = c("phnk", "phox", "pqqc", "ppk", "ppx"),
  no_pqqc = c("phnk", "phox", "gcd", "ppk", "ppx")
)

scalings <- c("z", "log1p_z", "minmax")

pindex_robust <- bind_rows(
  lapply(names(gene_sets), function(gs) {
    bind_rows(
      lapply(scalings, function(sc) {
        calc_pindex_general(
          df_long = dat_long,
          genes_use = gene_sets[[gs]],
          gene_set_name = gs,
          scaling = sc
        )
      })
    )
  })
)

cat("\nP-index robustness table dimension:\n")
print(dim(pindex_robust))

cat("\nVersion counts:\n")
print(table(pindex_robust$gene_set_name, pindex_robust$scaling))

# -------------------------
# 7) Read main table and merge AI
# -------------------------
dat_main <- read_excel(
  file_suppl,
  sheet = "Master_all_data"
) %>%
  mutate(
    PlotID = as.character(PlotID),
    AI = as.numeric(AI)
  ) %>%
  dplyr::select(PlotID, AI)

pindex_ai <- pindex_robust %>%
  left_join(dat_main, by = "PlotID")

# -------------------------
# 8) Fit AI ~ P_index for each robustness version
# -------------------------
fit_one <- function(df) {
  fit <- lm(AI ~ P_index, data = df)
  sm <- summary(fit)
  co <- sm$coefficients
  
  tibble(
    n = nrow(df),
    intercept = unname(co["(Intercept)", "Estimate"]),
    slope = unname(co["P_index", "Estimate"]),
    slope_se = unname(co["P_index", "Std. Error"]),
    t_value = unname(co["P_index", "t value"]),
    p_value = unname(co["P_index", "Pr(>|t|)"]),
    r2 = unname(sm$r.squared),
    adj_r2 = unname(sm$adj.r.squared)
  )
}

s8_model_res <- pindex_ai %>%
  group_by(gene_set_name, scaling) %>%
  group_modify(~ fit_one(.x)) %>%
  ungroup() %>%
  arrange(gene_set_name, scaling)

cat("\nS8 model results:\n")
print(s8_model_res)

# -------------------------
# 9) Prepare plotting table
# -------------------------
s8_model_res2 <- s8_model_res %>%
  mutate(
    ci_low  = slope - 1.96 * slope_se,
    ci_high = slope + 1.96 * slope_se,
    
    gene_set_name = factor(
      gene_set_name,
      levels = c("main_6", "no_gcd", "no_pqqc")
    ),
    
    scaling = factor(
      scaling,
      levels = c("z", "log1p_z", "minmax")
    ),
    
    gene_set_label = dplyr::case_when(
      gene_set_name == "main_6"  ~ "Main 6-gene set",
      gene_set_name == "no_gcd"  ~ "Without gcd",
      gene_set_name == "no_pqqc" ~ "Without pqqC"
    ),
    
    scaling_label = dplyr::case_when(
      scaling == "z"       ~ "Z-score",
      scaling == "log1p_z" ~ "log10(x+1) + Z-score",
      scaling == "minmax"  ~ "Min-max"
    ),
    
    panel_label = paste(gene_set_label, scaling_label, sep = " | ")
  ) %>%
  arrange(gene_set_name, scaling) %>%
  mutate(
    panel_label = factor(
      panel_label,
      levels = rev(c(
        "Main 6-gene set | Z-score",
        "Main 6-gene set | log10(x+1) + Z-score",
        "Main 6-gene set | Min-max",
        "Without gcd | Z-score",
        "Without gcd | log10(x+1) + Z-score",
        "Without gcd | Min-max",
        "Without pqqC | Z-score",
        "Without pqqC | log10(x+1) + Z-score",
        "Without pqqC | Min-max"
      ))
    )
  )

# -------------------------
# 10) Long table for one-ggplot faceted layout
# -------------------------
plot_slope <- s8_model_res2 %>%
  transmute(
    Panel = "Slope ± 95% CI",
    panel_label,
    scaling_label,
    x = slope,
    xmin = ci_low,
    xmax = ci_high
  )

plot_r2 <- s8_model_res2 %>%
  transmute(
    Panel = "Model fit",
    panel_label,
    scaling_label,
    x = r2
  )

plot_slope$Panel <- factor(plot_slope$Panel, levels = c("Slope ± 95% CI", "Model fit"))
plot_r2$Panel    <- factor(plot_r2$Panel,    levels = c("Slope ± 95% CI", "Model fit"))

# -------------------------
# 11) Theme and colours
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
      panel.spacing.x = unit(8, "mm"),
      axis.ticks.length = unit(2, "pt"),
      plot.margin = margin(4, 4, 4, 4, unit = "mm")
    )
}

pal_scaling <- c(
  "Z-score" = "#AEB6CC",
  "log10(x+1) + Z-score" = "#B2CBB2",
  "Min-max" = "#CBB3B3"
)

# -------------------------
# 12) Final S8 plot
# -------------------------
fig_s8 <- ggplot() +
  geom_vline(
    data = tibble(Panel = factor("Slope ± 95% CI", levels = c("Slope ± 95% CI", "Model fit"))),
    aes(xintercept = 0),
    linetype = 2,
    linewidth = 0.30,
    colour = "grey50"
  ) +
  geom_segment(
    data = plot_slope,
    aes(
      x = xmin,
      xend = xmax,
      y = panel_label,
      yend = panel_label,
      colour = scaling_label
    ),
    linewidth = 0.45
  ) +
  geom_point(
    data = plot_slope,
    aes(
      x = x,
      y = panel_label,
      colour = scaling_label
    ),
    size = 2.0
  ) +
  geom_col(
    data = plot_r2,
    aes(
      x = x,
      y = panel_label,
      fill = scaling_label
    ),
    width = 0.62,
    colour = "grey25",
    linewidth = 0.25,
    orientation = "y"
  ) +
  scale_colour_manual(values = pal_scaling) +
  scale_fill_manual(values = pal_scaling) +
  facet_grid(
    . ~ Panel,
    scales = "free_x"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_nature(base_size = 9) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8.8, colour = "grey25")
  )

print(fig_s8)

# -------------------------
# 13) Save outputs
# -------------------------
ggsave(
  filename = file.path(out_dir, "Figure_S7_Robustness_of_the_association_between_phoD_amplification_and_P-acquisition_potential.svg"),
  plot = fig_s8,
  width = 180,
  height = 105,
  units = "mm",
  device = svglite::svglite
)

ggsave(
  filename = file.path(out_dir, "Figure_S7_Robustness_of_the_association_between_phoD_amplification_and_P-acquisition_potential.png"),
  plot = fig_s8,
  width = 180,
  height = 105,
  units = "mm",
  dpi = 600
)

write_csv(
  s8_model_res2,
  file = file.path(out_dir, "Figure_S7_Pindex_robustness_stats.csv")
)

saveRDS(
  list(
    dat_long = dat_long,
    pindex_robust = pindex_robust,
    pindex_ai = pindex_ai,
    s8_model_res = s8_model_res,
    s8_model_res2 = s8_model_res2,
    fig_s8 = fig_s8
  ),
  file = file.path(out_dir, "Figure_S7_reproducible_objects.rds")
)

cat("\nSaved S8 outputs to:\n")
cat(out_dir, "\n")
