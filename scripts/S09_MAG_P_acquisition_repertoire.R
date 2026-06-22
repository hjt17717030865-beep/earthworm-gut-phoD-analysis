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
# Figure S9 | Complete genome-resolved P-acquisition repertoire
# Final logic:
#   A = complete module-level UpSet
#   B = MAG × P-gene presence/absence matrix
#
# This version avoids soil-vs-gut differential interpretation.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
dir_meta <- CODE_ROOT
out_dir <- OUTPUT_DIR
file_suppl <- DATA_XLSX
# ------------------------------------------------------------
# Read MAG module source table
# ------------------------------------------------------------
mag_detail <- readxl::read_excel(
  file_suppl,
  sheet = "P_MAG_summary_tax_compartment"
) %>%
  mutate(
    MAG = as.character(MAG),
    compartment = str_to_lower(as.character(compartment)),
    KO_list = replace_na(as.character(KO_list), "")
  ) %>%
  filter(compartment %in% c("soil", "gut"))

# ------------------------------------------------------------
# Read gene presence table
# ------------------------------------------------------------
mag_gene <- readxl::read_excel(
  file_suppl,
  sheet = "P_MAG_gene_presence"
) %>%
  mutate(
    MAG = as.character(MAG),
    compartment = str_to_lower(as.character(compartment)),
    P_score = as.integer(P_score)
  ) %>%
  filter(compartment %in% c("soil", "gut"))

gene_sets <- c("phoD", "phoX", "pqqC", "ppk", "ppx", "phnK")

miss_gene_cols <- setdiff(c("MAG", "compartment", gene_sets, "P_score"), names(mag_gene))
if (length(miss_gene_cols) > 0) {
  stop("Missing columns in P_MAG_gene_presence.tsv: ", paste(miss_gene_cols, collapse = ", "))
}

# ------------------------------------------------------------
# Style
# ------------------------------------------------------------
nature_col  <- "#3A5A78"
nature_grey <- "#D9DDE1"

# =========================
# Step now | fix layout only
# Replace the corresponding parts in your script
# =========================

theme_s11 <- function(base_size = 8.8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size, colour = "grey20"),
      axis.text  = element_text(size = base_size - 0.2, colour = "grey25"),
      axis.line  = element_line(linewidth = 0.30, colour = "grey25"),
      axis.ticks = element_line(linewidth = 0.30, colour = "grey25"),
      plot.title = element_text(
        size = base_size + 0.2,
        colour = "grey20",
        face = "plain",
        hjust = 0,
        margin = margin(b = 2, unit = "mm")
      ),
      plot.margin = margin(2, 3, 2, 3, unit = "mm")
    )
}
# ============================================================
# Helper | Manual UpSet
# ============================================================

make_manual_upset <- function(dat, set_cols, set_labels, title_text) {
  
  # Keep rows carrying at least one selected set
  dat2 <- dat %>%
    mutate(set_sum = rowSums(across(all_of(set_cols)), na.rm = TRUE)) %>%
    filter(set_sum > 0)
  
  # Set sizes
  set_counts <- dat2 %>%
    summarise(across(all_of(set_cols), ~sum(. == 1, na.rm = TRUE))) %>%
    pivot_longer(cols = everything(), names_to = "set", values_to = "set_size")
  
  set_order <- set_counts %>%
    arrange(desc(set_size), set) %>%
    pull(set)
  
  set_levels_y <- rev(set_order)
  
  # Exclusive intersections
  combo_df <- dat2 %>%
    mutate(
      combo_key = apply(
        select(., all_of(set_cols)),
        1,
        function(x) paste(set_cols[as.integer(x) == 1], collapse = "+")
      )
    ) %>%
    filter(combo_key != "") %>%
    count(combo_key, name = "intersection_size") %>%
    arrange(desc(intersection_size), combo_key) %>%
    mutate(
      intersection_id = paste0("I", row_number()),
      intersection_id = factor(intersection_id, levels = intersection_id)
    )
  
  # Matrix data
  matrix_df <- combo_df %>%
    select(combo_key, intersection_id) %>%
    crossing(set = set_order) %>%
    mutate(
      active = map2_lgl(combo_key, set, ~ .y %in% str_split(.x, "\\+")[[1]]),
      set = factor(set, levels = set_levels_y),
      y_pos = match(as.character(set), set_levels_y)
    )
  
  # Connecting lines
  segment_df <- matrix_df %>%
    filter(active) %>%
    group_by(intersection_id) %>%
    summarise(
      y_min = min(y_pos),
      y_max = max(y_pos),
      n_active = n(),
      .groups = "drop"
    ) %>%
    filter(n_active >= 2)
  
  # Top bars
  p_top <- ggplot(combo_df, aes(x = intersection_id, y = intersection_size)) +
    geom_col(width = 0.82, fill = nature_col) +
    geom_text(
      aes(label = intersection_size),
      vjust = -0.25,
      size = 2.5,
      colour = "grey20"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL, y = "Number of MAGs") +
    theme_s11(base_size = 8.4) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey92"),
      panel.grid.major.x = element_blank()
    )
  
  # Left set-size bars
  set_df <- set_counts %>%
    mutate(set = factor(set, levels = set_levels_y))
  
  p_set <- ggplot(set_df, aes(x = set_size, y = set)) +
    geom_col(width = 0.62, fill = nature_col) +
    scale_x_reverse(expand = expansion(mult = c(0.05, 0.08))) +
    labs(x = "MAGs with module", y = NULL) +
    theme_s11(base_size = 8.4) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey92"),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(0, 0, 8, 10, unit = "mm")
    )
  # Background stripes
  stripe_df <- tibble(y_pos = seq_along(set_levels_y)) %>%
    filter(y_pos %% 2 == 0)
  
  # Matrix
  p_matrix <- ggplot(matrix_df, aes(x = intersection_id, y = y_pos)) +
    geom_rect(
      data = stripe_df,
      aes(
        xmin = -Inf, xmax = Inf,
        ymin = y_pos - 0.43, ymax = y_pos + 0.43
      ),
      inherit.aes = FALSE,
      fill = "grey95",
      colour = NA
    ) +
    geom_segment(
      data = segment_df,
      aes(
        x = intersection_id, xend = intersection_id,
        y = y_min, yend = y_max
      ),
      inherit.aes = FALSE,
      linewidth = 0.70,
      colour = nature_col,
      lineend = "round"
    ) +
    geom_point(
      aes(fill = active, colour = active),
      shape = 21,
      size = 2.,
      stroke = 0.25
    ) +
    scale_fill_manual(
      values = c("FALSE" = nature_grey, "TRUE" = nature_col),
      guide = "none"
    ) +
    scale_colour_manual(
      values = c("FALSE" = nature_grey, "TRUE" = nature_col),
      guide = "none"
    ) +
    scale_y_continuous(
      breaks = seq_along(set_levels_y),
      labels = set_labels[set_levels_y],
      expand = expansion(mult = c(0.04, 0.04))
    ) +
    labs(x = NULL, y = NULL) +
    theme_s11(base_size = 8.7) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 8.3, colour = "grey20"),
      panel.grid.major.x = element_line(linewidth = 0.20, colour = "grey94"),
      panel.grid.major.y = element_blank()
    )
  
  library(cowplot)
  
  right_col <- cowplot::plot_grid(
    p_top    + theme(plot.margin = margin(2, 2, 4, 8, unit = "mm")),
    p_matrix + theme(plot.margin = margin(0, 2, 2, 0, unit = "mm")),
    ncol = 1, align = "v", axis = "lr",
    rel_heights = c(1, 1)
  )
  
  left_col <- cowplot::plot_grid(
    ggplot() + theme_void(),
    p_set + theme(plot.margin = margin(0, 0, 2, 2, unit = "mm")),
    ncol = 1, align = "v", axis = "lr",
    rel_heights = c(1, 1)
  )
  
  body <- cowplot::plot_grid(
    left_col, right_col,
    ncol = 2,
    rel_widths = c(0.45, 1)
    )
  
  p_out <- body
  
  return(p_out)
}
# ============================================================
# Panel A | Complete module-level UpSet
# ============================================================

module_cols <- c(
  "P_transport",
  "PolyP_turnover",
  "Organic_P_mineralization",
  "Phosphonate",
  "Inorganic_P_solubilization"
)

module_labels <- c(
  "P_transport" = "P transport",
  "PolyP_turnover" = "PolyP turnover",
  "Organic_P_mineralization" = "Organic P mineralization",
  "Phosphonate" = "Phosphonate",
  "Inorganic_P_solubilization" = "Inorganic P solubilization"
)

mag_module <- mag_detail %>%
  transmute(
    MAG,
    compartment,
    
    # phosphate transport / uptake
    P_transport = as.integer(str_detect(KO_list, "K02036|K02037|K02038|K02040|K02033")),
    
    # polyphosphate turnover
    PolyP_turnover = as.integer(str_detect(KO_list, "K00937|K01507")),
    
    # organic P mineralization
    Organic_P_mineralization = as.integer(str_detect(KO_list, "K01113|K01093|K01077")),
    
    # phosphonate metabolism
    Phosphonate = as.integer(str_detect(KO_list, "K06164|K06165|K05780")),
    
    # inorganic P solubilization / pqq-related
    Inorganic_P_solubilization = as.integer(str_detect(KO_list, "K00117|K06130"))
  )

# replace Panel A title with a shorter one

pA <- make_manual_upset(
  dat = mag_module,
  set_cols = module_cols,
  set_labels = module_labels,
  title_text = "A | Module co-occurrence in recovered MAGs"
)

# ============================================================
# Panel B | MAG × P-gene presence/absence matrix
# ============================================================

gene_labels <- c(
  "phoD" = "phoD",
  "phoX" = "phoX",
  "pqqC" = "pqqC",
  "ppk"  = "ppk",
  "ppx"  = "ppx",
  "phnK" = "phnK"
)

gene_dat <- mag_gene %>%
  mutate(
    P_score = as.integer(P_score),
    combo_key = apply(
      select(., all_of(gene_sets)),
      1,
      function(x) paste(gene_sets[as.integer(x) == 1], collapse = "+")
    )
  ) %>%
  filter(P_score > 0) %>%
  arrange(desc(P_score), combo_key, MAG) %>%
  mutate(
    MAG_order = paste0("MAG_", row_number()),
    MAG_order = factor(MAG_order, levels = rev(MAG_order))
  )

gene_long <- gene_dat %>%
  select(MAG_order, all_of(gene_sets)) %>%
  pivot_longer(
    cols = all_of(gene_sets),
    names_to = "gene",
    values_to = "present"
  ) %>%
  mutate(
    gene = factor(gene, levels = gene_sets),
    present = factor(present, levels = c(0, 1))
  )

# replace Panel B title with a shorter one

pB <- ggplot(gene_long, aes(x = gene, y = MAG_order, fill = present)) +
  geom_tile(colour = "white", linewidth = 0.20) +
  scale_fill_manual(
    values = c("0" = nature_grey, "1" = nature_col),
    labels = c("Absent", "Present"),
    name = NULL
  ) +
  scale_x_discrete(labels = gene_labels) +
  labs(
    x = NULL,
    y = "P-carrying MAGs ordered by P_score"
  ) +
  coord_cartesian(clip = "off") +
  theme_s11(base_size = 8.7) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(0.45, "lines")
  )
# ============================================================
# Combine
# ============================================================

# replace combined layout width

fig_s9 <- pA | pB

fig_s9 <- fig_s9 +
  plot_layout(widths = c(1.2, 1.0)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 9, face = "bold"),
    plot.tag.position = c(0, 1)
  )

print(fig_s9)

# ------------------------------------------------------------
# Export combined figure
# ------------------------------------------------------------
# replace combined export size only

ggsave(
  filename = file.path(out_dir, "Figure_S9_complete_genome_resolved_P_repertoire.svg"),
  plot = fig_s9,
  width = 192, height = 110, units = "mm",
  device = svglite::svglite,
  bg = "white"
)

ggsave(
  filename = file.path(out_dir, "Figure_S9_complete_genome_resolved_P_repertoire.png"),
  plot = fig_s9,
  width = 192, height = 110, units = "mm",
  dpi = 600,
  bg = "white"
)

message("Figure S19final version exported to: ", out_dir)
