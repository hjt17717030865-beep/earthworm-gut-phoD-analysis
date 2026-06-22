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

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(svglite)
  library(readxl)
})

# =========================
# 0) Paths
# =========================
out_dir <- OUTPUT_DIR
# =========================
# =========================
# 1) Read data
# =========================
plot_df <- readxl::read_excel(
  DATA_XLSX,
  sheet = "gut16S_PC1_loadings"
) %>%
  mutate(
    Loading = as.numeric(Loading)
  )
# =========================
# 2) Clean labels and make unique plotting labels
# =========================
plot_df2 <- plot_df %>%
  mutate(
    Direction = factor(Direction, levels = c("Negative", "Positive")),
    Taxon_label = trimws(as.character(Taxon_label)),
    Phylum = as.character(Phylum)
  ) %>%
  group_by(Direction, Taxon_label) %>%
  mutate(label_id = row_number()) %>%
  ungroup() %>%
  mutate(
    Taxon_plot = ifelse(
      ave(Taxon_label, Direction, Taxon_label, FUN = length) > 1,
      paste0(Taxon_label, "_", label_id),
      Taxon_label
    )
  )

plot_df2 <- plot_df2 %>%
  arrange(Direction, Loading) %>%
  mutate(Taxon_plot = make.unique(as.character(Taxon_plot))) %>%
  mutate(Taxon_plot = factor(Taxon_plot, levels = Taxon_plot))

# =========================
# 3) Colors
# =========================
phylum_cols <- c(
  "Actinobacteriota" = "#D39A92",
  "Chloroflexi"      = "#B8A873",
  "Crenarchaeota"    = "#8FAE7E",
  "Cyanobacteria"    = "#7FAFB4",
  "Desulfobacterota" = "#7FA6B8",
  "Nitrospirota"     = "#A79ACB",
  "Proteobacteria"   = "#C78FB8"
)

miss_cols <- setdiff(unique(plot_df2$Phylum), names(phylum_cols))
if (length(miss_cols) > 0) {
  phylum_cols <- c(phylum_cols, setNames(rep("grey70", length(miss_cols)), miss_cols))
}

# =========================
# 4) Plot
# =========================
p_s5 <- ggplot(plot_df2, aes(x = Taxon_plot, y = Loading, fill = Phylum)) +
  geom_col(width = 0.72, colour = NA) +
  coord_flip() +
  facet_wrap(~Direction, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = phylum_cols, drop = FALSE) +
  scale_y_continuous(
    limits = c(-0.047, 0.043),
    breaks = c(-0.04, -0.02, 0, 0.02, 0.04),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  labs(
    x = NULL,
    y = "PC1 loading",
    title = "Taxonomic contributors to gut16S_PC1"
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.title = element_text(size = 9, colour = "grey20"),
    axis.text.x = element_text(size = 8.5, colour = "grey25"),
    axis.text.y = element_text(size = 8.5, colour = "grey25"),
    axis.line = element_line(linewidth = 0.30, colour = "grey25"),
    axis.ticks = element_line(linewidth = 0.30, colour = "grey25"),
    axis.ticks.length = unit(2, "pt"),
    strip.background = element_blank(),
    strip.text = element_text(size = 9, face = "bold", colour = "grey20"),
    panel.spacing = unit(0.8, "lines"),
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8.5),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5, colour = "grey20"),
    plot.margin = margin(4, 6, 4, 4, unit = "mm")
  ) +
  guides(fill = guide_legend(override.aes = list(size = 3)))

print(p_s5)

# =========================
# 5) Export
# =========================
ggsave(
  filename = file.path(out_dir, "FigS5_gut16S_PC1_taxonomic_contributors.svg"),
  plot = p_s5,
  width = 170, height = 95, units = "mm",
  device = svglite::svglite
)

ggsave(
  filename = file.path(out_dir, "FigS5_gut16S_PC1_taxonomic_contributors.png"),
  plot = p_s5,
  width = 170, height = 95, units = "mm",
  dpi = 600
)
