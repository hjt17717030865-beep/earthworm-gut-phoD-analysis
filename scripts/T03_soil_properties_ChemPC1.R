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
# Supplementary Table 3
# Soil physicochemical properties and ChemPC1 summary
#
# Statistics:
#   lm(response ~ Duration * Regime)
#   Type-II ANOVA, consistent with Fig.2
#   Tukey-adjusted post hoc letters for Duration × Regime groups
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(car)
  library(emmeans)
  library(multcomp)
  library(openxlsx)
})

# -----------------------------
# -----------------------------
# 1) Paths
# -----------------------------
xlsx_path <- DATA_XLSX
sheet_name <- "Master_all_data"

out_dir <- OUTPUT_DIR
if (!file.exists(xlsx_path)) {
  stop("Input file not found. Please check xlsx_path: ", xlsx_path)
}

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

out_xlsx <- file.path(out_dir, "Supplementary_Table_S3_soil_summary.xlsx")

# -----------------------------
# 2) Read data
# -----------------------------
dat_raw <- read_excel(
  xlsx_path,
  sheet = sheet_name
) %>%
  as.data.frame()

# -----------------------------
# 3) Variables used in Table S3
#    column names must match your Excel exactly
# -----------------------------
var_info <- data.frame(
  variable = c(
    "pH",
    "TOC (%)",
    "water content (%)",
    "NH4+ (mg/kg)",
    "NO3- (mg/kg)",
    "TN (mg/kg)",
    "TP (mg/kg)",
    "Olsen-P (mg/kg)",
    "ChemPC1"
  ),
  column = c(
    "pH",
    "TOC(%)",
    "water content(%)",
    "NH4+(mg/kg)",
    "NO3-(mg/kg)",
    "TN(mg/kg)",
    "TP(mg/kg)",
    "Olsen-P(mg/kg)",
    "ChemPC1"
  ),
  digits = c(2, 2, 2, 2, 2, 2, 2, 2, 2),
  stringsAsFactors = FALSE
)

need_cols <- c("PlotID", "Duration", "Regime", var_info$column)
miss_cols <- setdiff(need_cols, names(dat_raw))

if (length(miss_cols) > 0) {
  stop(
    "Could not find these variables in the input table:\n",
    paste(miss_cols, collapse = ", "),
    "\n\nAvailable columns are:\n",
    paste(names(dat_raw), collapse = ", ")
  )
}

# -----------------------------
# 4) Clean data
# -----------------------------
dat <- dat_raw %>%
  mutate(
    PlotID = trimws(as.character(PlotID)),
    Duration = trimws(as.character(Duration)),
    Duration = gsub("\\s+", "", Duration),
    Duration = gsub("years|Years|year|Year|yr|YR", "y", Duration),
    Duration = ifelse(Duration == "5", "5y", Duration),
    Duration = ifelse(Duration == "8", "8y", Duration),
    Duration = ifelse(Duration == "10", "10y", Duration),
    Duration = ifelse(Duration == "3y", "5y", Duration),
    Duration = factor(Duration, levels = c("5y", "8y", "10y")),
    Regime = trimws(as.character(Regime)),
    Regime = factor(Regime, levels = c("CK", "NPK", "NPKOM"))
  )

# Convert selected variables to numeric safely
for (cc in var_info$column) {
  dat[[cc]] <- suppressWarnings(as.numeric(gsub(",", "", as.character(dat[[cc]]))))
}

dat <- dat %>%
  filter(!is.na(Duration), !is.na(Regime)) %>%
  droplevels()

# -----------------------------
# 5) Helper functions
# -----------------------------
fmt_num <- function(x, digits = 2) {
  ifelse(
    is.na(x),
    "NA",
    formatC(x, format = "f", digits = digits)
  )
}

clean_letters <- function(x) {
  gsub("\\s+", "", x)
}

make_one_variable <- function(var_label, col_name, digits = 2) {
  
  dat_var <- dat %>%
    filter(is.finite(.data[[col_name]])) %>%
    droplevels()
  
  fml <- as.formula(paste0("`", col_name, "` ~ Duration * Regime"))
  fit <- lm(fml, data = dat_var)
  
  # Type-II ANOVA, same logic as Fig.2
  anova_tab <- as.data.frame(car::Anova(fit, type = 2))
  anova_tab$term <- rownames(anova_tab)
  rownames(anova_tab) <- NULL
  
  anova_tab <- anova_tab %>%
    mutate(variable = var_label, .before = 1)
  
  # Estimated marginal means for 9 Duration × Regime groups
  emm <- emmeans::emmeans(fit, ~ Duration * Regime)
  
  # Tukey-adjusted compact letter display
  cld_tab <- as.data.frame(multcomp::cld(
    emm,
    adjust = "tukey",
    Letters = letters,
    sort = TRUE,
    reversed = TRUE
  )) %>%
    mutate(
      Duration = as.character(Duration),
      Regime = as.character(Regime),
      letter = clean_letters(.group)
    ) %>%
    dplyr::select(Duration, Regime, letter)
  
  # Numeric summaries
  summary_tab <- dat_var %>%
    group_by(Duration, Regime) %>%
    summarise(
      n = sum(is.finite(.data[[col_name]])),
      mean = mean(.data[[col_name]], na.rm = TRUE),
      sd = sd(.data[[col_name]], na.rm = TRUE),
      sem = sd / sqrt(n),
      .groups = "drop"
    ) %>%
    mutate(
      Duration = as.character(Duration),
      Regime = as.character(Regime)
    ) %>%
    left_join(cld_tab, by = c("Duration", "Regime")) %>%
    mutate(
      variable = var_label,
      value = paste0(
        fmt_num(mean, digits),
        " ± ",
        fmt_num(sem, digits),
        " ",
        letter
      )
    ) %>%
    dplyr::select(variable, Duration, Regime, n, mean, sd, sem, letter, value)
  
  # Tukey pairwise contrasts
  tukey_tab <- as.data.frame(pairs(emm, adjust = "tukey")) %>%
    mutate(variable = var_label, .before = 1)
  
  list(
    summary = summary_tab,
    anova = anova_tab,
    tukey = tukey_tab
  )
}

# -----------------------------
# 6) Run all variables
# -----------------------------
res_list <- lapply(seq_len(nrow(var_info)), function(i) {
  make_one_variable(
    var_label = var_info$variable[i],
    col_name  = var_info$column[i],
    digits    = var_info$digits[i]
  )
})

summary_long <- bind_rows(lapply(res_list, `[[`, "summary"))
anova_all    <- bind_rows(lapply(res_list, `[[`, "anova"))
tukey_all    <- bind_rows(lapply(res_list, `[[`, "tukey"))

# -----------------------------
# 7) Make submission-style wide table
# -----------------------------
group_n <- dat %>%
  group_by(Duration, Regime) %>%
  summarise(n_plot = n(), .groups = "drop") %>%
  mutate(
    Duration = as.character(Duration),
    Regime = as.character(Regime)
  )

table_s3 <- summary_long %>%
  dplyr::select(Duration, Regime, variable, value) %>%
  tidyr::pivot_wider(
    names_from = variable,
    values_from = value
  ) %>%
  dplyr::left_join(group_n, by = c("Duration", "Regime")) %>%
  dplyr::mutate(
    Duration = factor(Duration, levels = c("5y", "8y", "10y")),
    Regime = factor(Regime, levels = c("CK", "NPK", "NPKOM"))
  ) %>%
  dplyr::arrange(Duration, Regime) %>%
  dplyr::select(
    Duration,
    Regime,
    n = n_plot,
    dplyr::all_of(var_info$variable)
  )

# Numeric long table
summary_numeric <- summary_long %>%
  mutate(
    Duration = factor(Duration, levels = c("5y", "8y", "10y")),
    Regime = factor(Regime, levels = c("CK", "NPK", "NPKOM"))
  ) %>%
  arrange(variable, Duration, Regime)

# -----------------------------
# 8) Write Excel output
# -----------------------------
wb <- createWorkbook()

addWorksheet(wb, "Table_S3")
writeData(wb, "Table_S3", table_s3)

addWorksheet(wb, "Table_S3_numeric")
writeData(wb, "Table_S3_numeric", summary_numeric)

addWorksheet(wb, "TypeII_ANOVA")
writeData(wb, "TypeII_ANOVA", anova_all)

addWorksheet(wb, "Tukey_pairwise")
writeData(wb, "Tukey_pairwise", tukey_all)

addWorksheet(wb, "README")
readme <- data.frame(
  Item = c(
    "Purpose",
    "Model",
    "ANOVA",
    "Post hoc letters",
    "Values",
    "Letter interpretation"
  ),
  Description = c(
    "Supplementary Table 3 summarises soil physicochemical properties and ChemPC1 across Duration × Regime groups.",
    "Each variable was analysed using lm(response ~ Duration * Regime).",
    "Type-II ANOVA was performed using car::Anova(type = 2), consistent with Fig.2.",
    "Compact letters were generated from emmeans for the nine Duration × Regime combinations with Tukey adjustment.",
    "Values in Table_S3 are shown as mean ± s.e.m. followed by compact-letter display.",
    "Groups not sharing a letter differ significantly based on Tukey-adjusted pairwise comparisons."
  )
)
writeData(wb, "README", readme)

# Basic formatting
header_style <- createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "Bottom"
)

for (sh in names(wb)) {
  addStyle(
    wb,
    sheet = sh,
    style = header_style,
    rows = 1,
    cols = 1:ncol(readWorkbook(wb, sheet = sh)),
    gridExpand = TRUE
  )
  freezePane(wb, sheet = sh, firstRow = TRUE)
  setColWidths(wb, sheet = sh, cols = 1:50, widths = "auto")
}

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

cat("Saved Supplementary Table S3 to:\n", out_xlsx, "\n")
