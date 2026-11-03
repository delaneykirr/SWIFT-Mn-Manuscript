# ============================================================
# Alpha diversity plots and Generalized Additive Models
# ============================================================

# Written by Delaney Kirr and Allison Mason, with assistance from generative AI, 10/21/2026

# ---- Load libraries ----
library(vegan)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(qiime2R)
library(ggpubr)
library(mgcv)
library(modelsummary)

# ---- Load microbial data ----
meta <- read.csv(
  "gw-sw-metadata.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "NaN")
)

wp <- read_qza("level6-table.qza")
pdata <- as.data.frame(wp$data)

# QIIME 2 tables are commonly taxa x samples. Convert to samples x taxa.
otu_wide <- as.data.frame(t(pdata), check.names = FALSE)
otu_wide[] <- lapply(otu_wide, function(x) as.numeric(as.character(x)))

# Accept the usual metadata column names: sample-id, sample.id, or sample_id.
sample_id_column <- intersect(c("sample-id", "sample.id", "sample_id"), names(meta))[1]
if (is.na(sample_id_column)) {
  stop(
    "Could not find the sample ID column. Expected one of: sample-id, sample.id, sample_id.\n",
    "Columns found: ", paste(names(meta), collapse = ", ")
  )
}

names(meta)[names(meta) == sample_id_column] <- "sample.id"

# Normalize IDs while preserving their identity. This also handles dates such as
# 12-1-2025 versus 12-01-2025, which otherwise fail an exact join.
normalize_sample_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^\\ufeff", "", x)
  x <- gsub("[[:space:]]+", "", x)

  date_parts <- stringr::str_match(
    x,
    "^(.*-)([0-9]{1,2})-([0-9]{1,2})-([0-9]{4})$"
  )

  has_date <- !is.na(date_parts[, 1])
  x[has_date] <- sprintf(
    "%s%02d-%02d-%s",
    date_parts[has_date, 2],
    as.integer(date_parts[has_date, 3]),
    as.integer(date_parts[has_date, 4]),
    date_parts[has_date, 5]
  )
  x
}

meta <- meta %>%
  mutate(sample.id = normalize_sample_id(sample.id))

otu_wide$sample.id <- normalize_sample_id(rownames(otu_wide))

# Check the join before calculating diversity so mismatches are explicit.
missing_metadata <- setdiff(otu_wide$sample.id, meta$sample.id)
if (length(missing_metadata) > 0) {
  stop(
    "These QIIME 2 samples are missing from the metadata after ID normalization:\n",
    paste(missing_metadata, collapse = ", "),
    "\n\nCheck the metadata sample-id values and spelling."
  )
}

if (anyDuplicated(meta$sample.id)) {
  duplicate_ids <- unique(meta$sample.id[duplicated(meta$sample.id)])
  stop(
    "Duplicate sample IDs found in the metadata after normalization:\n",
    paste(duplicate_ids, collapse = ", ")
  )
}

# ---- Standardize metadata labels ----
# This accepts GW, GW_Media, SW, SW_Media and equivalent amended labels.
meta <- meta %>%
  mutate(
    sample.id = as.character(sample.id),
    redox_key = str_to_lower(str_trim(as.character(redox))),
    water_type_key = str_replace_all(
      str_to_upper(str_trim(as.character(water_type))),
      "[^A-Z]",
      ""
    ),
    redox_label = case_when(
      redox_key %in% c("oxidation", "oxidized", "aerobic") ~ "Aerobic",
      redox_key %in% c("reduction", "reduced", "anaerobic") ~ "Anaerobic",
      TRUE ~ NA_character_
    ),
    treatment = case_when(
      water_type_key %in% c("GW", "GROUNDWATER", "UNAMENDEDGW") ~
        "Unamended Groundwater",
      water_type_key %in% c("GWMEDIA", "AMENDEDGW", "AMENDEDGROUNDWATER") ~
        "Carbon-amended Groundwater",
      water_type_key %in% c("SW", "SWIFTWATER", "UNAMENDEDSW") ~
        "Unamended SWIFT Water",
      water_type_key %in% c("SWMEDIA", "AMENDEDSW", "AMENDEDSWIFTWATER") ~
        "Carbon-amended SWIFT Water",
      TRUE ~ NA_character_
    ),
    day = as.numeric(day)
  )

if (anyNA(meta$redox_label)) {
  stop("Unrecognized redox values: ", paste(unique(meta$redox[is.na(meta$redox_label)]), collapse = ", "))
}

if (anyNA(meta$treatment)) {
  stop("Unrecognized water_type values: ", paste(unique(meta$water_type[is.na(meta$treatment)]), collapse = ", "))
}

# ---- Calculate sample-level alpha-diversity metrics ----
# Remove the temporary sample.id column before vegan calculations.
otu_matrix <- otu_wide %>%
  select(-sample.id) %>%
  as.matrix()

rownames(otu_matrix) <- otu_wide$sample.id

chao1_values <- estimateR(otu_matrix)["S.chao1", ]

alpha_div_raw <- tibble(
  sample.id  = rownames(otu_matrix),
  shannon    = diversity(otu_matrix, index = "shannon"),
  simpson    = diversity(otu_matrix, index = "simpson"),
  invsimpson = diversity(otu_matrix, index = "invsimpson"),
  observed   = specnumber(otu_matrix),
  chao1      = as.numeric(chao1_values)
) %>%
  left_join(meta, by = "sample.id")

if (any(is.na(alpha_div_raw$redox_label))) {
  missing_samples <- alpha_div_raw$sample.id[is.na(alpha_div_raw$redox_label)]
  stop(
    "These samples are missing from the metadata or have invalid metadata: ",
    paste(missing_samples, collapse = ", ")
  )
}

# ---- Average across biological/technical replicates ----
# SD is set to zero when a group has only one non-missing replicate,
# so error bars still render instead of becoming NA.
summarise_mean_sd <- function(x) {
  n <- sum(!is.na(x))
  c(
    mean = if (n > 0) mean(x, na.rm = TRUE) else NA_real_,
    sd   = if (n > 1) sd(x, na.rm = TRUE) else 0
  )
}

alpha_div <- alpha_div_raw %>%
  group_by(redox_label, treatment, day) %>%
  summarise(
    shannon_stats    = list(summarise_mean_sd(shannon)),
    simpson_stats    = list(summarise_mean_sd(simpson)),
    invsimpson_stats = list(summarise_mean_sd(invsimpson)),
    observed_stats   = list(summarise_mean_sd(observed)),
    chao1_stats      = list(summarise_mean_sd(chao1)),
    .groups = "drop"
  ) %>%
  mutate(
    shannon_mean    = vapply(shannon_stats, `[[`, numeric(1), "mean"),
    shannon_sd      = vapply(shannon_stats, `[[`, numeric(1), "sd"),
    simpson_mean    = vapply(simpson_stats, `[[`, numeric(1), "mean"),
    simpson_sd      = vapply(simpson_stats, `[[`, numeric(1), "sd"),
    invsimpson_mean = vapply(invsimpson_stats, `[[`, numeric(1), "mean"),
    invsimpson_sd   = vapply(invsimpson_stats, `[[`, numeric(1), "sd"),
    observed_mean   = vapply(observed_stats, `[[`, numeric(1), "mean"),
    observed_sd     = vapply(observed_stats, `[[`, numeric(1), "sd"),
    chao1_mean      = vapply(chao1_stats, `[[`, numeric(1), "mean"),
    chao1_sd        = vapply(chao1_stats, `[[`, numeric(1), "sd"),
    treatment = factor(
      treatment,
      levels = c(
        "Unamended Groundwater",
        "Carbon-amended Groundwater",
        "Unamended SWIFT Water",
        "Carbon-amended SWIFT Water"
      )
    ),
    redox_label = factor(redox_label, levels = c("Aerobic", "Anaerobic")),
    redox_treatment = paste(redox_label, treatment, sep = " — ")
  ) %>%
  select(-ends_with("_stats"))

# ---- Consistent colors and shapes ----
water_colors <- c(
  "Unamended Groundwater"     = "#56B4E9",
  "Carbon-amended Groundwater" = "#005B8E",
  "Unamended SWIFT Water"      = "#E0B0CA",
  "Carbon-amended SWIFT Water" = "#A35A84"
)

redox_treatment_levels <- c(
  "Aerobic — Unamended Groundwater",
  "Aerobic — Carbon-amended Groundwater",
  "Aerobic — Unamended SWIFT Water",
  "Aerobic — Carbon-amended SWIFT Water",
  "Anaerobic — Unamended Groundwater",
  "Anaerobic — Carbon-amended Groundwater",
  "Anaerobic — Unamended SWIFT Water",
  "Anaerobic — Carbon-amended SWIFT Water"
)

# Shapes distinguish treatments while color also identifies treatment.
combined_shapes <- c(
  "Aerobic — Unamended Groundwater"      = 16,
  "Aerobic — Carbon-amended Groundwater" = 17,
  "Aerobic — Unamended SWIFT Water"      = 15,
  "Aerobic — Carbon-amended SWIFT Water" = 18,
  "Anaerobic — Unamended Groundwater"      = 16,
  "Anaerobic — Carbon-amended Groundwater" = 17,
  "Anaerobic — Unamended SWIFT Water"      = 15,
  "Anaerobic — Carbon-amended SWIFT Water" = 18
)

# ---- Generic plotting function ----
make_alpha_plot <- function(data, y_mean, y_sd, y_label, show_x_label = FALSE) {

  plot_data <- data %>%
    mutate(
      y_value = .data[[y_mean]],
      y_sd_value = replace_na(.data[[y_sd]], 0),
      ymin = y_value - y_sd_value,
      ymax = y_value + y_sd_value,
      redox_treatment = factor(redox_treatment, levels = redox_treatment_levels)
    )

  ggplot(
    plot_data,
    aes(
      x = day,
      y = y_value,
      group = redox_treatment,
      color = treatment,
      shape = redox_treatment
    )
  ) +
    geom_line(linewidth = 0.7, na.rm = TRUE) +
    geom_errorbar(
      aes(ymin = ymin, ymax = ymax),
      width = 0.15,
      linewidth = 0.4,
      na.rm = TRUE
    ) +
    geom_point(size = 3, na.rm = TRUE) +
    facet_wrap(~ redox_label) +
    scale_color_manual(
      values = water_colors,
      limits = names(water_colors),
      drop = FALSE,
      name = "Treatment"
    ) +
    scale_shape_manual(
      values = combined_shapes,
      limits = redox_treatment_levels,
      drop = FALSE,
      name = "Redox × treatment"
    ) +
    scale_x_continuous(breaks = sort(unique(plot_data$day))) +
    labs(
      x = if (show_x_label) "Day" else NULL,
      y = y_label
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8))
    )
}

# ---- Build all five plots ----
p_shannon <- make_alpha_plot(
  alpha_div, "shannon_mean", "shannon_sd", "Shannon Diversity"
)

p_simpson <- make_alpha_plot(
  alpha_div, "simpson_mean", "simpson_sd", "Simpson Index"
)

p_invsimpson <- make_alpha_plot(
  alpha_div, "invsimpson_mean", "invsimpson_sd", "Inverse Simpson"
)

p_observed <- make_alpha_plot(
  alpha_div, "observed_mean", "observed_sd", "Observed Taxa"
)

p_chao1 <- make_alpha_plot(
  alpha_div, "chao1_mean", "chao1_sd", "Chao1 Richness", show_x_label = TRUE
)

# ---- Combine into one figure ----
alpha_figure <- ggarrange(
  p_shannon,
  p_simpson,
  p_invsimpson,
  p_observed,
  p_chao1,
  ncol = 1,
  nrow = 5,
  common.legend = TRUE,
  legend = "right",
  labels = c("A", "B", "C", "D", "E")
)

print(alpha_figure)

#Export:
ggsave(
  "alpha_diversity_figure.png",
  alpha_figure,
  width = 12,
  height = 14,
  units = "in",
  dpi = 300
)

# ---- Generalized additive models (GAMs) ----

#Information about GAMs: https://r.qcbs.ca/workshop08/book-en/introduction-to-generalized-additive-mixed-models-gamms.html

#Parse treatment into amendment and water type
alpha_div_gamm <- alpha_div %>%
  mutate(
    amendment = case_when(
      grepl("Carbon-amended", treatment, ignore.case = TRUE) ~ "Carbon-amended",
      grepl("Unamended",      treatment, ignore.case = TRUE) ~ "Unamended"
    ),
    water_type = case_when(
      grepl("SWIFT", treatment, ignore.case = TRUE)       ~ "SWIFT Water",
      grepl("Groundwater", treatment, ignore.case = TRUE) ~ "Groundwater"
    )
  )

#Ensure categorical variables are factors
alpha_div_gamm$redox_label <- factor(alpha_div_gamm$redox_label)
alpha_div_gamm$amendment <- factor(alpha_div_gamm$amendment)
alpha_div_gamm$water_type <- factor(alpha_div_gamm$water_type)

#GAM for each alpha diversity metric

shannon_gam <- gam(shannon_mean ~ s(day, by = interaction(amendment, water_type), k = 5) + redox_label + amendment + water_type, data = alpha_div_gamm, method = "REML")
print(summary(shannon_gam), digits = 6)
gam.check(shannon_gam)
plot(shannon_gam, all.terms = TRUE, pages = 1)

simpson_gam <- gam(simpson_mean ~ s(day, by = interaction(amendment, water_type), k = 5) + redox_label + amendment + water_type, data = alpha_div_gamm, method = "REML")
print(summary(simpson_gam), digits = 6)
gam.check(simpson_gam)
plot(simpson_gam, all.terms = TRUE, pages = 1)

invsimpson_gam <- gam(invsimpson_mean ~ s(day, by = interaction(amendment, water_type), k = 5) + redox_label + amendment + water_type, data = alpha_div_gamm, method = "REML")
print(summary(invsimpson_gam), digits = 6)
gam.check(invsimpson_gam)
plot(invsimpson_gam, all.terms = TRUE, pages = 1)

observed_gam <- gam(observed_mean ~ s(day, by = interaction(amendment, water_type), k = 5) + redox_label + amendment + water_type, data = alpha_div_gamm, method = "REML")
print(summary(observed_gam), digits = 6)
gam.check(observed_gam)
plot(observed_gam, all.terms = TRUE, pages = 1)

chao1_gam <- gam(chao1_mean ~ s(day, by = interaction(amendment, water_type), k = 5) + redox_label + amendment + water_type, data = alpha_div_gamm, method = "REML")
print(summary(chao1_gam), digits = 6)
gam.check(chao1_gam)
plot(chao1_gam, all.terms = TRUE, pages = 1)