# Taxonomic Bar Plot of Families in all Microcosms

# Authored by Delaney Kirr, with assistance from Generative AI

# Last edited: 11/6/2026

# ---- Libraries ----
library(qiime2R)
library(tidyverse)
library(forcats)
library(pals)
library(ggplot2)
library(ggpubr)

# ---- Load data ----
meta <- read.csv("gw-sw-metadata.csv")
wp   <- read_qza("level6-table.qza")

pdata <- wp$data %>% as.data.frame()

# Relative abundance
for (i in 1:ncol(pdata)) {
  pdata[,i] <- pdata[,i] / sum(pdata[,i])
}

# ---- Pivot to long format ----
long_wp <- pdata %>%
  rownames_to_column(var = "Taxonomy") %>%
  pivot_longer(cols = -Taxonomy, names_to = "samples", values_to = "abundance")

# ---- Join metadata and average replicates ----
avg_wp <- long_wp %>%
  left_join(meta, by = c("samples" = "sample.id")) %>%
  group_by(Taxonomy, date, redox, water_type, day) %>%
  summarise(abundance = mean(abundance), .groups = "drop")

# ---- Extract Family, clean up f__ prefix ----
avg_wp <- avg_wp %>%
  mutate(
    Family = case_when(
      grepl("f__", Taxonomy) ~ sub(".*f__(.*?)(;.*|$)", "\\1", Taxonomy),
      grepl("o__", Taxonomy) ~ paste0("o__", sub(".*o__(.*?)(;.*|$)", "\\1", Taxonomy)),
      TRUE ~ "Unassigned"
    ),
    # Clean up blank family names
    Family = case_when(
      is.na(Family)              ~ "Unassigned",
      Family == ""               ~ "Unassigned",
      grepl("^\\s*$", Family)    ~ "Unassigned",
      TRUE                       ~ Family
    )
  )

# ---- Treatment labels ----
treatment_labels <- c(
  "oxidation_GW"       = "Aerobic Unamended Groundwater",
  "oxidation_GW_Media" = "Aerobic Carbon-amended Groundwater",
  "oxidation_SW"       = "Aerobic Unamended SWIFT Water",
  "oxidation_SW_Media" = "Aerobic Carbon-amended SWIFT Water",
  "reduction_GW"       = "Anaerobic Unamended Groundwater",
  "reduction_GW_Media" = "Anaerobic Carbon-amended Groundwater",
  "reduction_SW"       = "Anaerobic Unamended SWIFT Water",
  "reduction_SW_Media" = "Anaerobic Carbon-amended SWIFT Water"
)

avg_wp <- avg_wp %>%
  mutate(treatment_key = paste(redox, water_type, sep = "_"),
         treatment_label = treatment_labels[treatment_key])

# ---- Top 20 families globally ----
top_families <- avg_wp %>%
  group_by(Family) %>%
  summarise(total = sum(abundance, na.rm = TRUE)) %>%
  filter(Family != "Unassigned", Family != "") %>%
  arrange(desc(total)) %>%
  slice_head(n = 20) %>%
  pull(Family)

# ---- Collapse to top 20 + Other ----
plot_wp <- avg_wp %>%
  mutate(
    Family_Grouped = case_when(
      Family %in% top_families ~ Family,
      TRUE                     ~ "Other"
    ),
    Family_Grouped = fct_relevel(Family_Grouped, "Other", after = Inf)
  ) %>%
  group_by(treatment_label, redox, water_type, day, Family_Grouped) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  mutate(day = factor(day, levels = sort(unique(as.numeric(as.character(day))))))

# ---- Color palette ----
family_levels <- c(top_families, "Other", "Unassigned")
n_col         <- length(top_families)

my_colors <- setNames(
  as.character(pals::kelly(n_col)),
  top_families
)
my_colors["Other"]      <- "#D3D3D3"
my_colors["Unassigned"] <- "#D3D3D3"

# ---- Plot function ----
make_bar_plot <- function(redox_val, water_val, y_label = "", x_label = "") {
  
  label <- treatment_labels[paste(redox_val, water_val, sep = "_")]
  
  plot_wp %>%
    filter(redox == redox_val, water_type == water_val) %>%
    ggplot(aes(x = day, y = abundance, fill = Family_Grouped)) +
    geom_bar(stat = "identity", position = "stack", width = 0.85) +
    scale_fill_manual(values = my_colors, name = "Family") +
    scale_y_continuous(labels = scales::percent_format(),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = label, x = x_label, y = y_label) +
    theme_bw() +
    theme(
      axis.text.x     = element_text(angle = 0, hjust = 0.5),
      plot.title      = element_text(size = 10), 
      legend.key.size = unit(0.3, "cm"),
      legend.text     = element_text(size = 8),
      legend.title    = element_text(size = 8),
      legend.position = "bottom",
      panel.grid      = element_blank()
    ) +
    guides(fill = guide_legend(ncol = 8))
}

# ---- Generate all 8 plots ----
mnob_GW       <- make_bar_plot("oxidation", "GW",       y_label = "Relative Abundance")
mnob_GW_media <- make_bar_plot("oxidation", "GW_Media")
mnob_SW       <- make_bar_plot("oxidation", "SW")
mnob_SW_media <- make_bar_plot("oxidation", "SW_Media")
mnrb_GW       <- make_bar_plot("reduction", "GW",       y_label = "Relative Abundance")
mnrb_GW_media <- make_bar_plot("reduction", "GW_Media")
mnrb_SW       <- make_bar_plot("reduction", "SW")
mnrb_SW_media <- make_bar_plot("reduction", "SW_Media")

# ---- Combined figure ----
combo_barplot <- ggarrange(
  mnob_GW,       mnob_GW_media, mnob_SW,       mnob_SW_media,
  mnrb_GW,       mnrb_GW_media, mnrb_SW,       mnrb_SW_media,
  ncol          = 4,
  nrow          = 2,
  common.legend = TRUE,
  legend        = "bottom"
)

annotate_figure(combo_barplot,
                bottom = text_grob("Day", size = 10, vjust = -8),
                top = text_grob("Day", size = 10, vjust = 28)
)

#1300 width, 700 height
