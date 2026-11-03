# ============================================================
# Beta diversity plots and PERMANOVA
# ============================================================

# Written by Delaney Kirr, with assistance from generative AI, 10/21/2026

# ---- Load libraries ----
library(vegan)
library(ggplot2)
library(tidyverse)
library(qiime2R)

# ---- Load microbial data ----
meta <- read.csv("gw-sw-metadata.csv")
wp <- read_qza("level6-table.qza")

pdata <- wp$data %>% as.data.frame()

# Make relative abundance to avoid community size bias
for (i in 1:120) {
  pdata[,i] <- pdata[,i] / sum(pdata[,i])
}

# Transpose so rows = samples, columns = taxa
otu_wide <- t(pdata) %>% as.data.frame()

# Average triplicates using meta
otu_avg <- otu_wide %>%
  rownames_to_column(var = "sample.id") %>%
  left_join(meta %>% select(sample.id, redox, water_type, day, date), by = "sample.id") %>%
  group_by(redox, water_type, day, date) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  unite("bottle_ID", redox, water_type, day, date, sep = "_", remove = FALSE) %>%  # create a unique group ID
  column_to_rownames(var = "bottle_ID")

# Averaged metadata (one row per group)
metadata_avg <- otu_avg %>%
  select(redox, water_type, day, date)  # already one row per group after averaging

# Drop metadata columns before running NMDS (needs numeric taxa only)
otu_matrix <- otu_avg %>% select(-redox, -water_type, -day, -date)

# Run NMDS
set.seed(123)
nmds <- metaMDS(
  otu_matrix,
  distance = "bray",
  k = 2,
  trymax = 100
)
nmds$stress

# Extract coords and join metadata
nmds_coords <- as.data.frame(scores(nmds, display = "sites")) %>%
  rownames_to_column(var = "bottle_ID") %>%
  left_join(metadata_avg %>% rownames_to_column(var = "bottle_ID"), by = "bottle_ID")

ellipse_data <- nmds_coords #%>%
  #filter(!(water_type %in% c("GW_Media", "SW_Media") & day == 0)) # Uncomment this to exclude day 0 points from the ellipses for the amended treatments

# Static plot
ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2,
                        color = water_type,
                        shape = redox,
                        linetype = redox,
                        size = day)) +
  geom_point(aes(size = day), alpha = 0.8) +
  scale_size(range = c(3, 7), guide = "none") +
  #geom_text(aes(label = day), size = 2.5, vjust = 0.4, color = "black") + # Uncomment to add numbers for days in the middle of each point

  stat_ellipse(data = ellipse_data,
               aes(group = interaction(water_type, redox),
                   color = water_type,
                   linetype = redox)) +
  scale_color_manual(
    values = c("GW" = "#56B4E9",
               "GW_Media" = "#005b8e",
               "SW" = "#E0B0CA",
               "SW_Media" = "#A35A84"),
    labels = c("GW" = "Unamended Groundwater",
               "GW_Media" = "Carbon-amended Groundwater",
               "SW" = "Unamended SWIFT Water",
               "SW_Media" = "Carbon-amended SWIFT Water")
  ) +
  scale_shape_manual(
    values = c("oxidation" = 16, "reduction" = 17),
    labels = c("oxidation" = "Aerobic", "reduction" = "Anaerobic")
  ) +
  scale_linetype_manual(
    values = c("oxidation" = "solid", "reduction" = "dashed"),
    labels = c("oxidation" = "Aerobic", "reduction" = "Anaerobic")
  ) +
  # Merge shape and linetype into one legend under "Redox Condition"
  guides(
    color = guide_legend(title = "Water Type"),
    shape = guide_legend(title = "Redox Condition",
                         override.aes = list(
                             color = c("black", "black"),  # line and point outline
                             alpha = c(0.4, 0.4)           # dim the points slightly to look "grey"
                           )
                         ),
    linetype = guide_legend(title = "Redox Condition",
                            override.aes = list(
                              shape = c(16, 17),
                              color = c("black", "black"),
                              fill = c("black", "black")
                            ))
  ) +
  annotate("text", x = Inf, y = Inf,
           label = paste("Stress:", round(nmds$stress, 3)),
           hjust = 1.1, vjust = 1.5, size = 3.5) +
  theme_bw() +
  labs(title = "NMDS of Microbial Community Composition (Bray-Curtis)")

# ---- PERMANOVA ----

# Make sure predictors are categorical
metadata_avg <- metadata_avg %>%
  mutate(
    redox = factor(redox),
    water_type = factor(water_type),
    day = as.numeric(as.character(day))
  )

# Ensure metadata rows match community-matrix rows
metadata_avg <- metadata_avg[rownames(otu_matrix), , drop = FALSE]

# ---- Carbon-amended PERMANOVA ----
amended_meta <- metadata_avg %>%
  filter(water_type %in% c("GW_Media", "SW_Media")) %>%
  droplevels()

amended_otu <- otu_matrix[rownames(amended_meta), , drop = FALSE]

# Overall PERMANOVA
adonis2(
  amended_otu ~ day * redox * water_type,
  data = amended_meta,
  method = "bray",
  permutations = 999
)

# Sequential type 1 sum of squares PERMANOVA
amended_PERMANOVA <- adonis2(
  amended_otu ~ day * redox * water_type,
  data = amended_meta,
  method = "bray",
  permutations = 999,
  by = "term" #assesses each variable sequentially
)

print(amended_PERMANOVA, digits = 6)

# ---- Unamended PERMANOVA ----
unamended_meta <- metadata_avg %>%
  filter(water_type %in% c("GW", "SW")) %>%
  droplevels()

unamended_otu <- otu_matrix[rownames(unamended_meta), , drop = FALSE]

# Overall PERMANOVA
adonis2(
  unamended_otu ~ day * redox * water_type,
  data = amended_meta,
  method = "bray",
  permutations = 999
)

# Sequential type 1 sum of squares PERMANOVA
unamended_PERMANOVA <- adonis2(
  unamended_otu ~ day * redox * water_type,
  data = unamended_meta,
  method = "bray",
  permutations = 999,
  by = "term" #assesses each variable sequentially
)

print(unamended_PERMANOVA, digits = 6)
