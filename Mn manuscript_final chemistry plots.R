# Final Chemistry Panel Figure
# By Delaney Kirr, 10/23/2026, with help from generative AI

#-----------------------------------------------------------------------------#

# 1. Load Packages and Datasets

# Packages
library(dplyr)
library(ggplot2)
library(ggpubr)

# Load Data
mn2 <- read.csv("Mn_Formaldoxime_Oct26.csv")
fe2 <- read.csv("Fe_Ferrozine_Oct26.csv")

#-----------------------------------------------------------------------------#

# 2. Calculate mean and standard deviation across replicates (A, B, C)
mn2_averaged <- mn2 %>%
  group_by(redox, water_type, amendment, day) %>%
  summarise(
    mean_Mn_mg_L = mean(Mn_mg_L, na.rm = TRUE),
    sd_Mn_mg_L   = sd(Mn_mg_L, na.rm = TRUE),
    n_reps       = sum(!is.na(Mn_mg_L)),
    .groups = "drop"
  )

fe2_averaged <- fe2 %>%
  group_by(redox, water_type, amendment, day) %>%
  summarise(
    mean_Fe_mg_L = mean(Fe_mg_L, na.rm = TRUE),
    sd_Fe_mg_L   = sd(Fe_mg_L, na.rm = TRUE),
    n_reps       = sum(!is.na(Fe_mg_L)),
    .groups = "drop"
  )

#-----------------------------------------------------------------------------#

# 3. Aerobic Mn Plot

# Filter aerobic samples
aerobic_mn_data <- mn2_averaged %>%
  filter(redox == "Aerobic")

# Add a point at Day -1, 11 mg/L to represent initial concentration of Mn (II) added to aerobic treatments
day_minus_1_rows <- aerobic_mn_data %>%
  mutate(
    day = -1,
    mean_Mn_mg_L = 11,
    sd_Mn_mg_L = 0  # Sets error bar window to zero for the initial point
  )

# Append these baseline rows to the original dataset
aerobic_mn_data <- bind_rows(day_minus_1_rows, aerobic_mn_data)

aerobic_mn_plot <- ggplot(aerobic_mn_data, aes(x = day, y = mean_Mn_mg_L, color = amendment, linetype = water_type, shape = amendment)) +
  
# Plot lines for live treatments (Killed Controls excluded)
  geom_line(
    data = subset(aerobic_mn_data , amendment != "Killed Control"),
    aes(group = interaction(water_type, amendment)),
    linewidth = 1
  ) +
  
# Plot lines for Killed Controls (sampled less frequently)
  geom_line(
    data = subset(aerobic_mn_data , amendment == "Killed Control" & !is.na(mean_Mn_mg_L)),
    aes(group = interaction(water_type, amendment)),
    linewidth = 1
  ) +
  
# Add point layer on top of lines
  geom_point(size = 3) +
  
# Error bars
  geom_errorbar(
    aes(ymin = mean_Mn_mg_L - sd_Mn_mg_L, ymax = mean_Mn_mg_L + sd_Mn_mg_L),
    width = 0.15,
    linewidth = 0.7
  ) +
  
# Line marking carbon-amended treatment transition to anaerobic conditions
  geom_vline(
    aes(xintercept = 1, alpha = "Carbon-amended Treatments DO < 1 mg/L"), 
    colour = "black", 
    linetype = "twodash"
  ) +
  
# Color mapping scale
  scale_color_manual(
    values = c("Carbon-amended" = "#7B2D8B", "Unamended" = "#0072B2", "Killed Control" = "#E69F00")
  ) +
  
# Linetype scales 
  scale_linetype_manual(
    values = c("Groundwater" = "dashed", "SWIFT Water" = "solid"), 
  ) +
  
# Alpha scale controls opacity and ensures vline is in its own legend entry
  scale_alpha_manual(
    values = c("Amended Treatments DO < 1 mg/L" = 1)
  ) +
  
# Shapes
  scale_shape_manual(
    values = c("Carbon-amended" = 16, "Unamended" = 17, "Killed Control" = 18)
  ) +
  
# Axis adjustments
  scale_y_continuous(limits = c(0, 11.5)) +
  
# Text Labels 
  labs(
    x = "Day",
    y = "Mn (II) (mg/L)",
    color = "Treatment",
    shape = "Treatment",
    linetype = "Water Type"
  ) +
  
# Theme definitions 
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.box = "vertical",         
    legend.box.just = "center", 
    legend.position = "bottom"        
  )

aerobic_mn_plot

#-----------------------------------------------------------------------------#

# 4. Anaerobic Mn Plot
# Using secondary axis for the carbon-amended treatments

# Filter aerobic samples
anaerobic_mn_data <- mn2_averaged %>%
  filter(redox == "Anaerobic")

# Define scaling coefficient
coeff <- 20

# Create a adjusted Y variable for secondary axis data
# "Carbon-amended" stays raw, while "Unamended" and "Killed Control" are scaled up
anaerobic_mn_data <- anaerobic_mn_data %>%
  mutate(
    # Primary Y value mapping
    Mn_display = if_else(amendment == "Carbon-amended", mean_Mn_mg_L, mean_Mn_mg_L * coeff),
    
    ymin_display = case_when(
      is.na(sd_Mn_mg_L) ~ Mn_display, # No error bar if stdev is missing
      amendment == "Carbon-amended" ~ mean_Mn_mg_L - sd_Mn_mg_L,
      TRUE ~ (mean_Mn_mg_L - sd_Mn_mg_L) * coeff
    ),
    ymax_display = case_when(
      is.na(sd_Mn_mg_L) ~ Mn_display, # No error bar if stdev is missing
      amendment == "Carbon-amended" ~ mean_Mn_mg_L + sd_Mn_mg_L,
      TRUE ~ (mean_Mn_mg_L + sd_Mn_mg_L) * coeff
    )
  )

# Build the plot using the transformed y values
anaerobic_mn_plot <- ggplot(anaerobic_mn_data, aes(x = day, y = Mn_display, color = amendment, shape = amendment, linetype = water_type)) +
  
  # Plot lines for live treatments (Killed Controls excluded)
  geom_line(
    data = subset(anaerobic_mn_data , amendment != "Killed Control"),
    aes(group = interaction(water_type, amendment)),
    linewidth = 1
  ) +
  
  # Plot lines for Killed Controls (sampled less frequently)
  geom_line(
    data = subset(anaerobic_mn_data , amendment == "Killed Control" & !is.na(mean_Mn_mg_L)),
    aes(group = interaction(water_type, amendment)),
    linewidth = 1
  ) +
  
  geom_point(size = 3) +
  
  # Error bars scaled to match the corresponding primary or secondary alignment
  geom_errorbar(
    aes(
      ymin = ymin_display,
      ymax = ymax_display,
      group = interaction(amendment, water_type) 
    ),
    width = 0.15,
    linewidth = 0.7
  ) +
  
  # Color mapping scale
  scale_color_manual(
    values = c("Carbon-amended" = "#7B2D8B", "Unamended" = "#0072B2", "Killed Control" = "#E69F00")
  ) +
  
  # Linetype scales 
  scale_linetype_manual(
    values = c("Groundwater" = "dashed", "SWIFT Water" = "solid"), 
  ) +
  
  # Shapes
  scale_shape_manual(values = c("Carbon-amended" = 16, "Unamended" = 17, "Killed Control" = 18)) +
  
  # Dual-Axis Scaling Configuration
  scale_y_continuous(
    name = "Mn (II) (mg/L) - Carbon-amended", limits = c(-1, 85),
    sec.axis = sec_axis(~ . / coeff, name = "Mn (II) (mg/L) - Unamended & Killed")
  ) +
  
  labs(
    x = "Day",
    color = "Treatment",
    linetype = "Experiment",
    shape = "Treatment"
  ) +
  
  # Theme definitions 
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.box = "vertical",         
    legend.box.just = "center", 
    legend.position = "bottom" 
  ) 

anaerobic_mn_plot

#-----------------------------------------------------------------------------#

# 5. Aerobic Fe Plot

# Filter aerobic samples
aerobic_fe_data <- fe2_averaged %>%
  filter(redox == "Aerobic")

aerobic_fe_plot <- ggplot(aerobic_fe_data, aes(x = day, y = mean_Fe_mg_L, color = amendment, shape = amendment, linetype = water_type)) +
  
  # Plot lines for live treatments
  geom_line(
    data = subset(aerobic_fe_data , amendment != "Killed Control"),
    aes(group = interaction(water_type, amendment)),
    linewidth = 1
  ) +
  
  # Plot lines for Killed Controls (sampled less frequently)
  geom_line(
    data = subset(aerobic_fe_data , amendment == "Killed Control" & !is.na(mean_Fe_mg_L)),
    aes(group = interaction(water_type, amendment)),
    linewidth = 1
  ) +
  
  geom_point(size = 3) +
  
  geom_errorbar(
    aes(
      ymin = mean_Fe_mg_L - sd_Fe_mg_L,
      ymax = mean_Fe_mg_L + sd_Fe_mg_L
    ),
    width = 0.15,
    linewidth = 0.7
  ) +
  
  # Line marking carbon-amended treatment transition to anaerobic conditions
  geom_vline(
    aes(xintercept = 1, alpha = "Carbon-amended Treatments DO < 1 mg/L"), 
    colour = "black", 
    linetype = "twodash"
  ) +
  
  # Color mapping scale
  scale_color_manual(
    values = c("Carbon-amended" = "#7B2D8B", "Unamended" = "#0072B2", "Killed Control" = "#E69F00")
  ) +
  
  # Linetype scales 
  scale_linetype_manual(
    values = c("Groundwater" = "dashed", "SWIFT Water" = "solid"), 
  ) +
  
  # Alpha scale controls opacity and ensures vline is in its own legend entry
  scale_alpha_manual(
    values = c("Amended Treatments DO < 1 mg/L" = 1)
  ) +
  
  # Shapes
  scale_shape_manual(
    values = c("Carbon-amended" = 16, "Unamended" = 17, "Killed Control" = 18)
  ) +
  
  theme_minimal() +
  scale_y_continuous(limits = c(0, 6.5)) +
  scale_x_continuous(limits = c(-1, 21)) +
  
  labs(
    x = "Day",
    y = "Fe (II) (mg/L)",
    color = "Treatment",
    shape = "Treatment",
    linetype = "Experiment",
    alpha = NULL
  ) +
  
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

aerobic_fe_plot

#-----------------------------------------------------------------------------#

# 6. Anaerobic Fe Plot

# Filter aerobic samples
anaerobic_fe_data <- fe2_averaged %>%
  filter(redox == "Anaerobic")

anaerobic_fe_plot <- ggplot(anaerobic_fe_data, aes(x = day, y = mean_Fe_mg_L, color = amendment, shape = amendment, linetype = water_type)) +
  
  # Plot lines for live treatments
  geom_line(
    data = subset(anaerobic_fe_data , amendment != "Killed Control"),
    aes(group = interaction(water_type, amendment)),
    linewidth = 1
  ) +
  
  # Plot lines for Killed Controls (sampled less frequently)
  geom_line(
    data = subset(anaerobic_fe_data , amendment == "Killed Control" & !is.na(mean_Fe_mg_L)),
    aes(group = interaction(water_type, amendment)),
    linewidth = 1
  ) +
  
  geom_point(size = 3) +
  
  geom_errorbar(
    aes(
      ymin = mean_Fe_mg_L - sd_Fe_mg_L,
      ymax = mean_Fe_mg_L + sd_Fe_mg_L
    ),
    width = 0.15,
    linewidth = 0.7
  ) +
  
  # Color mapping scale
  scale_color_manual(
    values = c("Carbon-amended" = "#7B2D8B", "Unamended" = "#0072B2", "Killed Control" = "#E69F00")
  ) +
  
  # Linetype scales 
  scale_linetype_manual(
    values = c("Groundwater" = "dashed", "SWIFT Water" = "solid"), 
  ) +
  
  # Shapes
  scale_shape_manual(
    values = c("Carbon-amended" = 16, "Unamended" = 17, "Killed Control" = 18)
  ) +
  
  theme_minimal() +
  
  labs(
    x = "Day",
    y = "Fe (II) (mg/L)",
    color = "Treatment",
    shape = "Treatment",
    linetype = "Experiment",
    alpha = NULL
  ) +
  
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

anaerobic_fe_plot

#-----------------------------------------------------------------------------#

# 7. Create Combined Figure with all Four Plots

final_chemistry_plots <- ggarrange(anaerobic_mn_plot, aerobic_mn_plot, anaerobic_fe_plot, aerobic_fe_plot, ncol = 2, nrow = 2, align = "v", labels = c("A", "B", "C", "D"), common.legend = TRUE, legend = "bottom")

final_chemistry_plots

ggsave(filename = "final_chemistry_plots.png", plot = final_chemistry_plots, width = 10, height = 8.5, units = "in", dpi = 600)