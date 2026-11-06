# SWIFT-Mn-Manuscript
Microcosm experiments - plots and statistics
by Delaney Kirr, MS Thesis Research

"chemistry_plots.R" generates the final Mn (II) and Fe (II) plots.
"alpha_diversity_plots_stats.R" plots alpha diversity metrics (Shannon Diversity, Simpson Index, Inverse Simpson, Observed Taxa, and Chao1 Richness) and executes generalized additive models for each of these metrics.
"beta_diversity_plots_stats.R" creates a non-metric multidimentional scaling (NMDS) plot of microbial community compositions and executes permutation multivariate analysis of v ariance (PERMANOVA) analyses.
"taxonomic_barplot_family.R" makes a taxonomic bar plot of the families present in all microcosms.
"random_forest_model_plots.R" executes a random forest model (inputs = microbial community and Mn (II) data) and creates permutation importance score plots, scatter plots of relative abundance over time, and taxonomic bar plots of the taxa with the top 15 highest permutation importance scores.
