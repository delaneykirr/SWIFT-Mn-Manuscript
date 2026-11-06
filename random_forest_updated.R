#Random Forest Modeling of Microbial Communities Related to Mn (II) in Microcosm Experiments
#Authored by Rohit Shukla, 8/19/26
#Edited by Delaney Kirr, 8/20/26 and 11/6/2026

library(vegan)
library(ggplot2)
library(tidyverse)
library(qiime2R)
library(dplyr)
library(readr)
library(reticulate)

PY_PATH <- Sys.getenv("RETICULATE_PYTHON", unset = "/opt/anaconda3/bin/python3")
if (file.exists(PY_PATH)) use_python(PY_PATH, required = FALSE)

# =============================================================================
# 1. Settings
# =============================================================================

SEED <- 42L # everything stochastic derives from this
N_TRIALS <- 100L # Optuna trials per model, changeable
N_TRIALS_INNER <- 30L # Optuna trials inside each nested-CV fold (runs 16x), changeable
N_TREES <- 500L # not tuned; error plateaus well before this
TUNE_SEEDS <- 1:3 # OOB R2 is noisy at n = 16; average over 3 forests
MIN_OOB <- 3L # rows every tree must leave out (see section 5)
N_PERM_REPEATS <- 10L # shuffles per taxon for %IncMSE
SHAP_DIR_MIN <- 0.3 # weaker directions are drawn as "unclear"
OUTDIR <- "results_optuna"
FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")

for (d in c(OUTDIR, FIGDIR, TABDIR)) dir.create(d, showWarnings = FALSE)
set.seed(SEED)

# one theme for every figure; 300 dpi png plus a vector pdf
theme_pub <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, colour = "grey92"),
        strip.background = element_rect(fill = "grey95", colour = "grey70"),
        strip.text = element_text(face = "bold", size = 10),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, colour = "grey30"),
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"))

DIR_COLS <- c("higher with Mn" = "#1B7837", "lower with Mn" = "#762A83",
              "unclear" = "grey70")

save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIGDIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(FIGDIR, paste0(name, ".pdf")), p, width = w, height = h)
}

save_table <- function(df, name, digits = 3) {
  df <- df %>% mutate(across(where(is.numeric), ~ round(.x, digits)))
  write.csv(df, file.path(TABDIR, paste0(name, ".csv")), row.names = FALSE)
  df
}

# =============================================================================
# 2. Read in Data
# =============================================================================

read_qza_compat <- function(path) {
  if (requireNamespace("qiime2R", quietly = TRUE)) return(qiime2R::read_qza(path))
  py_run_string('
import zipfile, io, h5py, numpy as np
from scipy.sparse import csc_matrix
def _read_qza(path):
z = zipfile.ZipFile(path)
b = [n for n in z.namelist() if n.endswith("data/feature-table.biom")][0]
with z.open(b) as fh: buf = io.BytesIO(fh.read())
with h5py.File(buf, "r") as h:
obs = [x.decode() for x in h["observation/ids"][:]]
smp = [x.decode() for x in h["sample/ids"][:]]
m = csc_matrix((h["sample/matrix/data"][:], h["sample/matrix/indices"][:],
h["sample/matrix/indptr"][:]), shape=(len(obs), len(smp)))
return np.asarray(m.todense()).astype(float), obs, smp, b.split("/")[0]
')
  r <- py$`_read_qza`(path)
  m <- r[[1]]; dimnames(m) <- list(r[[2]], r[[3]])
  list(data = m, uuid = r[[4]])
}

# Load microbial data
meta <- read.csv("gw-sw-metadata.csv")
wp <- read_qza_compat("level6-table.qza")

pdata <- wp$data %>% as.data.frame()

# Make relative abundance
for (i in 1:120) {
  pdata[,i] <- pdata[,i] / sum(pdata[,i])
}

# Transpose so rows = samples, columns = taxa
otu_wide <- t(pdata) %>% as.data.frame()

# Average triplicates 
otu_avg <- otu_wide %>%
  rownames_to_column(var = "sample.id") %>%
  left_join(meta %>% select(sample.id, redox, water_type, day, date), by = "sample.id") %>%
  group_by(redox, water_type, day, date) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  unite("sample_key", redox, water_type, day, sep = "_", remove = FALSE) %>%
  column_to_rownames(var = "sample_key")

# Load geochemistry
mn2 <- read.csv("Mn_renamed.csv")

# =============================================================================
# 3. Prepare Data Frames
# =============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(randomForest)
library(ggplot2)

# Parse the "Sample" IDs in mn2 

parse_sample_id <- function(df) {
  df %>%
    mutate(
      redox = str_extract(Sample, "^(oxidation|reduction)"),
      day = as.numeric(str_extract(Sample, "\\d+$")),
      water_type = Sample %>%
        str_remove(paste0("^", redox, "_")) %>%
        str_remove(paste0("_", day, "$"))
    )
}

mn2 <- parse_sample_id(mn2)

# Trim all three data frames to day <= 16 (exponential phase) 
otu_avg <- otu_avg %>% filter(day <= 16)
mn2 <- mn2 %>% filter(day <= 16)

# Merge OTU table with Mn measurements 
# otu_avg already has redox / water_type / day columns, so join on those.

mn_data <- inner_join(
  otu_avg,
  mn2 %>% select(redox, water_type, day, Mn_mg_L),
  by = c("redox", "water_type", "day")
)

# Split each into oxidation / reduction subsets
mn_ox <- mn_data %>% filter(redox == "oxidation")
mn_red <- mn_data %>% filter(redox == "reduction")

# =============================================================================
# 4. Model inputs: taxa on one side, Mn on the other. Constant taxa go.
# =============================================================================

META_COLS <- c("redox", "water_type", "day", "date", "Mn_mg_L")

make_xy <- function(df, label) {
  X <- df %>% select(-any_of(META_COLS)) %>% as.matrix()
  storage.mode(X) <- "double"
  keep <- apply(X, 2, function(v) stats::var(v, na.rm = TRUE) > 0)
  # keep <- keep & (colMeans(X > 0) >= 0.25) # uncomment to add a prevalence filter
  X <- X[, keep, drop = FALSE]
  y <- df$Mn_mg_L
  cat(sprintf("%-10s n = %2d rows, %4d taxa retained of %4d, Mn %.2f - %.2f mg/L\n",
              label, nrow(X), ncol(X), length(keep),
              min(df$Mn_mg_L), max(df$Mn_mg_L)))
  list(X = X, y = y, label = label, meta = df[, intersect(META_COLS, names(df))])
}

cat("\n--- datasets ------------------------------------------------------\n")
DS <- list(oxidation = make_xy(mn_ox, "oxidation"),
           reduction = make_xy(mn_red, "reduction"))

# =============================================================================
# 5. Fitting the out-of-bag score 
# =============================================================================

fit_rf <- function(X, y, p, ntree = N_TREES, seed = SEED, importance = FALSE,
                   nPerm = 1L, keep.inbag = FALSE) {
  n <- nrow(X)
  ss <- max(2L, as.integer(round(p$sample_frac * n)))
  if (!p$replace) ss <- min(ss, n - MIN_OOB)
  set.seed(seed)
  randomForest(x = X, y = y, ntree = ntree,
               mtry = max(1L, min(ncol(X), as.integer(round(p$mtry_frac * ncol(X))))),
               nodesize = p$nodesize,
               maxnodes = if (is.null(p$maxnodes) || is.na(p$maxnodes)) NULL else p$maxnodes,
               replace = p$replace,
               sampsize = ss,
               importance = importance, nPerm = nPerm, keep.inbag = keep.inbag)
}

# average a few forests so Optuna optimises the model rather than a lucky seed
oob_r2 <- function(X, y, p, seeds = TUNE_SEEDS, ntree = N_TREES) {
  v <- mean(vapply(seeds, function(s) {
    rf <- fit_rf(X, y, p, ntree = ntree, seed = SEED + s)
    rf$rsq[length(rf$rsq)]
  }, numeric(1)))
  if (!is.finite(v)) -1 else v # score a bad trial rather than erroring out
}

oob_report <- function(X, y, p, seeds = 1:10, ntree = N_TREES) {
  r <- vapply(seeds, function(s) {
    rf <- fit_rf(X, y, p, ntree = ntree, seed = SEED + s)
    c(r2 = rf$rsq[length(rf$rsq)], rmse = sqrt(rf$mse[length(rf$mse)]))
  }, numeric(2))
  list(r2 = mean(r["r2", ]), r2_sd = sd(r["r2", ]),
       rmse = mean(r["rmse", ]), rmse_sd = sd(r["rmse", ]))
}

# Hyperparameter Tuning
optuna <- import("optuna", delay_load = TRUE)
optuna$logging$set_verbosity(optuna$logging$WARNING)

search_params <- function(X, y, n_trials = N_TRIALS, seed = SEED) {
  objective <- function(trial) {
    p <- list(
      mtry_frac = trial$suggest_float("mtry_frac", 0.02, 1.00),
      nodesize = trial$suggest_int("nodesize", 1L,
                                      max(2L, as.integer(ceiling(0.2 * nrow(X))))),
      sample_frac = trial$suggest_float("sample_frac", 0.4, 0.9),
      maxnodes = NA,
      replace = FALSE)
    oob_r2(X, y, p)
  }
  study <- optuna$create_study(
    direction = "maximize",
    sampler = optuna$samplers$TPESampler(seed = as.integer(seed)))
  study$optimize(objective, n_trials = as.integer(n_trials))
  bp <- study$best_params
  list(params = list(mtry_frac = bp$mtry_frac,
                     nodesize = as.integer(bp$nodesize),
                     sample_frac = bp$sample_frac,
                     maxnodes = NA_integer_,
                     replace = FALSE),
       value = study$best_value, study = study)
}

tune_one <- function(ds, n_trials = N_TRIALS) {
  X <- ds$X
  sp <- search_params(X, ds$y, n_trials)
  best <- sp$params
  cat(sprintf("%-10s best OOB R2 = %.3f (mtry_frac %.3f -> mtry %d of %d, nodesize %d, sample_frac %.2f)\n",
              ds$label, sp$value, best$mtry_frac,
              max(1L, as.integer(round(best$mtry_frac * ncol(X)))), ncol(X),
              best$nodesize, best$sample_frac))
  sp
}

cat("\n--- Optuna search (", N_TRIALS, " trials per model) -------------------\n", sep = "")
TUNED <- lapply(DS, tune_one)

# =============================================================================
# 6. Fit each tuned model and record its out-of-bag performance
# =============================================================================

# Shorten taxonomy strings and add rank 
short_label <- function(x) {
  rank_names <- c(d = "Domain", p = "Phylum", c = "Class", o = "Order",
                  f = "Family", g = "Genus", s = "Species", t = "Strain")
  
  parts <- str_split(str_replace_all(x, "\\.__|__", ";"), ";")
  vapply(parts, function(p) {
    n <- length(p)
    if (n < 2) return("unclassified")
    
    # split into (rank letter, taxon name) pairs by position
    letters <- p[seq(1, n - 1, by = 2)]
    names_  <- p[seq(2, n, by = 2)]
    
    valid <- nzchar(names_) & !names_ %in% c("uncultured", "NA")
    if (!any(valid)) return("unclassified")
    
    last <- max(which(valid))
    rank_letter <- letters[last]
    rank_label <- if (rank_letter %in% names(rank_names)) {
      rank_names[[rank_letter]]
    } else {
      toupper(rank_letter)  # fallback if an unexpected rank code shows up
    }
    
    paste0(names_[last], " (", rank_label, ")")
  }, character(1))
}

PERF <- list(); FITS <- list()

for (rx in names(DS)) {
  ds <- DS[[rx]]; X <- ds$X; y <- ds$y
  tuned_p <- TUNED[[rx]]$params
  
  # mean +- sd over 10 fresh forest seeds
  b <- oob_report(X, y, tuned_p)
  PERF[[rx]] <- tibble(
    redox = rx, n = nrow(X), n_taxa = ncol(X),
    oob_r2 = b$r2, oob_r2_sd = b$r2_sd, oob_rmse = b$rmse, oob_rmse_sd = b$rmse_sd,
    mtry = max(1L, as.integer(round(tuned_p$mtry_frac * ncol(X)))),
    nodesize = tuned_p$nodesize, sample_frac = tuned_p$sample_frac)
  
  cat(sprintf("%-10s OOB R2 = %.3f (+-%.3f), RMSE %.2f mg/L\n",
              rx, b$r2, b$r2_sd, b$rmse))
  
  # nPerm steadies the importance; keep.inbag is what section 8 needs
  rf <- fit_rf(X, y, tuned_p, importance = TRUE, nPerm = 5L, keep.inbag = TRUE)
  FITS[[rx]] <- rf
  
  if (all(importance(rf, type = 1, scale = FALSE)[, 1] == 0))
    warning(rx, ": every %IncMSE is exactly 0 -- too few out-of-bag rows per ",
            "tree for permutation importance. Raise MIN_OOB.", call. = FALSE)
  
}

# % IncMSE (OOB)

oob_pred_from <- function(rf, X, mask) {
  pa <- predict(rf, X, predict.all = TRUE)$individual
  rowSums(pa * mask) / rowSums(mask)
}

pct_incmse <- function(rf, X, y, n_rep = N_PERM_REPEATS, seed = SEED) {
  mask <- rf$inbag == 0
  ok <- rowSums(mask) > 0
  base <- mean((y[ok] - oob_pred_from(rf, X, mask)[ok])^2)
  
  if (!isTRUE(all.equal(base, rf$mse[length(rf$mse)], tolerance = 1e-6)))
    stop("OOB reconstruction (", signif(base, 6), ") does not match rf$mse (",
         signif(rf$mse[length(rf$mse)], 6), ") -- do not trust the %IncMSE below")
  
  used <- randomForest::varUsed(rf, count = TRUE) > 0
  pct <- numeric(ncol(X)); pct_sd <- numeric(ncol(X))
  set.seed(seed)
  for (j in which(used)) {
    v <- vapply(seq_len(n_rep), function(r) {
      Xp <- X; Xp[, j] <- Xp[sample.int(nrow(X)), j]
      mean((y[ok] - oob_pred_from(rf, Xp, mask)[ok])^2)
    }, numeric(1))
    pct[j] <- 100 * (mean(v) - base) / base
    pct_sd[j] <- 100 * stats::sd(v) / base
  }
  list(pct = pct, pct_sd = pct_sd, baseline_mse = base, used = used, n_used = sum(used))
}

# =============================================================================
# 7. SHAP - overall lower/higher Mn effects
# =============================================================================

#Chose to remove this for ease of discussion writing, but may use in future projects!

# shap_direction <- function(rf, X) {
#   if (!requireNamespace("treeshap", quietly = TRUE)) {
#     message("note: treeshap not installed - falling back to the marginal ",
#             "correlation for bar colour (install.packages('treeshap'))")
#     return(NULL)
#   }
#   Xdf <- as.data.frame(X)
#   uni <- treeshap::randomForest.unify(rf, Xdf)
#   ts <- treeshap::treeshap(uni, Xdf, verbose = FALSE)
#   sh <- as.matrix(ts$shaps)[, colnames(X), drop = FALSE] # match by name
#   
#   pred <- as.numeric(predict(rf, X))
#   gap <- max(abs((rowSums(sh) + mean(pred)) - pred))
#   if (gap > 1e-6)
#     stop("SHAP additivity violated by ", signif(gap, 3),
#          " -- the attributions do not reconstruct the predictions")
#   
#   list(direction = vapply(seq_len(ncol(sh)),
#                           function(j) suppressWarnings(
#                             cor(sh[, j], X[, j], method = "spearman")), numeric(1)),
#        additivity_gap = gap)
# }

# =============================================================================
# 8. The importance table and figure: %IncMSE bars, coloured by direction
# =============================================================================

# Version of the figure with SHAP coloring bars by direction

# cat("\n--- %IncMSE (", N_PERM_REPEATS, " permutations per taxon) ---------------\n", sep = "")
# INCMSE <- list()
# 
# for (rx in names(FITS)) {
#   rf <- FITS[[rx]]; X <- DS[[rx]]$X; y <- DS[[rx]]$y
#   pi <- pct_incmse(rf, X, y)
#   sd_ <- shap_direction(rf, X) # NULL if treeshap is missing
#   
#   tab <- tibble(
#     taxon = colnames(X),
#     label = short_label(colnames(X)),
#     pct_IncMSE = pi$pct,
#     pct_sd = pi$pct_sd,
#     used_by_a_tree = pi$used,
#     rho = apply(X, 2, function(v) suppressWarnings(cor(v, y, method = "spearman"))),
#     shap_dir = if (is.null(sd_)) NA_real_ else sd_$direction,
#     mean_relabun_pct = 100 * colMeans(X),
#     redox = rx) %>%
#     # colour comes from SHAP, falling back to the plain correlation if absent
#     mutate(direction_used = if (is.null(sd_)) rho else shap_dir,
#            direction = case_when(
#              is.na(direction_used) | abs(direction_used) < SHAP_DIR_MIN ~ "unclear",
#              direction_used > 0 ~ "higher with Mn",
#              TRUE ~ "lower with Mn")) %>%
#     arrange(desc(pct_IncMSE))
#   
#   INCMSE[[rx]] <- tab
#   
#   cat(sprintf("%-10s %d of %d taxa used by at least one tree, top %%IncMSE = %.1f%% (%s)\n",
#               rx, pi$n_used, ncol(X), tab$pct_IncMSE[1], tab$label[1]))
#   
# }
# 
# INCMSE <- bind_rows(INCMSE)
# 
# TOP <- INCMSE %>%
#   group_by(redox) %>% slice_max(pct_IncMSE, n = 15, with_ties = FALSE) %>%
#   mutate(rank = row_number()) %>% ungroup() %>%
#   mutate(key = paste(label, redox, sep = "__"))
# 
# fig_imp <- ggplot(TOP, aes(pct_IncMSE, reorder(key, pct_IncMSE), fill = direction)) +
#   geom_col(width = 0.75) +
#   geom_errorbarh(aes(xmin = pmax(0, pct_IncMSE - pct_sd), xmax = pct_IncMSE + pct_sd),
#                  height = 0.25, colour = "grey25", linewidth = 0.35) +
#   facet_wrap(~ redox, scales = "free") +
#   scale_y_discrete(labels = function(x) sub("__.*", "", x)) +
#   scale_fill_manual(values = DIR_COLS, name = NULL) +
#   labs(x = "Permuation Importance (%Inc MSE)", y = NULL) +
#        #caption = sprintf("bars: mean of %d permutations +- 1 SD; colour: direction of the taxon's effect on the prediction (TreeSHAP)",
#                          #N_PERM_REPEATS)) +
#   theme_pub + theme(plot.caption = element_text(size = 7.5, colour = "grey35", hjust = 0))
# save_fig(fig_imp, "fig_2_taxon_importance", 8.2, 4.6)

# The bar length says how much a genus matters. The colour says which way.
# 
# Green — "higher with Mn": in bottles where that genus is more abundant, the model predicts more dissolved Mn(II).
# Purple — "lower with Mn": where that genus is more abundant, the model predicts less.

#Version of the figure without SHAP coloring the bar direction

cat("\n--- %IncMSE (", N_PERM_REPEATS, " permutations per taxon) ---------------\n", sep = "")
INCMSE <- list()
for (rx in names(FITS)) {
  rf <- FITS[[rx]]; X <- DS[[rx]]$X; y <- DS[[rx]]$y
  pi <- pct_incmse(rf, X, y)
  
  tab <- tibble(
    taxon = colnames(X),
    label = short_label(colnames(X)),
    pct_IncMSE = pi$pct,
    pct_sd = pi$pct_sd,
    used_by_a_tree = pi$used,
    mean_relabun_pct = 100 * colMeans(X),
    redox = rx) %>%
    arrange(desc(pct_IncMSE))
  
  INCMSE[[rx]] <- tab
  
  cat(sprintf("%-10s %d of %d taxa used by at least one tree, top %%IncMSE = %.1f%% (%s)\n",
              rx, pi$n_used, ncol(X), tab$pct_IncMSE[1], tab$label[1]))
  
}
INCMSE <- bind_rows(INCMSE)
TOP <- INCMSE %>%
  group_by(redox) %>% slice_max(pct_IncMSE, n = 15, with_ties = FALSE) %>%
  mutate(rank = row_number()) %>% ungroup() %>%
  mutate(key = paste(label, redox, sep = "__"))

TOP <- TOP %>%
  mutate(redox = recode(redox, oxidation = "Aerobic", reduction = "Anaerobic"))

label_lookup <- setNames(TOP$label, TOP$key)
taxon_lookup <- setNames(TOP$taxon, TOP$key)

taxon_labels <- function(x) {
  labels <- unname(label_lookup[x])
  taxa <- unname(taxon_lookup[x])
  
  # Extract text between f__ and the next semicolon
  families <- sub(".*;f__", "", taxa)
  families <- sub(";.*", "", families)
  
  is_genus <- grepl("\\(Genus\\)$", labels)
  
  labels[is_genus] <- paste0(
    sub(" \\(Genus\\)$", "", labels[is_genus]),
    " (Genus, Family = ",
    families[is_genus],
    ")"
  )
  
  labels
}

fig_imp <- ggplot(TOP, aes(pct_IncMSE, reorder(key, pct_IncMSE))) +
  geom_col(width = 0.75, fill = "#C9A0DC") +
  geom_errorbarh(
    aes(
      xmin = pmax(0, pct_IncMSE - pct_sd),
      xmax = pct_IncMSE + pct_sd
    ),
    height = 0.25,
    colour = "grey25",
    linewidth = 0.35
  ) +
  facet_wrap(~redox, scales = "free") +
  scale_y_discrete(labels = taxon_labels) +
  labs(x = "Permutation Importance (%Inc MSE)", y = NULL) +
  theme_pub +
  theme(
    plot.caption = element_text(
      size = 7.5,
      colour = "grey35",
      hjust = 0
    )
  )
save_fig(fig_imp, "fig_2_taxon_importance", 12, 5)

# =============================================================================
# 9. Nested CV (cross-validation)
# =============================================================================

nested_cv <- function(ds, n_trials = N_TRIALS_INNER) {
  X <- ds$X; y <- ds$y; n <- nrow(X)
  oof <- rep(NA_real_, n); folds <- vector("list", n)
  
  for (i in seq_len(n)) {
    tr <- setdiff(seq_len(n), i)
    sp <- search_params(X[tr, , drop = FALSE], y[tr], n_trials, seed = SEED + i)
    f <- fit_rf(X[tr, , drop = FALSE], y[tr], sp$params, seed = SEED)
    oof[i] <- as.numeric(predict(f, X[i, , drop = FALSE]))
    folds[[i]] <- as_tibble(sp$params) %>%
      mutate(fold = i, inner_oob_r2 = sp$value,
             mtry = max(1L, as.integer(round(sp$params$mtry_frac * ncol(X)))))
  }
  
  ss_res <- sum((y - oof)^2); ss_tot <- sum((y - mean(y))^2)
  list(oof = oof, folds = bind_rows(folds),
       r2 = 1 - ss_res / ss_tot,
       rmse = sqrt(mean((y - oof)^2)),
       spearman = suppressWarnings(cor(y, oof, method = "spearman")))
}

PERF <- bind_rows(PERF)

cat("\n--- nested CV (16 folds, ", N_TRIALS_INNER, " Optuna trials inside each) ---\n", sep = "")
NESTED <- list(); NESTPRED <- list(); FOLDS <- list()

for (rx in names(DS)) {
  nc <- nested_cv(DS[[rx]])
  tuned_oob <- PERF %>% filter(redox == rx) %>% pull(oob_r2)
  NESTED[[rx]] <- tibble(
    redox = rx, nested_r2 = nc$r2, nested_rmse = nc$rmse, nested_spearman = nc$spearman,
    tuned_oob_r2 = tuned_oob, selection_optimism = tuned_oob - nc$r2,
    n_folds = nrow(DS[[rx]]$X), inner_trials = N_TRIALS_INNER)
  
  FOLDS[[rx]] <- nc$folds %>% mutate(redox = rx)
  NESTPRED[[rx]] <- tibble(redox = rx, observed = DS[[rx]]$y, predicted = nc$oof,
                           water_type = DS[[rx]]$meta$water_type, r2 = nc$r2)
  
  cat(sprintf("%-10s nested R2 = %.3f (RMSE %.2f, Spearman %.2f), OOB was %.3f, optimism %+.3f\n",
              rx, nc$r2, nc$rmse, nc$spearman, tuned_oob, tuned_oob - nc$r2))
  cat(sprintf(" per-fold mtry ranged %d-%d, nodesize %d-%d, sample_frac %.2f-%.2f\n",
              min(nc$folds$mtry), max(nc$folds$mtry),
              min(nc$folds$nodesize), max(nc$folds$nodesize),
              min(nc$folds$sample_frac), max(nc$folds$sample_frac)))
  
}

NESTED <- bind_rows(NESTED)
NESTPRED <- bind_rows(NESTPRED)
FOLDS <- bind_rows(FOLDS)

lab_nest <- NESTPRED %>% distinct(redox, r2) %>%
  mutate(txt = sprintf("nested R2 = %.2f", r2))

fig_nest <- ggplot(NESTPRED, aes(observed, predicted)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
  geom_point(aes(colour = water_type), size = 2.4, alpha = 0.9) +
  geom_text(data = lab_nest, aes(x = -Inf, y = Inf, label = txt),
            hjust = -0.15, vjust = 1.6, size = 3.2, inherit.aes = FALSE) +
  facet_wrap(~ redox, scales = "free") +
  scale_colour_brewer(palette = "Dark2", name = NULL) +
  labs(x = "observed Mn(II) (mg/L)", y = "nested out-of-fold prediction (mg/L)",
       caption = "each point predicted by a model whose hyperparameters never saw it") +
  theme_pub + theme(plot.caption = element_text(size = 7.5, colour = "grey35", hjust = 0))
save_fig(fig_nest, "fig_1_nested_cv", 7.2, 3.6)

# =============================================================================
# 10. Tables
# =============================================================================

t1 <- NESTED %>%
  transmute(Redox = redox, Cells = n_folds,
            Taxa = PERF$n_taxa[match(redox, PERF$redox)],
            `OOB R2` = tuned_oob_r2,
            `Nested R2` = nested_r2,
            `Nested RMSE (mg/L)` = nested_rmse,
            `Nested Spearman` = nested_spearman,
            `Selection optimism` = selection_optimism)
save_table(t1, "table_1_model_performance")

fold_rng <- FOLDS %>% group_by(redox) %>%
  summarise(mtry_lo = min(mtry), mtry_hi = max(mtry), .groups = "drop")

t2 <- PERF %>%
  transmute(Redox = redox, Trees = N_TREES, Sampler = "Optuna TPE",
            Trials = N_TRIALS, Objective = "out-of-bag R2",
            `mtry searched` = sprintf("1-%d", n_taxa),
            `mtry selected` = mtry,
            `Node size searched` = sprintf("1-%d", pmax(2, ceiling(0.2 * n))),
            `Node size selected` = nodesize,
            `Sample fraction searched` = "0.40-0.90",
            `Sample fraction selected` = sample_frac,
            Replacement = "without") %>%
  left_join(fold_rng %>% transmute(Redox = redox,
                                   `mtry across nested folds` = sprintf("%d-%d", mtry_lo, mtry_hi)), by = "Redox")
save_table(t2, "table_2_hyperparameters")

# t3 <- TOP %>%
#   transmute(Redox = redox, Rank = rank, Taxon = label,
#             `%IncMSE` = pct_IncMSE, `SD` = pct_sd,
#             Direction = direction, `Spearman rho` = rho,
#             `Mean relative abundance (%)` = mean_relabun_pct, Taxonomy = taxon)
# save_table(t3, "table_3_top_taxa")

# save_table(INCMSE %>%
#              transmute(Redox = redox, Taxon = label, `%IncMSE` = pct_IncMSE, SD = pct_sd,
#                        Direction = direction, `SHAP direction` = shap_dir,
#                        `Spearman rho` = rho, `Used by a tree` = used_by_a_tree,
#                        `Mean relative abundance (%)` = mean_relabun_pct, Taxonomy = taxon),
#            "table_S1_taxon_importance_all", digits = 4)

save_table(FOLDS %>%
             transmute(Redox = redox, Fold = fold, mtry, `mtry fraction` = mtry_frac,
                       `Node size` = nodesize, `Sample fraction` = sample_frac,
                       `Inner OOB R2` = inner_oob_r2),
           "table_S2_nested_cv_folds")

cat("\n--- model performance ----------------------------------------------\n")
print(as.data.frame(t1), row.names = FALSE, digits = 3)
cat("\n--- tuned hyperparameters ------------------------------------------\n")
print(as.data.frame(t2), row.names = FALSE, digits = 3)
cat("\nFigures -> ", FIGDIR, "/ Tables -> ", TABDIR, "/\n", sep = "")
print(t2)

# =============================================================================
# 11. Plotting Top 15 Significant Taxa as a Function of Mn (II) Concentrations
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)
library(NatParksPalettes)

# -----------------------------------------------------------------------------
# Choose the taxa to plot
# -----------------------------------------------------------------------------
# A random forest does not provide classical p-values. Here, "important" means
# highest %IncMSE. Set USE_TOP_N <- FALSE to plot every taxon in INCMSE.
USE_TOP_N <- TRUE
N_TOP <- 15L

if (USE_TOP_N) {
  taxa_by_redox <- list(
    oxidation = TOP %>%
      filter(redox %in% c("oxidation", "Aerobic")) %>%
      slice_max(order_by = pct_IncMSE, n = N_TOP, with_ties = FALSE) %>%
      pull(taxon),
    reduction = TOP %>%
      filter(redox %in% c("reduction", "Anaerobic")) %>%
      slice_max(order_by = pct_IncMSE, n = N_TOP, with_ties = FALSE) %>%
      pull(taxon)
  )
} else {
  taxa_by_redox <- list(
    oxidation = INCMSE %>%
      filter(redox %in% c("oxidation", "Aerobic")) %>%
      pull(taxon) %>% unique(),
    reduction = INCMSE %>%
      filter(redox %in% c("reduction", "Anaerobic")) %>%
      pull(taxon) %>% unique()
  )
}

taxa_by_redox <- lapply(taxa_by_redox, unique)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
water_labels <- c(
  GW       = "Unamended Groundwater",
  GW_Media = "Carbon-amended Groundwater",
  SW       = "Unamended SWIFT Water",
  SW_Media = "Carbon-amended SWIFT Water"
)

# Use one Death Valley colour per taxon, consistently across all four panels.
make_taxon_palette <- function(taxa) {
  taxa <- unique(taxa)
  n <- length(taxa)
  if (n == 0L) return(character(0))
  
  cols <- natparks.pals(
    name = "DeathValley",
    n = n,
    type = "continuous"
  )
  
  setNames(cols, taxa)
}

mn_red_size_limits <- mn_red %>%
  pull(Mn_mg_L) %>%
  .[is.finite(.)] %>%
  range(na.rm = TRUE)

mn_ox_size_limits <- mn_ox %>%
  pull(Mn_mg_L) %>%
  .[is.finite(.)] %>%
  range(na.rm = TRUE)

# -----------------------------------------------------------------------------
# Create one panel
# -----------------------------------------------------------------------------
make_single_mn_plot <- function(df, taxa, water_type_name, redox_label,
                                palette, size_limits, point_range = c(2.5, 8)) {
  taxa_present <- intersect(taxa, names(df))
  if (length(taxa_present) == 0L) return(NULL)
  
  plot_df <- df %>%
    filter(water_type == water_type_name) %>%
    select(water_type, day, Mn_mg_L, all_of(taxa_present)) %>%
    pivot_longer(
      cols = all_of(taxa_present),
      names_to = "taxon",
      values_to = "relative_abundance"
    ) %>%
    mutate(
      relative_abundance = relative_abundance * 100
    ) %>%
    filter(
      is.finite(day),
      is.finite(Mn_mg_L),
      is.finite(relative_abundance)
    ) %>%
    arrange(taxon, day)
  
  if (nrow(plot_df) == 0L) return(NULL)
  
  taxon_labels <- setNames(short_label(taxa_present), taxa_present)
  
  ggplot(
    plot_df,
    aes(
      x = day,
      y = relative_abundance,
      colour = taxon,
      group = taxon,
      size = Mn_mg_L
    )
  ) +
    geom_line(linewidth = 0.45, alpha = 0.65) +
    geom_point(alpha = 0.85) +
    scale_x_continuous(
      name = "Day",
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    scale_y_continuous(
      name = "Relative abundance (%)",
      labels = label_number(accuracy = 0.1, suffix = "%")
    ) +
    scale_size_continuous(
      name = "Mn (II) Concentration (mg/L)",
      limits = size_limits,
      trans = "sqrt",
      range = point_range
    ) +
    scale_colour_manual(
      values = palette,
      breaks = taxa_present,
      labels = taxon_labels,
      name = "Taxon",
      drop = FALSE
    ) +
    labs(
      title = water_labels[[water_type_name]]
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        linewidth = 0.25,
        colour = "grey92"
      ),
      plot.title = element_text(face = "bold", size = 10),
      legend.position = "right",
      legend.text = element_text(size = 9),
      legend.key.height = unit(0.35, "cm")
    )
}

# -----------------------------------------------------------------------------
# Build a four-panel figure
# -----------------------------------------------------------------------------
make_mn_abundance_panel <- function(df, taxa, redox_label, file_stub, size_limits, point_range = c(2.5, 8),
                                    output_dir = "rf_taxa_plots") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  palette <- make_taxon_palette(taxa)
  
  plots <- list(
    GW_Media = make_single_mn_plot(
      df, taxa, "GW_Media", redox_label, palette, size_limits, point_range
    ),
    SW_Media = make_single_mn_plot(
      df, taxa, "SW_Media", redox_label, palette, size_limits, point_range
    ),
    GW = make_single_mn_plot(
      df, taxa, "GW", redox_label, palette, size_limits, point_range
    ),
    SW = make_single_mn_plot(
      df, taxa, "SW", redox_label, palette, size_limits, point_range
    )
  )
  
  plots <- Filter(Negate(is.null), plots)
  
  p <- wrap_plots(plots, ncol = 2) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste(redox_label),
      theme = theme(
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 9, colour = "grey30", hjust = 0.5)
      )
    ) &
    theme(legend.position = "right")
  
  ggsave(
    file.path(output_dir, paste0(file_stub, ".png")),
    p, width = 13, height = 9, dpi = 300
  )
  ggsave(
    file.path(output_dir, paste0(file_stub, ".pdf")),
    p, width = 13, height = 9
  )
  
  p
}

# -----------------------------------------------------------------------------
# Generate the two four-panel figures
# -----------------------------------------------------------------------------
aerobic_mn_plot <- make_mn_abundance_panel(
  df = mn_ox,
  taxa = taxa_by_redox$oxidation,
  redox_label = "Aerobic treatments",
  file_stub = "fig_3_aerobic_taxa_vs_Mn",
  size_limits = mn_ox_size_limits,
  point_range = c(1.2, 4)
)

anaerobic_mn_plot <- make_mn_abundance_panel(
  df = mn_red,
  taxa = taxa_by_redox$reduction,
  redox_label = "Anaerobic treatments",
  file_stub = "fig_4_anaerobic_taxa_vs_Mn",
  size_limits = mn_red_size_limits,
  point_range = c(2.5, 8)
)

print(aerobic_mn_plot)
print(anaerobic_mn_plot)

cat("\nAerobic taxa plotted:\n")
print(short_label(taxa_by_redox$oxidation))

cat("\nAnaerobic taxa plotted:\n")
print(short_label(taxa_by_redox$reduction))


# =============================================================================
# 12. Relative Abundance Bar Plots of Taxa Identified by the Model
# =============================================================================

library(ggtext)
library(patchwork)
library(pals)

# Helper: extract genus, collapse unclassified to "Other"
extract_genus <- function(tax_string) {
  genus <- case_when(
    grepl("g__", tax_string) ~ str_extract(tax_string, "(?<=g__)([^;]+)"),
    TRUE ~ NA_character_
  ) %>% str_trim()
  
  case_when(
    is.na(genus)                    ~ "Other",
    genus == ""                     ~ "Other",
    grepl("^\\s*$", genus)          ~ "Other",
    grepl("^[0-9]", genus)          ~ "Other",   # unnamed/placeholder clades
    genus %in% c("uncultured", "NA") ~ "Other",
    TRUE                            ~ genus
  )
}

# Build global palette from your top-15 RF taxa (from TOP)
build_global_palette <- function(otu_data, top_taxa) {
  
  all_genera <- otu_data %>%
    select(all_of(top_taxa)) %>%
    pivot_longer(everything(), names_to = "Taxon", values_to = "Abundance") %>%
    mutate(Genus = extract_genus(Taxon)) %>%
    group_by(Genus) %>%
    summarise(total = sum(Abundance), .groups = "drop") %>%
    arrange(desc(total)) %>%
    filter(Genus != "Other") %>%
    pull(Genus)
  
  # Other always at bottom (first in rev = bottom of stacked bar)
  genus_order <- c("Other", rev(all_genera))
  
  n_col <- length(all_genera)
  kelly_colors <- unname(kelly(n_col + 1)[2:(n_col + 1)])
  pal <- setNames(kelly_colors, all_genera)
  pal["Other"] <- "#D3D3D3"
  
  list(top_taxa = top_taxa, genus_order = genus_order, palette = pal)
}

# Pull top-15 taxa per redox group straight from TOP
aerobic_top_taxa   <- TOP %>% filter(redox == "Aerobic")   %>% pull(taxon)
anaerobic_top_taxa <- TOP %>% filter(redox == "Anaerobic") %>% pull(taxon)

aerobic_global   <- build_global_palette(mn_ox,  aerobic_top_taxa)
anaerobic_global <- build_global_palette(mn_red, anaerobic_top_taxa)

# Check top taxa against rf model output
print(aerobic_top_taxa)
print(anaerobic_top_taxa)

# Core plot function (no y limit — set later per pair)
plot_rf_taxa <- function(otu_subset,
                         top_taxa,
                         genus_order,
                         palette,
                         treatment_name,
                         y_max = 1.0) {
  
  if (is.null(otu_subset) || nrow(otu_subset) == 0) return(NULL)
  
  meta_cols <- c("redox", "water_type", "day", "date")
  
  plot_df <- otu_subset %>%
    select(any_of(meta_cols), all_of(top_taxa)) %>%
    pivot_longer(
      cols      = all_of(top_taxa),
      names_to  = "Taxon",
      values_to = "Abundance"
    ) %>%
    mutate(
      Genus = extract_genus(Taxon),
      Genus = factor(Genus, levels = genus_order),
      day   = factor(day, levels = sort(unique(as.numeric(as.character(day)))))
    )
  
  ggplot(plot_df, aes(x = day, y = Abundance, fill = Genus)) +
    geom_bar(stat = "identity", position = "stack", width = 0.85) +
    scale_fill_manual(
      values = palette,
      name   = "Genus",
      drop   = FALSE
    ) +
    scale_y_continuous(
      limits = c(0, y_max),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(title = treatment_name, x = "Day", y = "Relative Abundance") +
    theme_bw(base_size = 16) +
    theme(
      axis.text.x     = element_text(angle = 0, hjust = 1, size = 16),
      axis.text.y     = element_text(size = 16),
      axis.title      = element_text(size = 16),
      plot.title      = element_text(size = 18, face = "bold"),
      legend.text     = element_text(size = 16),
      legend.title    = element_text(size = 16),
      legend.key.size = unit(0.5, "cm"),
      panel.grid      = element_blank()
    ) +
    guides(fill = guide_legend(nrow = 3, byrow = TRUE))
}

# Helper: compute shared y_max for a pair of water types 
get_y_max <- function(otu_data, water_types, top_taxa) {
  otu_data %>%
    filter(water_type %in% water_types) %>%
    select(water_type, day, all_of(top_taxa)) %>%
    pivot_longer(all_of(top_taxa), names_to = "Taxon", values_to = "Abundance") %>%
    group_by(water_type, day) %>%
    summarise(total = sum(Abundance), .groups = "drop") %>%
    pull(total) %>%
    max(na.rm = TRUE) %>%
    { min(., 1.0) * 1.05 }  # 5% headroom, cap at 1
}

# Treatment labels
aerobic_labels <- c(
  GW       = "Aerobic Unamended Groundwater",
  GW_Media = "Aerobic Carbon-amended Groundwater",
  SW       = "Aerobic Unamended SWIFT Water",
  SW_Media = "Aerobic Carbon-amended SWIFT Water"
)

anaerobic_labels <- c(
  GW       = "Anaerobic Unamended Groundwater",
  GW_Media = "Anaerobic Carbon-amended Groundwater",
  SW       = "Anaerobic Unamended SWIFT Water",
  SW_Media = "Anaerobic Carbon-amended SWIFT Water"
)

# Generate all plots 
make_plot_group <- function(otu_data, global, labels, redox_name) {
  
  top_taxa    <- global$top_taxa
  genus_order <- global$genus_order
  palette     <- global$palette
  
  y_max_unamended <- get_y_max(otu_data, c("GW", "SW"),             top_taxa)
  y_max_amended   <- get_y_max(otu_data, c("GW_Media", "SW_Media"), top_taxa)
  
  plots <- list(
    GW_Media = plot_rf_taxa(
      otu_subset     = otu_data %>% filter(water_type == "GW_Media"),
      top_taxa       = top_taxa,
      genus_order    = genus_order,
      palette        = palette,
      treatment_name = labels["GW_Media"],
      y_max          = y_max_amended
    ),
    SW_Media = plot_rf_taxa(
      otu_subset     = otu_data %>% filter(water_type == "SW_Media"),
      top_taxa       = top_taxa,
      genus_order    = genus_order,
      palette        = palette,
      treatment_name = labels["SW_Media"],
      y_max          = y_max_amended
    ),
    GW = plot_rf_taxa(
      otu_subset     = otu_data %>% filter(water_type == "GW"),
      top_taxa       = top_taxa,
      genus_order    = genus_order,
      palette        = palette,
      treatment_name = labels["GW"],
      y_max          = y_max_unamended
    ),
    SW = plot_rf_taxa(
      otu_subset     = otu_data %>% filter(water_type == "SW"),
      top_taxa       = top_taxa,
      genus_order    = genus_order,
      palette        = palette,
      treatment_name = labels["SW"],
      y_max          = y_max_unamended
    )
  )
  
  Filter(Negate(is.null), plots)
}

aerobic_plots   <- make_plot_group(mn_ox,  aerobic_global,   aerobic_labels,   "Aerobic")
anaerobic_plots <- make_plot_group(mn_red, anaerobic_global, anaerobic_labels, "Anaerobic")

# Combined panels: amended top row, unamended bottom row, shared legend on bottom
make_combined_panel <- function(plots, title) {
  
  top_row    <- list(plots[["GW_Media"]], plots[["SW_Media"]])
  bottom_row <- list(plots[["GW"]],       plots[["SW"]])
  
  top_row    <- Filter(Negate(is.null), top_row)
  bottom_row <- Filter(Negate(is.null), bottom_row)
  
  all_plots <- c(top_row, bottom_row)
  
  wrap_plots(all_plots, ncol = 2) +
    plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") +
    plot_annotation(
      title = title,
      theme = theme(
        plot.title  = element_text(size = 16, face = "bold", hjust = 0.5)
      )
    ) &
    theme(legend.position = "bottom") &
    guides(fill = guide_legend(nrow = 3, byrow = TRUE))
}

aerobic_combined   <- make_combined_panel(aerobic_plots,   "Top 15 RF Taxa — Aerobic Treatments")
anaerobic_combined <- make_combined_panel(anaerobic_plots, "Top 15 RF Taxa — Anaerobic Treatments")

# Save
output_dir <- "rf_taxa_plots"
dir.create(output_dir, showWarnings = FALSE)

for (wt in names(aerobic_plots)) {
  ggsave(file.path(output_dir, paste0("aerobic_", wt, "_rf_taxa.pdf")),
         aerobic_plots[[wt]], width = 10, height = 7)
  ggsave(file.path(output_dir, paste0("aerobic_", wt, "_rf_taxa.png")),
         aerobic_plots[[wt]], width = 10, height = 7, dpi = 300)
}
for (wt in names(anaerobic_plots)) {
  ggsave(file.path(output_dir, paste0("anaerobic_", wt, "_rf_taxa.pdf")),
         anaerobic_plots[[wt]], width = 10, height = 7)
  ggsave(file.path(output_dir, paste0("anaerobic_", wt, "_rf_taxa.png")),
         anaerobic_plots[[wt]], width = 10, height = 7, dpi = 300)
}

ggsave(file.path(output_dir, "aerobic_combined.pdf"),
       aerobic_combined, width = 18, height = 14)
ggsave(file.path(output_dir, "aerobic_combined.png"),
       aerobic_combined, width = 18, height = 14, dpi = 300)

ggsave(file.path(output_dir, "anaerobic_combined.pdf"),
       anaerobic_combined, width = 18, height = 14)
ggsave(file.path(output_dir, "anaerobic_combined.png"),
       anaerobic_combined, width = 18, height = 14, dpi = 300)

cat("All plots saved to:", output_dir, "\n")
