# =============================================================================
# 05_brt_species.R
# Boosted regression trees for each fish species: relative influence of the coral genera and
# exposure on species density, and the direction of each coral effect.
#
# Produces:  Figure 5  outputs/brt_species/Figure5_tall, Figure5_part1, Figure5_part2 (.pdf, .png)
#            Table S4  outputs/brt_species/BRT_Relative_Influence_Results_selected_corals.csv
#            Table S4 (formatted)  outputs/brt_species/Table_S4_species_relative_influence.csv
# Runtime:   many hours. Each species is saved in outputs/brt_species/cache/, so a stopped run can be restarted.
#
# For each species:
#   1. Fit a BRT (poisson) with the 21 coral genera and Exposure.
#   2. Refit using the predictors with relative influence > 0 (Exposure is always kept).
#   3. Refit on 100 bootstrap samples to get the mean and 95% CI of relative influence.
#   4. The direction of each coral effect is the slope of its partial-dependence curve.
#   5. Signed influence = bootstrap mean x direction, standardised across all species and predictors.
# =============================================================================

# Sets the working directory to the project folder (the folder holding R and data) if it isn't already.
# This works when the script is open in RStudio; otherwise run setwd("path to lizard_repo") first.
if (!file.exists("R/00_setup.R")) setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
source("R/00_setup.R")
source(file.path("R", "helpers_brt.R"))
source(file.path("R", "helpers_figure5.R"))
require_pkgs(c("gbm", "dismo", "ggplot2", "dplyr"))

res_dir <- "brt_species"
cache_dir <- dirname(out_path(res_dir, "cache", ".keep"))
set.seed(123)

# Quick test: Sys.setenv(BRT_SPECIES_LIMIT = "5") runs the first 5 species only;
# Sys.setenv(N_BOOT_OVERRIDE = "10") reduces the bootstrap refits.
species_limit <- as.integer(Sys.getenv("BRT_SPECIES_LIMIT", unset = NA))
if (nzchar(Sys.getenv("N_BOOT_OVERRIDE"))) N_BOOT <- as.integer(Sys.getenv("N_BOOT_OVERRIDE"))
n_cores <- max(1, parallel::detectCores() - 1)

# Only species recorded in at least this many transects are modelled
MIN_TRANSECTS_PRESENT <- 10
INCLUDE_EXPOSURE <- TRUE   # predictors = 21 coral genera + Exposure

lizard <- read_lizard(); lizard$Exposure <- factor(lizard$Exposure)
fish_species <- col_range(lizard, FISH_FIRST, FISH_LAST)
n_present <- vapply(fish_species, function(s) sum(suppressWarnings(as.numeric(lizard[[s]])) > 0, na.rm = TRUE), numeric(1))
fish_species <- fish_species[n_present >= MIN_TRANSECTS_PRESENT]
message(length(fish_species), " species recorded in >= ", MIN_TRANSECTS_PRESENT, " transects")
if (!is.na(species_limit)) fish_species <- head(fish_species, species_limit)

predictors_sp <- if (INCLUDE_EXPOSURE) c(CORAL_PREDICTORS, "Exposure") else CORAL_PREDICTORS
stopifnot(all(predictors_sp %in% names(lizard)))

fast_slope <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 2) return(NA_real_)
  vx <- var(x)
  if (!is.finite(vx) || vx == 0) return(NA_real_)
  cov(x, y) / vx
}

fit_species <- function(species) {
  cache_file <- file.path(cache_dir, paste0(safe_file(species), ".rds"))
  if (file.exists(cache_file)) return(readRDS(cache_file))

  message("BRT for ", species)
  df <- lizard[, c(species, predictors_sp), drop = FALSE]
  df <- df[stats::complete.cases(df), , drop = FALSE]
  df[[species]] <- suppressWarnings(as.numeric(df[[species]]))

  res <- NULL
  if (nrow(df) >= 10 && sum(df[[species]], na.rm = TRUE) > 0) {
    full <- fit_brt(df, species, predictors_sp, "poisson")
    if (!is.null(full)) {
      rel_full <- extract_rel_inf(full, predictors_sp)
      important <- rel_full$var[rel_full$rel.inf > 0]
      if (INCLUDE_EXPOSURE) important <- unique(c(important, "Exposure"))
      df_f <- df[, c(species, important), drop = FALSE]
      filt <- if (length(important) >= 2) fit_brt(df_f, species, important, "poisson") else NULL   # gbm needs >= 2 predictors
      if (!is.null(filt)) {
        rel <- extract_rel_inf(filt, important)
        ci  <- bootstrap_rel_inf(df_f, species, important, "poisson", n_boot = N_BOOT, n_cores = n_cores)
        m   <- match(rel$var, ci$var)

        direction <- setNames(rep(NA_real_, length(important)), important)
        for (pred in setdiff(important, "Exposure")) {
          pd <- tryCatch(gbm::plot.gbm(filt, i.var = pred, return.grid = TRUE), error = function(e) NULL)
          if (!is.null(pd) && nrow(pd) > 1) direction[pred] <- fast_slope(pd[[1]], exp(pd$y))
        }
        boot_mean <- ci$boot_mean[m]
        signed <- boot_mean * sign(direction[rel$var]); signed[is.na(signed)] <- 0

        res <- data.frame(Species = species, Predictor = rel$var, Influence = rel$rel.inf,
                          BootMean = boot_mean, Lower95 = ci$lower[m], Upper95 = ci$upper[m],
                          NSuccess = ci$n_success[m], Effect_Direction = direction[rel$var],
                          SignedInfluence = signed, stringsAsFactors = FALSE, row.names = NULL)
      }
    }
  }
  saveRDS(res, cache_file)   # NULL is cached too, so failed species are not retried every run
  res
}

all_res <- lapply(fish_species, fit_species)
mat <- do.call(rbind, all_res[!vapply(all_res, is.null, logical(1))])
message(length(unique(mat$Species)), " of ", length(fish_species), " species with converged models")

# Signed influence is standardised across all species and predictors
mat$SignedInfluence <- as.numeric(scale(mat$SignedInfluence))
write.csv(mat, out_path(res_dir, "BRT_Relative_Influence_Results_selected_corals.csv"), row.names = FALSE)
write.csv(species_influence_table(mat), out_path(res_dir, "Table_S4_species_relative_influence.csv"), row.names = FALSE)

# ---- Figure 5 ------------------------------------------------------------------
draw_figure5(mat, file.path(OUT_DIR, res_dir))
message("Done: Figure 5 and Table S4 written to ", file.path(OUT_DIR, res_dir))
