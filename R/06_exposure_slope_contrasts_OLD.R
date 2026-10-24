# OLD VERSION of script 06 (density values, minimum 2 non-zero records). Kept only to reproduce the old results for script 12.
# =============================================================================
# 06_exposure_slope_contrasts.R
# Does the relationship between a fish species and a coral genus differ among exposure zones?
# One INLA negative-binomial model per fish species x coral genus; exposure-specific slopes and their
# pairwise contrasts are calculated from posterior samples.
#
# Outputs (outputs/inla_exposure_slopes/):
#   all_exposure_slope_contrasts.csv                              every fitted contrast
#   significant_exposure_slope_differences_descending.csv         contrasts whose 95% CrI excludes zero
#   skipped_pairs.csv                                             fish x coral pairs not fitted, and why
#   Exposure_contrasts_heatmap.pdf
#
# Model:    fish ~ 1 + Exposure + coral_sc + coral_sc:Exposure + f(SiteTransect, iid), family nbinomial
#           (SiteTransect = one ID per transect position: site, depth stratum and transect number)
#           (coral_sc = coral cover standardised within each model: mean 0, SD 1)
# Run time: hours. Fitted models are saved and re-used on later runs.
# Needs:    INLA (not on CRAN: https://www.r-inla.org/download-install)
# =============================================================================

# Sets the working directory to the project folder (the folder holding R and data) if it isn't already.
# This works when the script is open in RStudio; otherwise run setwd("path to lizard_repo") first.
if (!file.exists("R/00_setup.R")) setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
source("R/00_setup.R")
require_pkgs(c("INLA"))
library(INLA)

res_dir <- "inla_exposure_slopes_old"   # OLD thresholds; writes to its own folder so it cannot mix with the new run
model_dir <- dirname(out_path(res_dir, "models", ".keep"))

# ---- Settings (as in the Methods) ---------------------------------------------
FAMILY           <- "nbinomial"
MIN_N_OBS        <- 8     # minimum observations per fish x coral dataset
MIN_NONZERO      <- 2     # minimum non-zero fish counts
MIN_CORAL_UNIQUE <- 2     # minimum distinct coral cover values
MIN_ROWS_LEVEL   <- 3     # per exposure level: minimum observations ...
MIN_CORAL_UNIQUE_LEVEL <- 2   # ... and distinct coral cover values; failing levels are dropped
NSAMP            <- 400   # posterior samples per model
set.seed(SEED_MASTER)

# ---- 1. Data ----------------------------------------------------------------
lizard <- read_lizard() |>
  mutate(Site = factor(Site), Exposure = factor(Exposure, levels = EXPOSURE_LEVELS_GRADIENT),
         Transect = factor(Transect), SiteTransect = interaction(Site, Depth, Transect, drop = TRUE))

fish_cols  <- col_range(lizard, FISH_FIRST, FISH_LAST)
coral_cols <- intersect(CORAL_PREDICTORS, names(lizard))
message(length(fish_cols), " fish columns x ", length(coral_cols), " coral genera")

# ---- 2. One nested dataset per fish x coral pair ------------------------------
nested <- lizard |>
  select(Site, Exposure, Depth, Transect, SiteTransect, all_of(fish_cols), all_of(coral_cols)) |>
  pivot_longer(all_of(fish_cols),  names_to = "fish",  values_to = "fish_abundance") |>
  pivot_longer(all_of(coral_cols), names_to = "coral", values_to = "coral_abundance") |>
  group_by(coral, fish) |> nest() |> ungroup() |>
  mutate(label = paste0(fish, " _ ", coral),
         data = purrr::map(data, \(d) d |>
           mutate(fish_abundance = suppressWarnings(as.integer(fish_abundance)),
                  coral_abundance = suppressWarnings(as.numeric(coral_abundance)),
                  Exposure = factor(Exposure)) |>
           filter(!is.na(fish_abundance), !is.na(coral_abundance), !is.na(Exposure), !is.na(SiteTransect)) |>
           droplevels()))

# ---- 3. Screening: overall, then per exposure level ---------------------------
level_ok <- function(d) {
  d |> group_by(Exposure) |>
    summarise(n = n(), coral_unique = n_distinct(coral_abundance), .groups = "drop") |>
    mutate(ok = n >= MIN_ROWS_LEVEL & coral_unique >= MIN_CORAL_UNIQUE_LEVEL)
}

screen_pair <- function(d) {
  if (nrow(d) < MIN_N_OBS)                       return(list(data = d, reason = "n_obs < min_n_obs"))
  if (sum(d$fish_abundance > 0) < MIN_NONZERO)   return(list(data = d, reason = "fish nonzero < min_nonzero"))
  if (n_distinct(d$coral_abundance) < MIN_CORAL_UNIQUE) return(list(data = d, reason = "coral_unique < min_coral_unique"))
  if (nlevels(d$Exposure) < 2)                   return(list(data = d, reason = "<2 Exposure levels"))

  keep <- as.character(level_ok(d)$Exposure[level_ok(d)$ok])
  if (length(keep) < 2) return(list(data = d, reason = "<2 Exposure levels pass per-level screens"))
  d <- d |> filter(as.character(Exposure) %in% keep) |> droplevels()
  if ("Exposed" %in% levels(d$Exposure)) d$Exposure <- relevel(d$Exposure, ref = "Exposed")
  d$coral_sc <- as.numeric(scale(d$coral_abundance))
  if (all(!is.finite(d$coral_sc))) return(list(data = d, reason = "coral_sc not finite"))
  list(data = d, reason = NA_character_)
}

screened <- purrr::map(nested$data, screen_pair)
nested$data        <- purrr::map(screened, "data")
nested$skip_reason <- purrr::map_chr(screened, "reason")
message(sum(is.na(nested$skip_reason)), " of ", nrow(nested), " pairs pass screening")

# ---- 4. Fit INLA models (cached) ----------------------------------------------
fit_pair <- function(d, label) {
  file <- file.path(model_dir, paste0("model_", safe_file(label), ".rds"))
  if (file.exists(file)) return(readRDS(file))
  # Time saver: if the new run (script 06) already fitted this pair on exactly the same data, reuse that model
  new_file <- file.path(OUT_DIR, "inla_exposure_slopes", "models", paste0("model_", safe_file(label), ".rds"))
  if (file.exists(new_file)) {
    m <- tryCatch(readRDS(new_file), error = function(e) NULL)
    nd <- tryCatch(m$.args$data, error = function(e) NULL)
    if (!is.null(nd) && nrow(nd) == nrow(d) && isTRUE(all.equal(as.numeric(nd$fish_abundance), as.numeric(d$fish_abundance))) &&
        isTRUE(all.equal(as.numeric(nd$coral_sc), as.numeric(d$coral_sc)))) { saveRDS(m, file); return(m) }
  }
  mod <- tryCatch(
    inla(fish_abundance ~ 1 + Exposure + coral_sc + coral_sc:Exposure + f(SiteTransect, model = "iid"),
         data = d, family = FAMILY, control.predictor = list(compute = FALSE),
         control.compute = list(config = TRUE), control.inla = list(strategy = "adaptive")),
    error = function(e) e)
  if (inherits(mod, "error")) { message("INLA error for ", label, ": ", conditionMessage(mod)); return(NULL) }
  saveRDS(mod, file); mod
}

# ---- 5. Exposure-specific slopes and pairwise contrasts ------------------------
# slope(base level) = beta_coral;  slope(other level) = beta_coral + beta_coral:level
contrast_pair <- function(mod, d) {
  fixed_names <- rownames(mod$summary.fixed)
  if (!("coral_sc" %in% fixed_names)) return(list(tab = NULL, reason = "coral_sc coefficient missing"))

  samp <- tryCatch(inla.posterior.sample(NSAMP, mod), error = function(e) e)
  if (inherits(samp, "error")) return(list(tab = NULL, reason = paste("posterior.sample error:", conditionMessage(samp))))

  lat <- rownames(samp[[1]]$latent)
  rows <- lat[grepl(paste(fixed_names, collapse = "|"), lat, fixed = FALSE)]
  if (length(rows) != length(fixed_names))
    return(list(tab = NULL, reason = "could not match fixed effects in posterior samples"))
  B <- do.call(rbind, lapply(samp, \(s) s$latent[rows, 1]))
  colnames(B) <- fixed_names

  levs <- levels(d$Exposure); base <- levs[1]
  int_col <- function(level) {
    cand <- c(paste0("coral_sc:Exposure", level), paste0("Exposure", level, ":coral_sc"))
    cand <- cand[cand %in% colnames(B)]
    if (length(cand)) cand[1] else NA_character_
  }
  slopes <- setNames(vector("list", length(levs)), levs)
  slopes[[base]] <- B[, "coral_sc"]
  for (L in setdiff(levs, base)) { nm <- int_col(L); if (!is.na(nm)) slopes[[L]] <- B[, "coral_sc"] + B[, nm] }
  slopes <- slopes[!vapply(slopes, is.null, logical(1))]
  if (length(slopes) < 2) return(list(tab = NULL, reason = "<2 slopes estimable (interaction dropped)"))

  tab <- purrr::map_dfr(utils::combn(names(slopes), 2, simplify = FALSE), function(pr) {
    dr <- slopes[[pr[1]]] - slopes[[pr[2]]]
    tibble(exposure_levels_used = paste(levs, collapse = ";"), contrast = paste0(pr[1], "_minus_", pr[2]),
           mean = mean(dr), sd = sd(dr), q025 = as.numeric(quantile(dr, 0.025)), q975 = as.numeric(quantile(dr, 0.975)),
           prob_gt0 = mean(dr > 0), n_obs = nrow(d), n_nonzero = sum(d$fish_abundance > 0))
  })
  list(tab = tab, reason = NA_character_)
}

ok <- which(is.na(nested$skip_reason))
results <- vector("list", nrow(nested)); reasons <- nested$skip_reason
for (k in seq_along(ok)) {
  i <- ok[k]
  if (k %% 25 == 0) message("Pair ", k, " / ", length(ok))
  mod <- fit_pair(nested$data[[i]], nested$label[i])
  if (is.null(mod)) { reasons[i] <- "INLA error"; next }
  cp <- contrast_pair(mod, nested$data[[i]])
  results[[i]] <- cp$tab; if (!is.na(cp$reason)) reasons[i] <- cp$reason
}

contrasts_all <- purrr::map2_dfr(seq_len(nrow(nested)), results, \(i, tab)
  if (is.null(tab)) NULL else mutate(tab, fish = nested$fish[i], coral = nested$coral[i], .before = 1))
write.csv(contrasts_all, out_path(res_dir, "all_exposure_slope_contrasts.csv"), row.names = FALSE)
write.csv(tibble(fish = nested$fish, coral = nested$coral, reason = reasons) |> filter(!is.na(reason)),
          out_path(res_dir, "skipped_pairs.csv"), row.names = FALSE)

# ---- 6. Contrasts whose 95% CrI excludes zero ----------------------------------
significant <- contrasts_all |>
  mutate(significant = q025 > 0 | q975 < 0, abs_difference = abs(mean),
         posterior_support = pmax(prob_gt0, 1 - prob_gt0),
         direction = case_when(mean > 0 ~ "first exposure in contrast has steeper slope",
                               mean < 0 ~ "second exposure in contrast has steeper slope", TRUE ~ "no difference")) |>
  filter(significant) |> arrange(desc(abs_difference)) |>
  select(fish, coral, contrast, exposure_levels_used, mean, abs_difference, sd, q025, q975,
         prob_gt0, posterior_support, n_obs, n_nonzero, direction)
write.csv(significant, out_path(res_dir, "significant_exposure_slope_differences_descending.csv"), row.names = FALSE)

# ---- 7. Heatmap of all slope contrasts -----------------------------------------
if (nrow(contrasts_all) > 0) {
  heat <- contrasts_all |> mutate(fish_id = match(fish, sort(unique(fish)))) |>
    ggplot(aes(fish_id, coral, fill = mean)) + geom_tile() + facet_wrap(~ contrast) +
    scale_fill_gradient2(low = "red", mid = "white", high = "green", name = "Slope\ndifference") +
    scale_x_continuous(breaks = NULL) + labs(x = "Fish species (alphabetical)", y = NULL) +
    theme_minimal() + theme(axis.text.y = element_text(size = 6, face = "italic"))
  ggsave(out_path(res_dir, "Exposure_contrasts_heatmap.pdf"), heat, width = 12, height = 8)
}
message("Done: ", nrow(contrasts_all), " contrasts (", nrow(significant), " with 95% CrI excluding zero)")
