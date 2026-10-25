# =============================================================================
# 07c_temporal_slopes_inla_figureS3.R
# Figure S3: how does each fish-coral slope change through the disturbance-recovery cycle?
# For every fish species x coral genus pair that passes the same minimum-data screening as script 06,
# and in each exposure zone, an INLA model estimates the fish-coral slope in three periods:
#   2011 (pre-disturbance), 2015/2017 (disturbance; average of the two years' slopes), 2020 (post-disturbance).
# The slopes are drawn as one caterpillar plot per exposure zone, and each pair is given a pattern
# (e.g. Neutral, Positive-neutral, Negative-positive-negative) from the sign of its three slopes.
#
# Model (one per pair and exposure zone):
#   fish ~ 1 + Year + coral_2011 + coral_2015 + coral_2017 + coral_2020 + f(SiteTransect, iid), family nbinomial
#   coral_<year> = coral cover standardised within the model (mean 0, SD 1), set to 0 in other years,
#                  so each coefficient is the fish-coral slope in that year
#   SiteTransect = one ID per transect position (site, depth stratum and transect number)
#   Slopes are on the log scale per 1% of coral cover (the standardised slope divided by the SD of coral cover)
#   Response: fish density (individuals per 1000 m2; integer-valued, counts x 4 or 20), as in script 06
#
# Which pairs are fitted (settings below; the first block is identical to script 06):
#   every fish species x coral genus in the data, with
#     all zones together:  >= MIN_N_OBS observations, >= MIN_NONZERO non-zero fish records, >= 2 coral values
#     within the zone:     >= MIN_NONZERO_LEVEL non-zero records, >= MIN_ROWS_LEVEL observations, >= 2 coral values
#     within each period:  >= MIN_NONZERO_PERIOD non-zero records, >= MIN_ROWS_LEVEL observations, >= 2 coral values
#
# Outputs (outputs/temporal_slopes_inla/):
#   period_slopes_long.csv                  slope, 95% interval and symbol for each period of every model
#   pattern_by_pair.csv                     one row per pair and zone: symbols, pattern, stable or not
#   pattern_summary_by_zone.csv             number and percentage of pairs in each pattern, per zone
#   screening_summary.csv                   how many pairs pass each rule, per zone
#   skipped_pairs.csv                       pairs not fitted, and why
#   comparison_with_hand_compiled.csv       pairs also in the original hand-compiled Figure S3 (if the file is present)
#   Figure_S3a_..._exposed.pdf, S3b oblique, S3c lagoon, S3d sheltered
# Run time: hours (fitted models are not stored; each model's slope draws are saved, so a stopped run can be restarted).
# Needs:    INLA (not on CRAN: https://www.r-inla.org/download-install), ggplot2
# =============================================================================

# Sets the working directory to the project folder (the folder holding R and data) if it isn't already.
# This works when the script is open in RStudio; otherwise run setwd("path to lizard_repo") first.
if (!file.exists("R/00_setup.R")) setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
source("R/00_setup.R")
if (!exists("TEST_NO_INLA")) { require_pkgs(c("INLA")); library(INLA) }
library(ggplot2)

res_dir   <- "temporal_slopes_inla"
cache_dir <- dirname(out_path(res_dir, "slope_draws", ".keep"))

# ---- Settings ------------------------------------------------------------------
FAMILY            <- "nbinomial"
MIN_N_OBS         <- 8      # as script 06
MIN_NONZERO       <- 10     # as script 06
MIN_CORAL_UNIQUE  <- 2      # as script 06
MIN_NONZERO_LEVEL <- 5      # as script 06 (here: per exposure zone)
MIN_ROWS_LEVEL    <- 3      # as script 06
MIN_CORAL_UNIQUE_LEVEL <- 2 # as script 06
MIN_NONZERO_PERIOD <- 3     # new: non-zero records needed in each of the three periods within the zone
NSAMP             <- 1000   # posterior samples per model
EXPOSURES <- c("Exposed", "Oblique", "Lagoon", "Sheltered")   # order of Figures S3a-d
set.seed(SEED_MASTER)

# ---- 1. Data ----------------------------------------------------------------
lizard <- read_lizard() |>
  mutate(Year = factor(Year), Exposure = as.character(Exposure),
         SiteTransect = interaction(Site, Depth, Transect, drop = TRUE))
fish_cols  <- col_range(lizard, FISH_FIRST, FISH_LAST)
coral_cols <- intersect(CORAL_PREDICTORS, names(lizard))
period_of  <- c("2011" = "2011", "2015" = "2015/2017", "2017" = "2015/2017", "2020" = "2020")
PERIODS    <- c("2011", "2015/2017", "2020")

# ---- 2. Screening (same rules as script 06, applied to each pair and zone) -----------
screen <- function(fish, coral) {
  d <- tibble(Exposure = lizard$Exposure, period = unname(period_of[as.character(lizard$Year)]),
              y = as.integer(round(lizard[[fish]])), x = suppressWarnings(as.numeric(lizard[[coral]]))) |>
    filter(!is.na(y), !is.na(x))
  ok_all <- nrow(d) >= MIN_N_OBS && sum(d$y > 0) >= MIN_NONZERO && n_distinct(d$x) >= MIN_CORAL_UNIQUE
  purrr::map_dfr(EXPOSURES, function(z) {
    dz <- d |> filter(Exposure == z)
    pz <- dz |> group_by(period) |> summarise(n = n(), nz = sum(y > 0), u = n_distinct(x), .groups = "drop")
    reason <- if (!ok_all) "pair fails overall rules (non-zero records, coral values)"
      else if (nrow(dz) < MIN_ROWS_LEVEL || sum(dz$y > 0) < MIN_NONZERO_LEVEL || n_distinct(dz$x) < MIN_CORAL_UNIQUE_LEVEL)
        "zone fails per-zone rules"
      else if (nrow(pz) < 3 || any(pz$n < MIN_ROWS_LEVEL) || any(pz$nz < MIN_NONZERO_PERIOD) || any(pz$u < MIN_CORAL_UNIQUE_LEVEL))
        "a period fails per-period rules"
      else NA_character_
    tibble(fish = fish, coral = coral, zone = z, n_nonzero_zone = sum(dz$y > 0), skip_reason = reason)
  })
}
grid <- tidyr::crossing(fish = fish_cols, coral = coral_cols)
screening <- purrr::map2_dfr(grid$fish, grid$coral, screen)
write.csv(dplyr::filter(screening, !is.na(skip_reason)), out_path(res_dir, "skipped_pairs.csv"), row.names = FALSE)
screening_summary <- screening |> group_by(zone) |>
  summarise(pairs_considered = n(), pairs_fitted = sum(is.na(skip_reason)),
            fail_overall = sum(grepl("overall", skip_reason)), fail_zone = sum(grepl("per-zone", skip_reason)),
            fail_period = sum(grepl("per-period", skip_reason)), .groups = "drop")
write.csv(screening_summary, out_path(res_dir, "screening_summary.csv"), row.names = FALSE)
print(screening_summary)
jobs <- screening |> filter(is.na(skip_reason))
message(nrow(jobs), " fish x coral x zone models to fit")

# ---- 3. Helpers ------------------------------------------------------------------
hdi95 <- function(x, mass = 0.95) {
  x <- sort(x); n <- length(x); k <- floor(mass * n)
  i <- which.min(x[(k + 1):n] - x[1:(n - k)])
  c(lower = x[i], upper = x[i + k])
}
symbol <- function(lower, upper) ifelse(lower > 0, "+", ifelse(upper < 0, "-", "n"))
# Pattern from the three symbols: repeated neighbours are merged, e.g. "+nn" -> Positive-neutral, "n+n" -> Neutral-positive-neutral
pattern_name <- function(s) {
  s <- s[c(TRUE, s[-1] != s[-length(s)])]
  paste(c("+" = "Positive", "-" = "Negative", "n" = "Neutral")[s], collapse = "-") |> tolower() |>
    (\(x) paste0(toupper(substring(x, 1, 1)), substring(x, 2)))()
}

# Slope draws (NSAMP x 3 periods, per 1% coral cover, log scale) for one pair in one zone; NULL if the model fails
fit_slopes <- function(fish, coral, zone) {
  d <- lizard |> filter(Exposure == zone) |>
    transmute(Year, SiteTransect, y = as.integer(round(.data[[fish]])), coral = suppressWarnings(as.numeric(.data[[coral]]))) |>
    filter(!is.na(y), !is.na(coral)) |> droplevels()
  d$SiteTransect_ID <- as.integer(factor(d$SiteTransect))
  sc <- scale(d$coral); sd_coral <- attr(sc, "scaled:scale"); d$coral_sc <- as.numeric(sc)
  yrs <- levels(d$Year)
  for (y in yrs) d[[paste0("coral_", y)]] <- ifelse(d$Year == y, d$coral_sc, 0)
  slope_terms <- paste0("coral_", yrs)
  form <- as.formula(paste("y ~ 1 + Year +", paste(slope_terms, collapse = " + "), "+ f(SiteTransect_ID, model = 'iid')"))
  mod <- tryCatch(inla(form, data = d, family = FAMILY, control.predictor = list(compute = FALSE),
                       control.compute = list(config = TRUE), control.inla = list(strategy = "adaptive")),
                  error = function(e) e)
  if (inherits(mod, "error")) { message("  INLA error: ", conditionMessage(mod)); return(NULL) }
  samp <- tryCatch(inla.posterior.sample(NSAMP, mod), error = function(e) e)
  if (inherits(samp, "error")) { message("  sampling error: ", conditionMessage(samp)); return(NULL) }
  lat  <- sub(":[0-9]+$", "", rownames(samp[[1]]$latent))
  rows <- vapply(slope_terms, function(t) which(lat == t)[1], integer(1))
  if (anyNA(rows)) { message("  slope terms not found"); return(NULL) }
  B <- do.call(rbind, lapply(samp, \(s) s$latent[rows, 1])) / sd_coral
  colnames(B) <- yrs
  cbind(`2011` = B[, "2011"],
        `2015/2017` = rowMeans(B[, intersect(c("2015", "2017"), yrs), drop = FALSE]),
        `2020` = B[, "2020"])
}

# ---- 4. Fit every selected pair in every zone (draws are cached) ------------------------
results <- vector("list", nrow(jobs))
for (k in seq_len(nrow(jobs))) {
  j <- jobs[k, ]
  f <- file.path(cache_dir, paste0(safe_file(paste(j$fish, j$coral, j$zone)), ".rds"))
  if (file.exists(f)) draws <- readRDS(f) else {
    message("Model ", k, " / ", nrow(jobs), ": ", j$fish, " ~ ", j$coral, " (", j$zone, ")")
    draws <- fit_slopes(j$fish, j$coral, j$zone)
    saveRDS(draws, f)
  }
  if (is.null(draws)) next
  h <- apply(draws, 2, hdi95)
  results[[k]] <- tibble(coral = j$coral, fish = j$fish, zone = j$zone, period = colnames(draws),
                         slope = colMeans(draws), lower95 = h["lower", ], upper95 = h["upper", ]) |>
    mutate(symbol = symbol(lower95, upper95))
}
long <- bind_rows(results) |>
  mutate(zone = factor(zone, EXPOSURES), period = factor(period, PERIODS))
write.csv(long, out_path(res_dir, "period_slopes_long.csv"), row.names = FALSE)

# ---- 5. Pattern of each pair ------------------------------------------------------------
by_pair <- long |> arrange(zone, coral, fish, period) |>
  group_by(zone, coral, fish) |>
  summarise(symbols = paste(symbol, collapse = " "), pattern = pattern_name(symbol), .groups = "drop") |>
  mutate(stable = pattern %in% c("Neutral", "Positive", "Negative"))
write.csv(by_pair, out_path(res_dir, "pattern_by_pair.csv"), row.names = FALSE)
pattern_summary <- by_pair |> count(zone, pattern) |> group_by(zone) |>
  mutate(percent = round(100 * n / sum(n), 1)) |> ungroup()
write.csv(pattern_summary, out_path(res_dir, "pattern_summary_by_zone.csv"), row.names = FALSE)
stable_by_zone <- by_pair |> group_by(zone) |> summarise(pairs = n(), n_stable = sum(stable), percent_stable = round(100 * n_stable / pairs, 1))
write.csv(stable_by_zone, out_path(res_dir, "stable_by_zone.csv"), row.names = FALSE)
print(stable_by_zone)

# Comparison with the original hand-compiled figure, if its file is present
hand_file <- file.path("data", "figure_S3_hand_compiled_slopes.csv")
if (file.exists(hand_file)) {
  hand <- read.csv(hand_file, stringsAsFactors = FALSE) |> distinct(zone, coral, fish, hand_class)
  cmp <- inner_join(by_pair, hand, by = c("zone", "coral", "fish")) |> mutate(same_pattern = pattern == hand_class)
  write.csv(cmp, out_path(res_dir, "comparison_with_hand_compiled.csv"), row.names = FALSE)
  message(nrow(cmp), " pairs are in both; ", round(100 * mean(cmp$same_pattern), 1), "% have the same pattern")
}

# ---- 6. Figure S3a-d ------------------------------------------------------------------------
fig_names <- c(Exposed = "Figure_S3a_Fish-coral_slope_comparison_exposed.pdf",
               Oblique = "Figure_S3b_Fish-coral_slope_comparison_oblique.pdf",
               Lagoon  = "Figure_S3c_Fish-coral_slope_comparison_lagoon.pdf",
               Sheltered = "Figure_S3d_Fish-coral_slope_comparison_sheltered.pdf")
for (z in EXPOSURES) {
  d <- long |> filter(zone == z)
  if (!nrow(d)) next
  d <- d |> mutate(fish_lab = gsub("\\.", " ", fish)) |>
    mutate(coral = factor(coral, levels = rev(sort(unique(coral)))))
  n_pairs <- n_distinct(d$coral, d$fish)
  p <- ggplot(d, aes(x = slope, y = fish_lab)) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0, linewidth = 0.3) +
    geom_point(size = 0.8) +
    facet_grid(coral ~ period, scales = "free_y", space = "free_y", switch = "y") +
    labs(x = "Slope (log scale, per 1% coral cover)", y = "Fish",
         title = paste0("Fish-coral slopes through time: ", tolower(z), " zone")) +
    theme_classic(base_size = 8) +
    theme(axis.text.y = element_text(face = "italic", size = 4),
          strip.text.y.left = element_text(angle = 0, face = "italic", hjust = 1),
          strip.background = element_blank(), strip.placement = "outside",
          panel.spacing.y = unit(0.15, "lines"), plot.title = element_text(size = 9))
  ggsave(out_path(res_dir, fig_names[[z]]), p, width = 8, height = max(5, 1.6 + 0.09 * n_pairs + 0.15 * n_distinct(d$coral)),
         limitsize = FALSE)
}
message("Done. Results are in outputs/", res_dir, "/")
