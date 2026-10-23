# =============================================================================
# 07_temporal_slopes_brms.R
# How does each fish-coral relationship change through the disturbance-recovery cycle?
# For every fish species x coral genus pair in data/fish_coral_pairs_figure6.csv, and in each exposure zone,
# a Bayesian model (brms) estimates the fish-coral slope in three periods:
#   2011 (pre-disturbance), 2015/2017 (disturbance; the average of the two years' slopes), 2020 (post-disturbance).
# The three slopes are turned into a trajectory class, which gives the table used to draw Figure 6.
#
# Model (one per pair and exposure zone):   fish density ~ coral_sc * Year + (1 | SiteTransect)
#   fish density = the density values in the data file (counts per 1000 m2), zero-inflated Poisson
#   coral_sc     = coral cover standardised within each model (mean 0, SD 1)
#   SiteTransect = site and transect number
#   Priors: intercept normal(-2.3, 1.5); coefficients normal(0, 1); random-effect SD student-t(3, 0, 1.5); zero-inflation beta(1, 1)
# Script 07b repeats the analysis on raw counts with a negative binomial model, as a check.
#
# Produces (outputs/temporal_slopes/):
#   fish_coral_relationship_change_table.csv   the table read by script 08 (same layout as the original table)
#   period_slopes_long.csv                     slope, 95% interval and symbol for each period of every model
#   comparison_with_original_table.csv         cell-by-cell comparison with the table in data/
#   threshold_sensitivity.csv                  cells retained at stricter minimum-data rules (for the Methods / reviewer reply)
#   unclassified_patterns.csv                  patterns not covered by the rules below (check by hand)
# Run time: many hours (88 pairs x 4 exposure zones). Each model's slopes are saved, so a stopped run can be restarted.
# Needs:    brms, emmeans
# =============================================================================

# Sets the working directory to the project folder (the folder holding R and data) if it isn't already.
# This works when the script is open in RStudio; otherwise run setwd("path to lizard_repo") first.
if (!file.exists("R/00_setup.R")) setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
source("R/00_setup.R")
require_pkgs(c("brms", "emmeans", "readr"))
library(brms); library(emmeans); library(readr)

# ---- Settings ------------------------------------------------------------------
MIN_NONZERO_ZONE   <- 5                 # a model needs at least this many records of the species in the exposure zone (as in script 06)
MIN_NONZERO_PERIOD <- 1                 # ... and at least this many in each of the three periods
SPARSE_PERIOD      <- 3                 # cells with fewer records than this in any period are flagged as 'sparse' (see sensitivity file)
CHAINS <- 4; ITER <- 3000; WARMUP <- 1000
SEED   <- SEED_MASTER
EXPOSURES    <- c("Lagoon", "Sheltered", "Oblique", "Exposed")
TRAJECTORIES <- c("Stable", "Stronger", "Weaker", "Change", "Bounce back")

res_dir   <- "temporal_slopes"
cache_dir <- dirname(out_path(res_dir, "model_slopes", ".keep"))
set.seed(SEED)

# ---- 1. Data ------------------------------------------------------------------
pairs <- read.csv(file.path("data", "fish_coral_pairs_figure6.csv"), stringsAsFactors = FALSE)
lizard <- read_lizard()
lizard <- lizard |>
  mutate(Year = factor(Year), Exposure = as.character(Exposure),
         SiteTransect = paste(Site, Transect))
stopifnot(all(pairs$fish %in% names(lizard)), all(pairs$coral %in% names(lizard)))

# ---- 2. Helpers ------------------------------------------------------------------
# 95% highest-density interval of a vector of posterior draws
hdi95 <- function(x, mass = 0.95) {
  x <- sort(x); n <- length(x); k <- floor(mass * n)
  i <- which.min(x[(k + 1):n] - x[1:(n - k)])
  c(lower = x[i], upper = x[i + k])
}
# "+" if the interval is above zero, "-" if below zero, "n" if it spans zero
symbol <- function(x) { h <- hdi95(x); if (h["lower"] > 0) "+" else if (h["upper"] < 0) "-" else "n" }

template_fit <- NULL
# Slopes of the three periods (matrix of posterior draws) for one pair in one exposure zone; NULL if it cannot be fitted
fit_slopes <- function(fish, coral, zone) {
  d <- lizard |> filter(Exposure == zone) |>
    transmute(Year, SiteTransect, y = as.integer(round(.data[[fish]])), coral = as.numeric(.data[[coral]])) |>
    filter(!is.na(y), !is.na(coral)) |> droplevels()
  period <- c("2011" = "2011", "2015" = "2015/2017", "2017" = "2015/2017", "2020" = "2020")[as.character(d$Year)]
  nz <- tapply(d$y > 0, period, sum)
  if (sum(d$y > 0) < MIN_NONZERO_ZONE || any(nz < MIN_NONZERO_PERIOD) || length(unique(d$coral)) < 3) return(NULL)
  d$coral_sc <- as.numeric(scale(d$coral))
  priors <- prior(normal(-2.3, 1.5), class = "Intercept") + prior(normal(0, 1), class = "b") +
            prior(student_t(3, 0, 1.5), class = "sd") + prior(beta(1, 1), class = "zi")
  fam <- zero_inflated_poisson()
  # The Stan model is compiled once (first fit) and reused for every later pair, which saves most of the run time
  fit <- tryCatch({
    if (is.null(template_fit)) {
      f1 <- brm(y ~ coral_sc * Year + (1 | SiteTransect), data = d, family = fam, prior = priors,
                chains = CHAINS, iter = ITER, warmup = WARMUP, cores = min(CHAINS, parallel::detectCores()),
                control = list(adapt_delta = 0.95, max_treedepth = 15), seed = SEED, refresh = 0, silent = 2)
      template_fit <<- f1; f1
    } else {
      update(template_fit, newdata = d, chains = CHAINS, iter = ITER, warmup = WARMUP,
             cores = min(CHAINS, parallel::detectCores()), seed = SEED, refresh = 0, silent = 2)
    }},
    error = function(e) { message("  brms failed: ", conditionMessage(e)); NULL })
  if (is.null(fit)) return(NULL)
  em <- emtrends(fit, ~ Year, var = "coral_sc")                    # slope for each year, link scale
  draws <- as.matrix(coda::as.mcmc(em))
  colnames(draws) <- levels(d$Year)
  cbind(`2011` = draws[, "2011"], `2015/2017` = rowMeans(draws[, c("2015", "2017")]), `2020` = draws[, "2020"])
}

# Trajectory class from the three symbols (and, for sign-consistent patterns, the change in slope size)
#   Stable       the same symbol in all three periods, and no clear change in the size of the slope
#   Stronger     same sign throughout and the slope gets larger (|slope 2020| > |slope 2011|, interval excludes 0),
#                or a neutral slope that becomes + or -
#   Weaker       same sign throughout and the slope gets smaller, or a + or - slope that becomes neutral
#   Change       the slope is + in 2011 and - in 2020, or the reverse ("+ to -", "- to +")
#   Bounce back  the same symbol in 2011 and 2020 but a different one in 2015/2017 (e.g. "+ - +", "n + n")
classify <- function(draws) {
  s <- apply(draws, 2, symbol); s1 <- s[1]; s2 <- s[2]; s3 <- s[3]
  if (s1 == s2 && s2 == s3) {
    if (s1 == "n") return(list(class = "Stable", text = "n", symbols = s))
    size_change <- hdi95(abs(draws[, 3]) - abs(draws[, 1]))
    cls <- if (size_change["lower"] > 0) "Stronger" else if (size_change["upper"] < 0) "Weaker" else "Stable"
    return(list(class = cls, text = s1, symbols = s))
  }
  if (s1 == s3) return(list(class = "Bounce back", text = paste(s1, s2, s3), symbols = s))
  if (s1 != "n" && s3 != "n") return(list(class = "Change", text = paste(s1, "to", s3), symbols = s))
  if (s1 == "n") return(list(class = "Stronger", text = s3, symbols = s))
  if (s3 == "n") return(list(class = "Weaker", text = s1, symbols = s))
  list(class = NA_character_, text = paste(s, collapse = " "), symbols = s)
}

# ---- 3. Fit every pair in every exposure zone -------------------------------------
jobs <- tidyr::crossing(pairs, zone = EXPOSURES) |> arrange(match(paste(coral, fish), paste(pairs$coral, pairs$fish)), match(zone, EXPOSURES))
results <- vector("list", nrow(jobs))
for (k in seq_len(nrow(jobs))) {
  j <- jobs[k, ]
  f <- file.path(cache_dir, paste0(safe_file(paste(j$fish, j$coral, j$zone)), ".rds"))
  if (file.exists(f)) { draws <- readRDS(f) } else {
    message("Model ", k, " / ", nrow(jobs), ": ", j$fish, " ~ ", j$coral, " (", j$zone, ")")
    draws <- fit_slopes(j$fish, j$coral, j$zone)
    saveRDS(draws, f)                       # NULL is saved too, so pairs that cannot be fitted are not retried
  }
  if (is.null(draws)) next
  zd <- lizard[lizard$Exposure == j$zone, ]
  min_period_n <- min(tapply(zd[[j$fish]] > 0, c("2011" = "a", "2015" = "b", "2017" = "b", "2020" = "c")[as.character(zd$Year)], sum))
  cl <- classify(draws)
  h  <- apply(draws, 2, hdi95)
  results[[k]] <- tibble(coral = j$coral, fish = j$fish, zone = j$zone, period = colnames(draws),
                         slope = colMeans(draws), lower95 = h["lower", ], upper95 = h["upper", ],
                         symbol = cl$symbols, trajectory = cl$class, cell = cl$text,
                         zone_records = sum(zd[[j$fish]] > 0), min_period_records = min_period_n)
}
long <- bind_rows(results)
write.csv(long, out_path(res_dir, "period_slopes_long.csv"), row.names = FALSE)
# Sensitivity to the minimum-data rule: how many cells remain, and how many are not 'Stable', at stricter thresholds
sens <- long |> distinct(coral, fish, zone, trajectory, zone_records, min_period_records) |>
  tidyr::crossing(min_records_per_period = c(1, 3, 5)) |> filter(min_period_records >= min_records_per_period) |>
  group_by(min_records_per_period) |> summarise(cells = n(), not_stable = sum(trajectory != "Stable", na.rm = TRUE),
                                               percent_not_stable = round(100 * not_stable / cells, 1))
write.csv(sens, out_path(res_dir, "threshold_sensitivity.csv"), row.names = FALSE)
write.csv(filter(distinct(long, coral, fish, zone, trajectory, cell), is.na(trajectory)),
          out_path(res_dir, "unclassified_patterns.csv"), row.names = FALSE)

# ---- 4. Table in the layout read by script 08 --------------------------------------
cells <- long |> distinct(coral, fish, zone, trajectory, cell) |> filter(!is.na(trajectory))
layout <- expand.grid(trajectory = TRAJECTORIES, zone = EXPOSURES, stringsAsFactors = FALSE)   # 20 columns
rows <- pairs |> rowwise() |> mutate(vals = list(vapply(seq_len(nrow(layout)), \(i) {
  m <- cells$cell[cells$coral == coral & cells$fish == fish & cells$zone == layout$zone[i] & cells$trajectory == layout$trajectory[i]]
  if (length(m)) m[1] else ""
}, character(1)))) |> ungroup()
pretty <- \(x) gsub("\\.", " ", x)
body <- do.call(rbind, lapply(seq_len(nrow(rows)), \(i) c(
  if (i == 1 || rows$coral[i] != rows$coral[i - 1]) pretty(rows$coral[i]) else "", pretty(rows$fish[i]), rows$vals[[i]])))
header1 <- c("", "", unlist(lapply(EXPOSURES, \(e) c(e, rep("", length(TRAJECTORIES) - 1)))))
header2 <- c("Coral genus", "Fish", rep(TRAJECTORIES, length(EXPOSURES)))
write.table(rbind(header1, header2, body), out_path(res_dir, "fish_coral_relationship_change_table.csv"),
            sep = ",", row.names = FALSE, col.names = FALSE, na = "", qmethod = "double")

# ---- 5. Comparison with the original table ------------------------------------------
read_table_long <- function(file) {
  raw <- read_csv(file, col_names = FALSE, show_col_types = FALSE, col_types = cols(.default = "c"))
  zone <- tidyr::fill(tibble(z = as.character(unlist(raw[1, 3:22]))) |> mutate(z = na_if(z, "")), z)$z
  traj <- as.character(unlist(raw[2, 3:22]))
  b <- raw[-c(1, 2), 1:22]; names(b) <- c("coral", "fish", paste(zone, traj, sep = "|"))
  b$coral <- tidyr::fill(tibble(coral = na_if(b$coral, "")), coral)$coral
  b |> filter(!is.na(fish)) |> tidyr::pivot_longer(-c(coral, fish), values_drop_na = TRUE) |>
    tidyr::separate(name, c("zone", "trajectory"), sep = "\\|") |>
    transmute(coral = gsub(" ", ".", coral), fish = gsub(" ", ".", fish), zone, trajectory, cell = str_squish(value))
}
original <- read_table_long(CHANGE_TABLE_FILE)
cmp <- full_join(rename(original, original_trajectory = trajectory, original_cell = cell),
                 rename(distinct(cells, coral, fish, zone, trajectory, cell), new_trajectory = trajectory, new_cell = cell),
                 by = c("coral", "fish", "zone")) |>
  mutate(same_class = original_trajectory == new_trajectory, same_cell = same_class & original_cell == new_cell)
write.csv(cmp, out_path(res_dir, "comparison_with_original_table.csv"), row.names = FALSE)
message("Cells in the original table: ", nrow(original), "; in the new table: ", nrow(cells),
        "; present in both: ", sum(!is.na(cmp$original_cell) & !is.na(cmp$new_cell)),
        "; same trajectory class: ", sum(cmp$same_class, na.rm = TRUE), "; identical cell: ", sum(cmp$same_cell, na.rm = TRUE))
message("Done: table written to ", file.path(OUT_DIR, res_dir))
