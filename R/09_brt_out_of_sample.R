# =============================================================================
# 09_brt_out_of_sample.R
# Out-of-sample accuracy of the BRTs, tested on transects that were not used to fit the model.
#
# Two validation schemes:
#   random_70_30        20 random 70/30 train/test splits, stratified by Exposure x Year
#   leave_one_site_out  each site held out in turn
# Each model is fitted on the training transects only (same settings as 04 and 05, including the
# predictor selection step for species), then used to predict the held-out transects.
# Metrics: RMSE, MAE, predictive R2 (1 - SSE/SST, relative to the training mean) and Spearman rho.
#
# Produces:  outputs/brt_out_of_sample/BRT_out_of_sample_table.csv (test-set RMSE and R2 per response)
#            BRT_out_of_sample_summary_by_response.csv, BRT_out_of_sample_species_overview.csv,
#            BRT_out_of_sample_all_splits.csv
# Runtime:   many hours for all species. Test first with SPECIES_LIMIT <- 5.
# =============================================================================

source(file.path("R", "00_setup.R"))
source(file.path("R", "helpers_brt.R"))
require_pkgs(c("gbm", "dismo"))

# ---- Settings ------------------------------------------------------------------
REPEATS          <- 20     # number of random 70/30 splits
TRAIN_FRAC       <- 0.70
SPECIES_LIMIT    <- NA     # e.g. 5 to test quickly; NA = all species
MIN_NONZERO      <- 10     # species must be present in at least this many transects
RUN_TOTALS       <- FALSE  # TRUE also validates total density, richness and biomass
INCLUDE_EXPOSURE <- TRUE   # species models use the 21 coral genera + Exposure
SEED             <- 123

res_dir <- "brt_out_of_sample"
set.seed(SEED)

dat <- read_lizard(); dat$Exposure <- factor(dat$Exposure)
fish_species <- col_range(dat, FISH_FIRST, FISH_LAST)
if (!is.na(SPECIES_LIMIT)) fish_species <- head(fish_species, SPECIES_LIMIT)

# ---- helpers ------------------------------------------------------------------
pred_brt <- function(fit, newdata) as.numeric(predict(fit, newdata, n.trees = fit$gbm.call$best.trees, type = "response"))

rel_inf_vars <- function(fit) { s <- summary(fit, plotit = FALSE); as.character(s$var[s$rel.inf > 0]) }

pois_dev <- function(y, mu) { mu <- pmax(mu, 1e-9); 2 * sum(ifelse(y > 0, y * log(y / mu), 0) - (y - mu)) }

# Metrics for predictions `p` of observations `y`; `ref` = mean of the TRAINING response (null model)
metrics <- function(y, p, ref, family) {
  sse <- sum((y - p)^2); sst <- sum((y - ref)^2)
  out <- tibble(n = length(y), RMSE = sqrt(mean((y - p)^2)), MAE = mean(abs(y - p)),
                R2_pred = if (sst > 0) 1 - sse / sst else NA_real_,
                Spearman = suppressWarnings(if (length(unique(p)) > 1 && length(unique(y)) > 1) cor(y, p, method = "spearman") else NA_real_))
  if (family == "poisson") out$Dev_explained <- { d0 <- pois_dev(y, rep(ref, length(y))); if (d0 > 0) 1 - pois_dev(y, p) / d0 else NA_real_ }
  out
}

# One train/test evaluation. Returns a tibble (Set = "train"/"test") or NULL.
# `response_fun` / `inverse` handle the log1p transform used for the total metrics.
eval_split <- function(df, response, predictors, family, train_idx, test_idx, two_stage) {
  tr <- df[train_idx, , drop = FALSE]; te <- df[test_idx, , drop = FALSE]
  fit <- fit_brt(tr, response, predictors, family); if (is.null(fit)) return(NULL)
  used <- predictors
  if (two_stage) {                                   # reduced model, selected from TRAINING data only
    used <- rel_inf_vars(fit); used <- used[used %in% predictors]
    if (INCLUDE_EXPOSURE && "Exposure" %in% predictors) used <- unique(c(used, "Exposure"))   # Exposure always kept
    if (length(used) < 2) return(NULL)
    tr2 <- tr[, c(response, used), drop = FALSE]
    fit <- fit_brt(tr2, response, used, family); if (is.null(fit)) return(NULL)
  }
  ref <- mean(tr[[response]])
  bind_rows(mutate(metrics(tr[[response]], pred_brt(fit, tr[, used, drop = FALSE]), ref, family), Set = "train", pred = NA),
            mutate(metrics(te[[response]], pred_brt(fit, te[, used, drop = FALSE]), ref, family), Set = "test", pred = NA)) |>
    mutate(n_trees = fit$gbm.call$best.trees) |> dplyr::select(-pred) |>
    (\(m) list(metrics = m, test_obs = te[[response]], test_pred = pred_brt(fit, te[, used, drop = FALSE]), ref = ref))()
}

# Stratified 70/30 split (strata = Exposure x Year) -> list of train/test indices
random_splits <- function(df, n, frac) {
  strata <- interaction(df$Exposure, df$Year, drop = TRUE)
  map(seq_len(n), \(i) { tr <- unlist(lapply(split(seq_len(nrow(df)), strata), \(ix) if (length(ix) < 2) ix else sample(ix, max(1, round(frac * length(ix)))))); list(train = sort(tr), test = setdiff(seq_len(nrow(df)), tr)) })
}
site_splits <- function(df) map(unique(df$Site), \(s) list(train = which(df$Site != s), test = which(df$Site == s)))

run_response <- function(response, predictors, family, transform = identity, back = identity, two_stage = FALSE, label = response) {
  df <- dat[, unique(c(response, predictors, "Site", "Year", "Exposure")), drop = FALSE]
  df <- df[stats::complete.cases(df), ]   # Exposure is also used to stratify the splits
  df[[response]] <- transform(as.numeric(df[[response]])); df$Exposure <- factor(df$Exposure, levels = levels(dat$Exposure))
  # A. repeated random 70/30
  sp <- random_splits(df, REPEATS, TRAIN_FRAC)
  A <- map_dfr(seq_along(sp), \(i) { r <- eval_split(df, response, predictors, family, sp[[i]]$train, sp[[i]]$test, two_stage)
                                     if (is.null(r)) NULL else mutate(r$metrics, Scheme = "random_70_30", Split = i) })
  # B. leave-one-site-out (pooled out-of-fold predictions)
  ss <- site_splits(df); obs <- c(); prd <- c(); refs <- c(); B_train <- list()
  for (s in ss) { r <- eval_split(df, response, predictors, family, s$train, s$test, two_stage)
                  if (!is.null(r)) { obs <- c(obs, r$test_obs); prd <- c(prd, r$test_pred); B_train[[length(B_train) + 1]] <- dplyr::filter(r$metrics, Set == "train")
                                     refs <- c(refs, rep(r$ref, length(r$test_obs))) } }
  B <- if (length(obs)) bind_rows(mutate(metrics(obs, prd, refs, family), Set = "test", n_trees = NA_real_),
                                  mutate(summarise(bind_rows(B_train), across(c(RMSE, MAE, R2_pred, Spearman, n_trees), mean), n = mean(n)), Set = "train")) |>
        mutate(Scheme = "leave_one_site_out", Split = NA_integer_) else NULL
  res <- bind_rows(A, B)
  if (nrow(res) == 0) return(NULL)
  mutate(res, Response = label, Family = family, Scale = if (identical(transform, identity)) "response" else "log1p(response)")
}

# ---- 1. Total fish metrics (gaussian, log1p, as in Fig 4) ------------------------
pred_total <- c(intersect(c(CORAL_PREDICTORS, "Coral.genera"), names(dat)), "Exposure")
res_total <- if (!RUN_TOTALS) NULL else map_dfr(intersect(c("Total.density", "Species.richness", "Total.biomass"), names(dat)), \(r) {
  message("Total metric: ", r); run_response(r, pred_total, "gaussian", transform = log1p, two_stage = FALSE) })

# ---- 2. Species (poisson, two-stage, 21 coral genera + Exposure, as in Fig 5) --------------------------------
pred_sp <- if (INCLUDE_EXPOSURE) c(CORAL_PREDICTORS, "Exposure") else CORAL_PREDICTORS   # species models
keep_sp <- fish_species[vapply(fish_species, \(s) sum(as.numeric(dat[[s]]) > 0, na.rm = TRUE) >= MIN_NONZERO, logical(1))]
message(length(keep_sp), " species pass the screen (present in >= ", MIN_NONZERO, " transects)")
res_sp <- map_dfr(keep_sp, \(s) { message("Species: ", s); run_response(s, pred_sp, "poisson", two_stage = TRUE) })

all_res <- bind_rows(res_total, res_sp)
write.csv(all_res, out_path(res_dir, "BRT_out_of_sample_all_splits.csv"), row.names = FALSE)

# ---- 3. Summary tables -------------------------------------------------------------
summ <- all_res |> group_by(Response, Family, Scheme, Set) |>
  summarise(n_splits = n(), across(c(RMSE, MAE, R2_pred, Spearman), list(mean = \(x) mean(x, na.rm = TRUE), sd = \(x) sd(x, na.rm = TRUE)), .names = "{.col}_{.fn}"),
            Dev_explained_mean = if ("Dev_explained" %in% names(all_res)) mean(Dev_explained, na.rm = TRUE) else NA_real_, .groups = "drop")
write.csv(summ, out_path(res_dir, "BRT_out_of_sample_summary_by_response.csv"), row.names = FALSE)

manuscript <- summ |> dplyr::filter(Set == "test") |>
  mutate(`RMSE (mean ± SD)` = ifelse(is.na(RMSE_sd), sprintf("%.2f", RMSE_mean), sprintf("%.2f ± %.2f", RMSE_mean, RMSE_sd)),
         `Predictive R2 (mean ± SD)` = ifelse(is.na(R2_pred_sd), sprintf("%.2f", R2_pred_mean), sprintf("%.2f ± %.2f", R2_pred_mean, R2_pred_sd))) |>
  dplyr::select(Response, Scheme, `RMSE (mean ± SD)`, `Predictive R2 (mean ± SD)`)
write.csv(manuscript, out_path(res_dir, "BRT_out_of_sample_table.csv"), row.names = FALSE)

# Species overview: how many species are predicted better than the training-mean null (R2_pred > 0)?
sp_overview <- summ |> dplyr::filter(Set == "test", Family == "poisson") |> group_by(Scheme) |>
  summarise(species = n(), median_R2_pred = median(R2_pred_mean, na.rm = TRUE), q25 = quantile(R2_pred_mean, .25, na.rm = TRUE),
            q75 = quantile(R2_pred_mean, .75, na.rm = TRUE), n_better_than_null = sum(R2_pred_mean > 0, na.rm = TRUE))
write.csv(sp_overview, out_path(res_dir, "BRT_out_of_sample_species_overview.csv"), row.names = FALSE)
print(manuscript); print(sp_overview)
message("Done. Results in ", file.path(OUT_DIR, res_dir))
