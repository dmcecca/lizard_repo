# =============================================================================
# 11_depth_sensitivity.R
# Does pooling the shallow (3-5 m) and deep (8-10 m) transects matter?
# Refits the total-density, richness and biomass models with Depth and Depth x Exposure added,
# and compares them with the pooled models from script 01:
#   pooled model:  Response ~ Year * Exposure + (1 | Site) + (1 | SiteDepthTransect)
#   depth model:   Response ~ Year * Exposure + Depth * Exposure + (1 | Site) + (1 | SiteDepthTransect)
#
# Produces (outputs/depth_sensitivity/):
#   Depth_model_comparison      LOO comparison (elpd_diff of the depth model minus the pooled model)
#   Depth_effects_by_exposure   deep / shallow rate ratio [95% CrI] within each exposure zone and overall
# Run time: as long as script 01 (hours). The depth models are saved and re-used on later runs.
# Run script 01 first, so that the pooled models exist.
# =============================================================================

# Sets the working directory to the project folder (the folder holding R and data) if it isn't already.
# This works when the script is open in RStudio; otherwise run setwd("path to lizard_repo") first.
if (!file.exists("R/00_setup.R")) setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
source("R/00_setup.R")
source(file.path("R", "helpers_glmm.R"))
require_pkgs(c("flextable", "officer", "loo"))

res_dir <- "depth_sensitivity"
set.seed(SEED_MASTER)

dat         <- read_lizard()
df_base     <- prep_glmm_data(dat)
metric_cols <- col_range(dat, METRIC_FIRST, METRIC_LAST)
df_base$Depth <- factor(df_base$Depth)

pooled_dir <- file.path(OUT_DIR, "glmm_fish_metrics", "brms_fish_models_site")
depth_dir  <- dirname(out_path(res_dir, "brms_depth_models_site", ".keep"))

# ---- 1. Fit (or load) the depth models, using the same family as the pooled model ----
fit_depth_model <- function(resp_col) {
  file <- file.path(depth_dir, paste0(safe_file(resp_col), ".rds"))
  if (file.exists(file)) return(readRDS(file))
  df <- df_base |>
    transmute(Year, Exposure, Depth, Site = factor(Site), SiteDepthTransect, Response = as.numeric(.data[[resp_col]])) |>
    filter(!is.na(Response))
  f <- Response ~ Year * Exposure + Depth * Exposure + (1 | Site) + (1 | SiteDepthTransect)
  if (all(df$Response > 0)) {
    fit <- brm(f, data = df, family = Gamma(link = "log"), prior = priors_gamma,
               chains = 4, iter = 3000, warmup = 1000, control = ctrl, seed = SEED_MASTER)
  } else {
    fit <- brm(bf(f, hu ~ 1), data = df, family = hurdle_lognormal(), prior = priors_hurdle_logn,
               chains = 4, iter = 3000, warmup = 1000, control = ctrl, seed = SEED_MASTER)
  }
  saveRDS(fit, file); fit
}

pooled_files <- file.path(pooled_dir, paste0(safe_file(metric_cols), ".rds"))
if (!all(file.exists(pooled_files))) stop("Pooled models not found. Run 01_fish_summary_glmms.R first.", call. = FALSE)
pooled <- setNames(lapply(pooled_files, readRDS), metric_cols)
depth  <- setNames(lapply(metric_cols, fit_depth_model), metric_cols)

# ---- 2. Model comparison (leave-one-out cross-validation) ----------------------------
# elpd difference > 0 favours the depth model; a difference smaller than 2 SE is treated as no clear difference
comparison <- purrr::map_dfr(metric_cols, \(m) {
  l_pooled <- loo::loo(pooled[[m]]); l_depth <- loo::loo(depth[[m]])
  d  <- l_depth$pointwise[, "elpd_loo"] - l_pooled$pointwise[, "elpd_loo"]
  est <- sum(d); se <- sqrt(length(d) * var(d))
  tibble(Response = gsub("\\.", " ", m), `elpd difference (depth - pooled)` = round(est, 1), SE = round(se, 1),
         Conclusion = if (est > 2 * se) "depth model better" else if (est < -2 * se) "pooled model better" else "no clear difference")
})
save_table(comparison, file.path(res_dir, "Depth_model_comparison"),
           caption = "Leave-one-out comparison of models with and without Depth and Depth x Exposure.")

# ---- 3. Deep vs shallow rate ratios (deep / shallow) ---------------------------------
depth_ratios <- function(model, label) {
  one <- function(spec) {
    ct    <- contrast(emmeans(model, spec, type = "link", re_formula = NA), method = "revpairwise")
    draws <- as.data.frame(do.call(rbind, as.mcmc(ct)))
    lab   <- broom::tidy(ct) |> mutate(.col = colnames(draws))
    purrr::imap_dfr(draws, \(v, nm) tibble(
      .col = nm, RR = exp(mean(v)), Lower = exp(quantile(v, 0.025, names = FALSE)),
      Upper = exp(quantile(v, 0.975, names = FALSE)), `Pr(RR>1)` = mean(v > 0))) |>
      left_join(lab, by = ".col")
  }
  zones <- one(~ Depth | Exposure) |> transmute(Exposure = as.character(Exposure), RR, Lower, Upper, `Pr(RR>1)`)
  all_  <- one(~ Depth) |> transmute(Exposure = "All zones", RR, Lower, Upper, `Pr(RR>1)`)
  bind_rows(zones, all_) |> mutate(Response = label, .before = 1)
}
ratios <- purrr::imap_dfr(depth, \(m, nm) depth_ratios(m, gsub("\\.", " ", nm))) |>
  mutate(across(c(RR, Lower, Upper), \(x) round(x, 2)), `Pr(RR>1)` = round(`Pr(RR>1)`, 2),
         `RR [95% CrI]` = sprintf("%.2f [%.2f, %.2f]", RR, Lower, Upper)) |>
  select(Response, Exposure, `RR [95% CrI]`, `Pr(RR>1)`)
save_table(ratios, file.path(res_dir, "Depth_effects_by_exposure"),
           caption = "Deep / shallow rate ratio (RR) with 95% credible interval and Pr(RR>1), within each exposure zone and averaged over zones.")
message("Done: depth sensitivity tables written to ", file.path(OUT_DIR, res_dir))
