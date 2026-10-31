# =============================================================================
# 12_site_effect_sensitivity.R
# Does adding a site-level random intercept change the conclusions for total density, richness and biomass?
# Refits the models from script 01 with transect positions nested within sites, and compares them:
#   pooled model (script 01):  Response ~ Year * Exposure + (1 | SiteDepthTransect)
#   site model:                Response ~ Year * Exposure + (1 | Site) + (1 | SiteDepthTransect)
# The site term models the correlation among different transects at the same site, which the pooled model does not.
#
# Produces (outputs/site_effect_sensitivity/):
#   Site_model_comparison     LOO comparison (elpd_diff of the site model minus the pooled model)
#   Site_variance             SD of the site and transect-position random intercepts in the site model
#   Contrast_agreement        for the contrasts in Table 2 (years within each exposure zone, exposure zones within each year):
#                             how many have strong evidence (exceedance probability > 0.95) in each model, and how many differ
#   Contrast_comparison       every contrast with its exceedance probability in both models
# Run time: as long as script 01 (hours). The site models are saved and re-used on later runs.
# Run script 01 first, so that the pooled models exist.
# =============================================================================

# Sets the working directory to the project folder (the folder holding R and data) if it isn't already.
# This works when the script is open in RStudio; otherwise run setwd("path to lizard_repo") first.
if (!file.exists("R/00_setup.R")) setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
source("R/00_setup.R")
source(file.path("R", "helpers_glmm.R"))
require_pkgs(c("flextable", "officer", "loo"))

res_dir <- "site_effect_sensitivity"
set.seed(SEED_MASTER)

dat         <- read_lizard()
df_base     <- prep_glmm_data(dat)
metric_cols <- col_range(dat, METRIC_FIRST, METRIC_LAST)

pooled_dir <- file.path(OUT_DIR, "glmm_fish_metrics", "brms_fish_models")
site_dir   <- dirname(out_path(res_dir, "brms_site_models", ".keep"))

# ---- 1. Fit (or load) the site models, using the same family as the pooled model ----
fit_site_model <- function(resp_col) {
  file <- file.path(site_dir, paste0(safe_file(resp_col), ".rds"))
  if (file.exists(file)) return(readRDS(file))
  df <- df_base |>
    transmute(Year, Exposure, Site = factor(Site), SiteDepthTransect, Response = as.numeric(.data[[resp_col]])) |>
    filter(!is.na(Response))
  f <- Response ~ Year * Exposure + (1 | Site) + (1 | SiteDepthTransect)
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
site   <- setNames(lapply(metric_cols, fit_site_model), metric_cols)

# ---- 2. Model comparison (leave-one-out cross-validation) ----------------------------
# elpd difference > 0 favours the site model; a difference smaller than 2 SE is treated as no clear difference
comparison <- purrr::map_dfr(metric_cols, \(m) {
  l_pooled <- loo::loo(pooled[[m]]); l_site <- loo::loo(site[[m]])
  d   <- l_site$pointwise[, "elpd_loo"] - l_pooled$pointwise[, "elpd_loo"]
  est <- sum(d); se <- sqrt(length(d) * var(d))
  tibble(Response = gsub("\\.", " ", m), `elpd difference (site - pooled)` = round(est, 1), SE = round(se, 1),
         Conclusion = if (est > 2 * se) "site model better" else if (est < -2 * se) "pooled model better" else "no clear difference")
})
save_table(comparison, file.path(res_dir, "Site_model_comparison"),
           caption = "Leave-one-out comparison of models with and without a site random intercept.")

# ---- 3. Size of the site and transect-position random intercepts ----------------------
variance <- purrr::imap_dfr(site, \(m, nm) {
  s <- brms::posterior_summary(m, variable = c("sd_Site__Intercept", "sd_SiteDepthTransect__Intercept"))
  tibble(Response = gsub("\\.", " ", nm),
         `SD site [95% CrI]` = sprintf("%.2f [%.2f, %.2f]", s[1, "Estimate"], s[1, "Q2.5"], s[1, "Q97.5"]),
         `SD transect position [95% CrI]` = sprintf("%.2f [%.2f, %.2f]", s[2, "Estimate"], s[2, "Q2.5"], s[2, "Q97.5"]))
})
save_table(variance, file.path(res_dir, "Site_variance"),
           caption = "Standard deviation of the random intercepts (log scale) in the models with a site term.")

# ---- 4. Contrasts in both models -----------------------------------------------------
# Exceedance probability = the larger of Pr(contrast > 0) and Pr(contrast < 0), as in Table 2
contrast_table <- function(model, spec, family) {
  ct    <- contrast(emmeans(model, spec, type = "link", re_formula = NA), method = "pairwise")
  draws <- as.data.frame(do.call(rbind, as.mcmc(ct)))
  lab   <- broom::tidy(ct) |> mutate(.col = colnames(draws))
  purrr::imap_dfr(draws, \(v, nm) tibble(.col = nm, p_gt0 = mean(v > 0))) |>
    left_join(lab, by = ".col") |>
    mutate(Family = family, EP = pmax(p_gt0, 1 - p_gt0))
}
all_contrasts <- function(model) {
  bind_rows(contrast_table(model, ~ Year | Exposure, "Years within each exposure zone"),
            contrast_table(model, ~ Exposure | Year, "Exposure zones within each year"))
}
evidence <- function(ep) cut(ep, c(-Inf, 0.9, 0.95, Inf), labels = c("none", "possible", "strong"), right = FALSE)

comp <- purrr::imap_dfr(metric_cols, \(m, i) {
  a <- all_contrasts(pooled[[m]]); b <- all_contrasts(site[[m]])
  by_cols <- intersect(c("Exposure", "Year", "contrast", "Family"), intersect(names(a), names(b)))
  inner_join(select(a, all_of(by_cols), EP_pooled = EP), select(b, all_of(by_cols), EP_site = EP), by = by_cols) |>
    mutate(Response = gsub("\\.", " ", m), .before = 1,
           Evidence_pooled = evidence(EP_pooled), Evidence_site = evidence(EP_site),
           Changed = Evidence_pooled != Evidence_site)
})
write.csv(mutate(comp, across(c(EP_pooled, EP_site), \(x) round(x, 3))),
          out_path(res_dir, "Contrast_comparison.csv"), row.names = FALSE)

agreement <- comp |> group_by(Response, Family) |>
  summarise(Contrasts = n(),
            `Strong evidence, pooled model` = sum(Evidence_pooled == "strong"),
            `Strong evidence, site model`   = sum(Evidence_site == "strong"),
            `Contrasts whose evidence category differs` = sum(Changed),
            `Strong in pooled but not in site model` = sum(Evidence_pooled == "strong" & Evidence_site != "strong"),
            .groups = "drop")
save_table(agreement, file.path(res_dir, "Contrast_agreement"),
           caption = "Number of contrasts with strong evidence (exceedance probability > 0.95) with and without a site random intercept.")
print(agreement, n = Inf, width = Inf)
message("Done: site-effect sensitivity tables written to ", file.path(OUT_DIR, res_dir))
