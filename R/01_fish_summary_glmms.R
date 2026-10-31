# =============================================================================
# 01_fish_summary_glmms.R
# Temporal trends in total fish density, species richness and total biomass in each exposure zone
# (Bayesian GLMMs, brms).
#
# Outputs (outputs/glmm_fish_metrics/):
#   EMM_fish_metrics.pdf                  Figure 2
#   Fish_GLMM_contrasts_with_EP           Table 2 (pairwise contrasts; .csv and .docx)
#   Fish_GLMM_results_with_EP             population-level effects (.csv and .docx)
#   Fish_model_info.csv                   family used for each response
#
# Model:    Response ~ Year * Exposure + (1 | Site) + (1 | SiteDepthTransect)
#           Gamma (log link) if all values are > 0, otherwise hurdle-lognormal.
#           4 chains x 3000 iterations (1000 warm-up).
# Run time: hours. Fitted models are saved and re-used on later runs.
# =============================================================================

# Sets the working directory to the project folder (the folder holding R and data) if it isn't already.
# This works when the script is open in RStudio; otherwise run setwd("path to lizard_repo") first.
if (!file.exists("R/00_setup.R")) setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
source("R/00_setup.R")
source("R/helpers_glmm.R")
require_pkgs(c("flextable", "officer"))   # only used to write .docx tables

set.seed(SEED_MASTER)
res_dir <- "glmm_fish_metrics"

# ---- 1. Data ----------------------------------------------------------------
dat         <- read_lizard()
df_base     <- prep_glmm_data(dat)
metric_cols <- col_range(dat, METRIC_FIRST, METRIC_LAST)   # Total.density, Species.richness, Total.biomass

# ---- 2. Fit (or load cached) models ------------------------------------------
model_dir <- dirname(out_path(res_dir, "brms_fish_models_site", ".keep"))
cached    <- file.path(model_dir, paste0(safe_file(metric_cols), ".rds"))
if (all(file.exists(cached))) {
  message("Loading cached models from ", model_dir)
  models <- setNames(lapply(cached, readRDS), metric_cols)
} else {
  models <- fit_all_glmms(metric_cols, df_base, model_dir,
                          info_file = out_path(res_dir, "Fish_model_info.csv"), id_col = "FishMetric")
}

# ---- 3. Tables ---------------------------------------------------------------
save_table(build_fixef_table(models), file.path(res_dir, "Fish_GLMM_results_with_EP"),
           caption = "Bayesian fish GLMM results with 95% CrI and exceedance probabilities.")
save_table(build_contrast_table(models), file.path(res_dir, "Fish_GLMM_contrasts_with_EP"),   # Table 2
           caption = "Pairwise contrasts by Year and Exposure; RR [95% CrI]; EP = Pr(RR > 1).")

# Condensed contrasts among main categories (Exposure averaged over Year; Year averaged over Exposure)
save_main_effect_table(build_main_effect_long(models), file.path(res_dir, "Fish_GLMM_main_effect_contrasts"),
  labels = c(Total.density = "Total density", Species.richness = "Species richness", Total.biomass = "Total biomass"),
  caption = MAIN_EFFECT_CAPTION)

# ---- 4. Figure 2: estimated marginal means -----------------------------------
panel_spec <- list(
  list(col = "Total.density",    title = "a) Total density",    ylab = "Total density"),
  list(col = "Species.richness", title = "b) Species richness", ylab = "Species richness"),
  list(col = "Total.biomass",    title = "c) Total biomass",    ylab = "Total biomass"))

panels <- lapply(panel_spec, function(p) {
  if (is.null(models[[p$col]])) stop("Model missing for ", p$col)
  plot_emm(get_emm_df(models[[p$col]], p$ylab), p$title, p$ylab)
})
legend_g <- cowplot::get_legend(panels[[1]] + theme(legend.position = "right"))
stack <- cowplot::plot_grid(plotlist = lapply(panels, \(p) p + theme(legend.position = "none")),
                            ncol = 1, align = "v", axis = "lr")
fig2 <- cowplot::plot_grid(stack, legend_g, ncol = 2, rel_widths = c(0.8, 0.2))
ggsave(out_path(res_dir, "EMM_fish_metrics.pdf"), fig2, width = 8.2, height = 10.5, dpi = 600)
message("Done: Figure 2 and Table 2 written to ", file.path(OUT_DIR, res_dir))
