# =============================================================================
# 02_trophic_glmms.R
# Temporal trends in the biomass of each trophic group in each exposure zone (same GLMM as script 01).
#
# Outputs (outputs/glmm_trophic/):
#   EMM_trophic_groups_ALL_3col.pdf       Figure S1
#   Trophic_GLMM_contrasts_with_EP        Table S2 (pairwise contrasts; .csv and .docx)
#   Trophic_GLMM_results_with_EP          population-level effects (.csv and .docx)
#   Trophic_model_info.csv                family used for each group
#
# Run time: hours. Fitted models are saved and re-used on later runs.
# =============================================================================

source("R/00_setup.R")
source("R/helpers_glmm.R")
require_pkgs(c("flextable", "officer"))

set.seed(SEED_MASTER)
res_dir <- "glmm_trophic"

# ---- 1. Data ----------------------------------------------------------------
dat          <- read_lizard()
df_base      <- prep_glmm_data(dat)
trophic_cols <- col_range(dat, TROPHIC_FIRST, TROPHIC_LAST)   # 18 columns, B_Benthic.invertivore ... B_Corallivore
message(length(trophic_cols), " trophic-group columns: ", paste(trophic_cols, collapse = ", "))

# ---- 2. Fit (or load cached) models ------------------------------------------
model_dir <- dirname(out_path(res_dir, "brms_trophic_models", ".keep"))
cached    <- file.path(model_dir, paste0(safe_file(trophic_cols), ".rds"))
if (all(file.exists(cached))) {
  message("Loading cached models from ", model_dir)
  models <- setNames(lapply(cached, readRDS), trophic_cols)
} else {
  models <- fit_all_glmms(trophic_cols, df_base, model_dir,
                          info_file = out_path(res_dir, "Trophic_model_info.csv"), id_col = "TrophicGroup")
}

# ---- 3. Tables ---------------------------------------------------------------
save_table(build_fixef_table(models), file.path(res_dir, "Trophic_GLMM_results_with_EP"),
           caption = "Bayesian trophic-group GLMM results with 95% CrI and exceedance probabilities.")
save_table(build_contrast_table(models), file.path(res_dir, "Trophic_GLMM_contrasts_with_EP"),   # Table S2
           caption = "Pairwise contrasts by Year and Exposure; RR [95% CrI]; EP = Pr(RR > 1).")

# ---- 3b. Condensed contrasts among main categories (same layout as the fish-metric table) ----------
# Column order and labels follow Figure S1 (defined in section 4 below, so build them here first).
fig_s1_labels <- c(B_Benthic.invertivore = "Benthic invertivore", B_Carnivore = "Carnivore", B_Omnivore = "Omnivore",
  B_Piscivore = "Piscivore", B_Planktivore = "Planktivore", B_Grazer = "Grazer - All", B_Detritivore = "Grazer - Detritivore",
  B_Cropper = "Grazer - Cropper", B_Browser = "Grazer - Browser", B_Scraper = "Grazer - Scraper", B_Excavator = "Excavator",
  B_Farmer = "Farmer - All", B_Intensive.farmer = "Intensive farmer", B_Extensive.farmer = "Extensive farmer",
  B_Indeterminate.farmer = "Indeterminate farmer", B_Corallivore = "Corallivore - All",
  B_Obligate.corallivore = "Obligate corallivore", B_Facultative.corallivore = "Facultative corallivore")
save_main_effect_table(build_main_effect_long(models), file.path(res_dir, "Trophic_GLMM_main_effect_contrasts"),
  labels = fig_s1_labels[names(fig_s1_labels) %in% names(models)], caption = MAIN_EFFECT_CAPTION)

# ---- 4. Figure S1: panel order and labels --------------------------------------
# Published panel order (a-r), keyed on the data column name without the "B_" prefix.
fig_s1 <- tibble::tribble(
  ~column,                       ~label,
  "B_Benthic.invertivore",       "Benthic invertivore",     # a
  "B_Carnivore",                 "Carnivore",               # b
  "B_Omnivore",                  "Omnivore",                # c
  "B_Piscivore",                 "Piscivore",               # d
  "B_Planktivore",               "Planktivore",             # e
  "B_Grazer",                    "Grazer - All",            # f
  "B_Detritivore",               "Grazer - Detritivore",    # g
  "B_Cropper",                   "Grazer - Cropper",        # h
  "B_Browser",                   "Grazer - Browser",        # i
  "B_Scraper",                   "Grazer - Scraper",        # j
  "B_Excavator",                 "Excavator",               # k
  "B_Farmer",                    "Farmer - All",            # l
  "B_Intensive.farmer",          "Intensive farmer",        # m
  "B_Extensive.farmer",          "Extensive farmer",        # n
  "B_Indeterminate.farmer",      "Indeterminate farmer",    # o
  "B_Corallivore",               "Corallivore - All",       # p
  "B_Obligate.corallivore",      "Obligate corallivore",    # q
  "B_Facultative.corallivore",   "Facultative corallivore") # r

missing <- setdiff(fig_s1$column, names(models))
if (length(missing)) stop("Trophic columns not found for Fig S1: ", paste(missing, collapse = ", "))
extra <- setdiff(names(models), fig_s1$column)
if (length(extra)) message("Models not in published Fig S1 (skipped in the figure): ", paste(extra, collapse = ", "))

# ---- 5. Figure S1: estimated marginal means ----------------------------------
panels <- purrr::map2(fig_s1$column, seq_len(nrow(fig_s1)), function(col, i) {
  plot_emm(get_emm_df(models[[col]], fig_s1$label[i]),
           panel_title = paste0(letters[i], ") ", fig_s1$label[i]),
           y_lab = "Estimated mean biomass", compact = TRUE)
})
names(panels) <- fig_s1$column

purrr::iwalk(panels, ~ ggsave(out_path(res_dir, "panels", paste0("EMM_", safe_file(.y), ".png")),
                              .x, width = 5.2, height = 4.1, dpi = 300))

legend_g   <- cowplot::get_legend(panels[[1]] + theme(legend.position = "right"))
panel_grid <- cowplot::plot_grid(plotlist = lapply(panels, \(p) p + theme(legend.position = "none")),
                                 ncol = 3, align = "hv")
fig <- cowplot::plot_grid(panel_grid, legend_g, ncol = 2, rel_widths = c(0.88, 0.12))
ggsave(out_path(res_dir, "EMM_trophic_groups_ALL_3col.pdf"), fig,
       width = 12, height = 3.1 * ceiling(length(panels) / 3), dpi = 600)
message("Done: Figure S1 and Table S2 written to ", file.path(OUT_DIR, res_dir))
