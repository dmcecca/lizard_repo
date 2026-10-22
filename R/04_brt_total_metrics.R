# =============================================================================
# 04_brt_total_metrics.R
# Boosted regression trees (BRT): relative influence of coral genera and exposure on total fish
# density, species richness and total biomass.
#
# Outputs (outputs/brt_total/):
#   Figure4_relative_influence.(pdf|png)                         Figure 4
#   Table_relative_influence_all_predictors_with_CI.csv          influence, bootstrap mean and 95% CI, in-sample R2
#
# Model:      gaussian BRT (dismo::gbm.step) on log(1 + response).
#             Predictors: 21 coral genera, generic coral richness (Coral.genera) and Exposure.
# Settings:   see helpers_brt.R (tree complexity 3, learning rate 0.002, bag fraction 0.75, up to 2000 trees;
#             number of trees chosen by 10-fold cross-validation).
# Run time:   minutes to an hour.
# =============================================================================

# Sets the working directory to the project folder (the folder holding R and data) if it isn't already.
# This works when the script is open in RStudio; otherwise run setwd("path to lizard_repo") first.
if (!file.exists("R/00_setup.R")) setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
source("R/00_setup.R")
source("R/helpers_brt.R")

res_dir <- "brt_total"
set.seed(123)

lizard <- read_lizard(); lizard$Exposure <- factor(lizard$Exposure)
predictors_all <- c(intersect(c(CORAL_PREDICTORS, "Coral.genera"), names(lizard)), "Exposure")
response_vars  <- intersect(c("Total.density", "Species.richness", "Total.biomass"), names(lizard))
stopifnot(length(response_vars) > 0, length(predictors_all) > 1)

compute_r2 <- function(actual, predicted) {
  ss_tot <- sum((actual - mean(actual, na.rm = TRUE))^2, na.rm = TRUE)
  if (is.na(ss_tot) || ss_tot == 0) return(NA_real_)
  1 - sum((actual - predicted)^2, na.rm = TRUE) / ss_tot
}

results <- purrr::map_dfr(response_vars, function(response) {
  message("Gaussian BRT for ", response)
  df <- na.omit(lizard[, c(response, predictors_all), drop = FALSE])
  if (nrow(df) < 10) return(NULL)
  df[[response]] <- log1p(as.numeric(df[[response]])); df$Exposure <- factor(df$Exposure)
  fit <- fit_brt(df, response, predictors_all, "gaussian")
  if (is.null(fit)) return(NULL)
  fitted_vals <- tryCatch(as.numeric(predict(fit, df, n.trees = fit$gbm.call$best.trees, type = "response")),
                          error = function(e) rep(NA_real_, nrow(df)))
  extract_rel_inf(fit, predictors_all) |>
    left_join(bootstrap_rel_inf(df, response, predictors_all, "gaussian"), by = "var") |>
    mutate(Response = response, R2 = compute_r2(df[[response]], fitted_vals))
})

tbl <- results |>
  transmute(Response = factor(Response, levels = c("Total.density", "Species.richness", "Total.biomass"),
                              labels = c("Total density", "Species richness", "Total biomass")),
            Predictor = var, rel.inf, boot_mean, lower_95 = lower, upper_95 = upper, n_success, R2) |>
  arrange(Response, desc(boot_mean))
write.csv(tbl, out_path(res_dir, "Table_relative_influence_all_predictors_with_CI.csv"), row.names = FALSE)

# ---- Figure 4 ----------------------------------------------------------------
# Order predictors within each panel by bootstrap mean (key = panel + predictor, labelled by predictor only)
plot_dat <- tbl |> filter(is.finite(boot_mean), is.finite(lower_95), is.finite(upper_95)) |>
  arrange(Response, boot_mean) |>
  mutate(Predictor_plot = factor(paste(Response, Predictor, sep = "___"), levels = unique(paste(Response, Predictor, sep = "___"))))
ref_dat <- plot_dat |> group_by(Response) |> summarise(ref_line = 100 / n(), .groups = "drop")  # equal-influence line

fig4 <- ggplot(plot_dat, aes(Predictor_plot, boot_mean)) +
  geom_segment(aes(xend = Predictor_plot, y = 0, yend = boot_mean), linewidth = 5, colour = "darkgrey") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.18, linewidth = 0.55) +
  geom_point(size = 1.9) +
  geom_hline(data = ref_dat, aes(yintercept = ref_line), linetype = "dashed", linewidth = 0.6) +
  coord_flip() + facet_wrap(~ Response, scales = "free_y", ncol = 3) + scale_x_discrete(labels = function(x) sub("^.*___", "", x)) +
  labs(x = NULL, y = "Bootstrap mean relative influence (%)") +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 10),
        panel.border = element_rect(linewidth = 0.6, fill = NA),
        axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8), legend.position = "none")
ggsave(out_path(res_dir, "Figure4_relative_influence.pdf"), fig4, width = 11, height = 6)
ggsave(out_path(res_dir, "Figure4_relative_influence.png"), fig4, width = 11, height = 6, dpi = 600)
message("Done: Figure 4 written to ", file.path(OUT_DIR, res_dir))
