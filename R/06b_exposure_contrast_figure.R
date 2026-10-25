# =============================================================================
# 06b_exposure_contrast_figure.R
# Figure and tables summarising script 06: do fish-coral slopes differ among exposure zones?
# Uses the contrasts for the fish x coral pairs identified as influential in the BRTs (one contrast = one fish species x coral genus x pair of exposure zones;
# a contrast is "supported" when its 95% credible interval excludes zero).
#
# Input:   outputs/inla_exposure_slopes/all_exposure_slope_contrasts.csv   (from 06)
# Outputs (outputs/exposure_contrast_summary/):
#   Exposure_contrasts_summary.pdf / .png     Figure: (a) share of contrasts supported for each pair of zones,
#                                             (b) share supported for each coral genus, (c) size of the slope differences
#   Table_by_zone_contrast.csv                per pair of zones: contrasts, supported, % supported, median size of supported differences
#   Table_by_coral_genus.csv                  per coral genus: contrasts, supported, % supported
#   Exposure_contrast_summary_numbers.csv     the numbers quoted in the text
# Slope differences are on the log scale per SD of coral cover (coral cover is standardised within each model).
# =============================================================================

# Sets the working directory to the project folder (the folder holding R and data) if it isn't already.
# This works when the script is open in RStudio; otherwise run setwd("path to lizard_repo") first.
if (!file.exists("R/00_setup.R")) setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
source("R/00_setup.R")
require_pkgs(c("patchwork")); library(patchwork)

res_dir <- "exposure_contrast_summary"
in_file <- out_path("inla_exposure_slopes", "all_exposure_slope_contrasts.csv")
if (!file.exists(in_file)) stop("Run R/06_exposure_slope_contrasts.R first (it makes ", in_file, ").", call. = FALSE)

pair_list <- brt_pairs()
con <- read.csv(in_file, stringsAsFactors = FALSE) |>
  semi_join(pair_list, by = c("fish", "coral")) |>
  mutate(supported = q025 > 0 | q975 < 0,
         zones = gsub("_minus_", " vs ", contrast),
         coral = gsub("\\.", " ", coral))

# ---- Tables ---------------------------------------------------------------------------
by_zone <- con |> group_by(zones) |>
  summarise(contrasts = n(), pairs = n_distinct(fish, coral), n_supported = sum(supported),
            percent_supported = round(100 * n_supported / contrasts, 1),
            median_abs_difference_supported = round(median(abs(mean[supported])), 2), .groups = "drop") |>
  arrange(desc(percent_supported))
by_coral <- con |> group_by(coral) |>
  summarise(contrasts = n(), species = n_distinct(fish), n_supported = sum(supported),
            percent_supported = round(100 * n_supported / contrasts, 1), .groups = "drop") |>
  arrange(desc(percent_supported))
write.csv(by_zone,  out_path(res_dir, "Table_by_zone_contrast.csv"), row.names = FALSE)
write.csv(by_coral, out_path(res_dir, "Table_by_coral_genus.csv"), row.names = FALSE)

pairs_all <- con |> distinct(fish, coral)
numbers <- tibble(
  measure = c("fish x coral pairs modelled", "fish species", "coral genera", "contrasts", "contrasts supported (95% CrI excludes zero)",
              "percent supported", "pairs with at least one supported contrast", "percent of pairs",
              "median |difference|, all contrasts", "median |difference|, supported contrasts"),
  value = c(nrow(pairs_all), n_distinct(con$fish), n_distinct(con$coral), nrow(con), sum(con$supported),
            round(100 * mean(con$supported), 1),
            n_distinct(con$fish[con$supported], con$coral[con$supported]),
            round(100 * n_distinct(con$fish[con$supported], con$coral[con$supported]) / nrow(pairs_all), 1),
            round(median(abs(con$mean)), 2), round(median(abs(con$mean[con$supported])), 2)))
write.csv(numbers, out_path(res_dir, "Exposure_contrast_summary_numbers.csv"), row.names = FALSE)
print(numbers, n = Inf); print(by_zone)

# ---- Figure --------------------------------------------------------------------------------
base <- theme_classic(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 9))
pa <- ggplot(by_zone, aes(percent_supported, reorder(zones, percent_supported))) +
  geom_col(fill = "grey35", width = 0.7) +
  geom_text(aes(label = paste0(n_supported, " / ", contrasts)), hjust = -0.1, size = 2.8) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(title = "a  Exposure zones compared", x = "Contrasts supported (%)", y = NULL) + base
pb <- ggplot(by_coral, aes(percent_supported, reorder(coral, percent_supported))) +
  geom_col(fill = "grey35", width = 0.7) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "b  Coral genus", x = "Contrasts supported (%)", y = NULL) + base +
  theme(axis.text.y = element_text(face = "italic", size = 7))
lim <- as.numeric(quantile(abs(con$mean), 0.99))
pc <- con |> mutate(group = ifelse(supported, "95% CrI excludes zero", "95% CrI includes zero")) |>
  filter(abs(mean) <= lim) |>
  ggplot(aes(mean, fill = group)) +
  geom_histogram(bins = 80, position = "identity", alpha = 0.8) +
  scale_fill_manual(values = c("95% CrI excludes zero" = "#c0392b", "95% CrI includes zero" = "grey65"), name = NULL) +
  scale_y_sqrt() +
  labs(title = "c  Size of the slope differences", x = "Slope difference (log scale per SD of coral cover)", y = "Contrasts (square-root scale)") +
  base + theme(legend.position = "bottom")
fig <- (pa / pc) | pb
fig <- fig + plot_layout(widths = c(1.2, 1))
ggsave(out_path(res_dir, "Exposure_contrasts_summary.pdf"), fig, width = 9, height = 6)
ggsave(out_path(res_dir, "Exposure_contrasts_summary.png"), fig, width = 9, height = 6, dpi = 300)
message("Done. Results are in outputs/", res_dir, "/")
