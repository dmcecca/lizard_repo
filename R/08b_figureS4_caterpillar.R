# =============================================================================
# 08b_figureS4_caterpillar.R
# Figure S4: fish-coral slopes in 2011, 2015/2017 and 2020, one caterpillar plot per exposure zone,
# and the pattern of each fish-coral pair (e.g. Neutral, Positive-neutral, Negative-positive-negative).
# Uses the slopes estimated by script 07, so Figure S4 and Figure 6 come from the same models and the same pairs.
#
# Input:   outputs/temporal_slopes/period_slopes_long.csv   (made by 07_temporal_slopes_brms.R)
# Outputs (outputs/figureS4/):
#   Figure_S4a_Fish-coral_slope_comparison_exposed.pdf, S4b oblique, S4c lagoon, S4d sheltered
#   pattern_by_pair.csv             one row per pair and zone: the three symbols, the pattern, stable or not
#   pattern_summary_by_zone.csv     number and percentage of pairs in each pattern, per zone
#   stable_by_zone.csv              per zone: pairs with the same symbol in all three periods, and pairs in the Stable class of Figure 6
#   comparison_with_hand_compiled.csv   pattern of each pair against the original hand-compiled figure (if data/figure_S3_hand_compiled_slopes.csv is present)
# Patterns: each period's slope is + or - if its 95% interval excludes zero, otherwise n; repeated neighbours are merged
#           (e.g. "+ n n" = Positive-neutral). Stable = Neutral, Positive or Negative in all three periods.
# =============================================================================

# Sets the working directory to the project folder (the folder holding R and data) if it isn't already.
# This works when the script is open in RStudio; otherwise run setwd("path to lizard_repo") first.
if (!file.exists("R/00_setup.R")) setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
source("R/00_setup.R")
library(ggplot2)

res_dir   <- "figureS4"
in_file   <- out_path("temporal_slopes", "period_slopes_long.csv")
if (!file.exists(in_file)) stop("Run R/07_temporal_slopes_brms.R first (it makes ", in_file, ").", call. = FALSE)

EXPOSURES <- c("Exposed", "Oblique", "Lagoon", "Sheltered")     # order of Figures S4a-d
PERIODS   <- c("2011", "2015/2017", "2020")

long <- read.csv(in_file, stringsAsFactors = FALSE) |>
  mutate(zone = factor(zone, EXPOSURES), period = factor(period, PERIODS),
         symbol = ifelse(lower95 > 0, "+", ifelse(upper95 < 0, "-", "n")))

# ---- Pattern of each pair ---------------------------------------------------------
pattern_name <- function(s) {
  s <- s[c(TRUE, s[-1] != s[-length(s)])]
  x <- tolower(paste(c("+" = "Positive", "-" = "Negative", "n" = "Neutral")[s], collapse = "-"))
  paste0(toupper(substring(x, 1, 1)), substring(x, 2))
}
by_pair <- long |> arrange(zone, coral, fish, period) |> group_by(zone, coral, fish) |>
  summarise(symbols = paste(symbol, collapse = " "), pattern = pattern_name(symbol), .groups = "drop") |>
  mutate(stable = pattern %in% c("Neutral", "Positive", "Negative"))
# Trajectory class from script 07 (the table behind Figure 6), which also takes the size of the slope into account
if ("trajectory" %in% names(long))
  by_pair <- by_pair |> left_join(distinct(long, zone, coral, fish, trajectory), by = c("zone", "coral", "fish"))
write.csv(by_pair, out_path(res_dir, "pattern_by_pair.csv"), row.names = FALSE)

pattern_summary <- by_pair |> count(zone, pattern) |> group_by(zone) |>
  mutate(percent = round(100 * n / sum(n), 1)) |> ungroup()
write.csv(pattern_summary, out_path(res_dir, "pattern_summary_by_zone.csv"), row.names = FALSE)

stable_by_zone <- by_pair |> group_by(zone) |>
  summarise(pairs = n(), n_same_sign_all_periods = sum(stable), percent_same_sign = round(100 * n_same_sign_all_periods / pairs, 1),
            n_stable_class = if ("trajectory" %in% names(by_pair)) sum(trajectory == "Stable") else NA_integer_,
            percent_stable_class = round(100 * n_stable_class / pairs, 1))
write.csv(stable_by_zone, out_path(res_dir, "stable_by_zone.csv"), row.names = FALSE)
print(stable_by_zone)

hand_file <- file.path("data", "figure_S3_hand_compiled_slopes.csv")
if (file.exists(hand_file)) {
  hand <- read.csv(hand_file, stringsAsFactors = FALSE) |> distinct(zone, coral, fish, hand_class)
  cmp <- inner_join(by_pair, hand, by = c("zone", "coral", "fish")) |> mutate(same_pattern = pattern == hand_class)
  write.csv(cmp, out_path(res_dir, "comparison_with_hand_compiled.csv"), row.names = FALSE)
  message(nrow(cmp), " pairs are in both; ", round(100 * mean(cmp$same_pattern), 1), "% have the same pattern")
}

# ---- Figures S4a-d -------------------------------------------------------------------
fig_names <- c(Exposed = "Figure_S4a_Fish-coral_slope_comparison_exposed.pdf",
               Oblique = "Figure_S4b_Fish-coral_slope_comparison_oblique.pdf",
               Lagoon  = "Figure_S4c_Fish-coral_slope_comparison_lagoon.pdf",
               Sheltered = "Figure_S4d_Fish-coral_slope_comparison_sheltered.pdf")
for (z in EXPOSURES) {
  d <- long |> filter(zone == z)
  if (!nrow(d)) next
  d <- d |> mutate(fish_lab = gsub("\\.", " ", fish), coral = factor(coral, levels = rev(sort(unique(coral)))))
  p <- ggplot(d, aes(x = slope, y = fish_lab)) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0, linewidth = 0.3) +
    geom_point(size = 0.8) +
    facet_grid(coral ~ period, scales = "free_y", space = "free_y", switch = "y") +
    labs(x = "Slope (log scale, per standard deviation of coral cover)", y = "Fish") +
    theme_classic(base_size = 8) +
    theme(axis.text.y = element_text(face = "italic", size = 5),
          strip.text.y.left = element_text(angle = 0, face = "italic", hjust = 1),
          strip.background = element_blank(), strip.placement = "outside",
          panel.spacing.y = unit(0.15, "lines"))
  n_pairs <- n_distinct(d$coral, d$fish)
  ggsave(out_path(res_dir, fig_names[[z]]), p, width = 8,
         height = max(5, 1.6 + 0.09 * n_pairs + 0.15 * n_distinct(d$coral)), limitsize = FALSE)
}
message("Done. Results are in outputs/", res_dir, "/")

