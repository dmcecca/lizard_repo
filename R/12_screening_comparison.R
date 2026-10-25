# =============================================================================
# 12_screening_comparison.R
# How robust are the fish-coral exposure contrasts to the data thresholds used to select models?
# Compares an earlier run of script 06 (lenient thresholds) with a new run (stricter thresholds), for the fish x coral pairs
# in data/fish_coral_pairs_figure6.csv:
#   how many significant contrasts came from species recorded in fewer than 10 transects,
#   and how many contrasts from the earlier run are still found, with the same direction, in the new run.
#
# Before running: keep the earlier outputs in  outputs/inla_exposure_slopes_old/
#                 and the new outputs in       outputs/inla_exposure_slopes/
# Produces:  outputs/screening_comparison/screening_comparison_summary.csv
# =============================================================================

if (!file.exists("R/00_setup.R")) setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
source("R/00_setup.R")

old_file <- file.path(OUT_DIR, "inla_exposure_slopes_old", "all_exposure_slope_contrasts.csv")
new_file <- file.path(OUT_DIR, "inla_exposure_slopes",     "all_exposure_slope_contrasts.csv")
if (!all(file.exists(c(old_file, new_file)))) stop("Both all_exposure_slope_contrasts.csv files are needed; see the header.", call. = FALSE)

brt_pairs <- read.csv(file.path("data", "fish_coral_pairs_figure6.csv"), stringsAsFactors = FALSE)
flag <- function(x) mutate(x, significant = q025 > 0 | q975 < 0, direction = sign(mean))
old <- flag(semi_join(read.csv(old_file), brt_pairs, by = c("fish", "coral")))
new <- flag(semi_join(read.csv(new_file), brt_pairs, by = c("fish", "coral")))

key <- c("fish", "coral", "contrast")
both <- inner_join(select(old, all_of(key), old_sig = significant, old_dir = direction),
                   select(new, all_of(key), new_sig = significant, new_dir = direction), by = key)
old_sig <- both |> filter(old_sig)

summary_tab <- tibble(
  Measure = c("Contrasts in earlier run", "Significant contrasts in earlier run",
              "  of which from species recorded in fewer than 10 transects", "Contrasts in new run",
              "Significant contrasts in new run", "Earlier significant contrasts also present in new run",
              "  of which still significant, same direction", "  of which same direction (significant or not)"),
  Value = c(nrow(old), sum(old$significant), sum(old$significant & old$n_nonzero < 10), nrow(new), sum(new$significant),
            nrow(old_sig), sum(old_sig$new_sig & old_sig$new_dir == old_sig$old_dir), sum(old_sig$new_dir == old_sig$old_dir)))
print(as.data.frame(summary_tab))
write.csv(summary_tab, out_path("screening_comparison", "screening_comparison_summary.csv"), row.names = FALSE)
