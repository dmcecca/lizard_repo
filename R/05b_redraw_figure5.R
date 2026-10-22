# =============================================================================
# 05b_redraw_figure5.R
# Redraws Figure 5 from the saved results of 05_brt_species.R (no model refitting; takes seconds).
# Change the sizes below to adjust the figure.
#
# Reads:     outputs/brt_species/BRT_Relative_Influence_Results_selected_corals.csv
# Produces:  outputs/brt_species/Figure5_tall, Figure5_part1, Figure5_part2 (.pdf, .png)
# =============================================================================

source(file.path("R", "00_setup.R"))
source(file.path("R", "helpers_figure5.R"))
require_pkgs(c("ggplot2", "dplyr"))

ROW_HEIGHT   <- 0.20   # inches per species row (larger = more spaced out)
SPECIES_TEXT <- 9      # font size of the species names
WIDTH        <- 11     # inches

results <- read.csv(out_path("brt_species", "BRT_Relative_Influence_Results_selected_corals.csv"), stringsAsFactors = FALSE)
draw_figure5(results, file.path(OUT_DIR, "brt_species"), row_height = ROW_HEIGHT, species_text = SPECIES_TEXT, width = WIDTH)
message("Done: figures written to ", file.path(OUT_DIR, "brt_species"))
