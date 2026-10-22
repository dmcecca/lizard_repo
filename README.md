# Fish–coral relationships at Lizard Island

R code for the analyses and figures in the manuscript on relationships between fish assemblages and coral genera through a disturbance–recovery cycle at Lizard Island (Great Barrier Reef).

## How to run

1. Install R and RStudio, and download this repository (green **Code** button → Download ZIP, or clone it with GitHub Desktop).
2. Put the data file in the `data` folder (see `data/README.md`).
3. Double-click `lizard-fish-coral.Rproj` to open the project in RStudio. This sets the working directory, so no paths need editing.
4. Open a script from the `R` folder and click **Source**. Outputs are written to a new `outputs` folder.

Scripts can be run in any order, except 05b (needs 05) and 10 (needs 01 and 02).

## Scripts

| Script | What it does | Output |
|---|---|---|
| `01_fish_summary_glmms.R` | Bayesian GLMMs of total density, richness and biomass | Figure 2, Table 2 |
| `02_trophic_glmms.R` | Bayesian GLMMs of trophic groups | Figure S1, Table S2 (trophic) |
| `03_community_composition.R` | PERMANOVA, dispersion tests and SIMPER | Figure 3 (nMDS), Table 3, Table S3 |
| `04_brt_total_metrics.R` | BRTs for total density, richness and biomass | Figure 4 |
| `05_brt_species.R` | BRTs for each fish species | Figure 5, Table S5 |
| `05b_redraw_figure5.R` | Redraws Figure 5 from saved results (seconds) | Figure 5 |
| `06_exposure_slope_contrasts.R` | Differences between exposure zones in species–coral slopes | contrast tables, Figure S2 |
| `07_temporal_slope_template.R` | Template for year-specific slopes of one species–coral pair | slope tables |
| `08_figure6_relationship_change.R` | Draws Figure 6 from the relationship-change table | Figure 6 |
| `09_brt_out_of_sample.R` | Out-of-sample accuracy of the BRTs | validation tables, Table S4 |
| `10_condensed_glmm_tables.R` | Condensed contrast tables | Table S2 |

`00_setup.R` and the `helpers_*.R` files hold settings and functions used by the scripts above.

## Run times

The Bayesian models (01, 02) and species BRTs (05, 09) take hours. Fitted models and per-species results are saved, so an interrupted run restarts where it stopped.

## R packages

brms, emmeans, INLA, gbm, dismo, vegan, dplyr, tidyr, purrr, tibble, stringr, ggplot2, cowplot, readr, flextable, officer, coda, broom. Each script checks for the packages it needs and says which are missing.
