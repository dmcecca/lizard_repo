# Fish–coral relationships at Lizard Island

R code for the analyses and figures in the manuscript on relationships between fish assemblages and coral genera through a disturbance–recovery cycle at Lizard Island (Great Barrier Reef).

## How to run

1. Install R and RStudio, and download this repository (green **Code** button → Download ZIP, or clone it with GitHub Desktop).
2. Put the data file in the `data` folder (see `data/README.md`).
3. Open a script from the `R` folder in RStudio and click **Source**. Each script sets the working directory to the project folder itself, so no paths need editing. (Outside RStudio, run `setwd("path to lizard_repo")` first.) Outputs are written to a new `outputs` folder.

Scripts can be run in any order, except:

- 05b needs 05
- 10 needs 01 and 02
- 11 needs 01
- 08 and 08b need 07

## Scripts

| Script | What it does | Output |
|---|---|---|
| `01_fish_summary_glmms.R` | Bayesian GLMMs of total fish density, richness and biomass, with random intercepts for site and for transect position | Figure 2, Table 2 |
| `02_trophic_glmms.R` | Bayesian GLMMs of trophic-group biomass (same random-effect structure as 01) | Figure S2, Table S2 |
| `03_community_composition.R` | PERMANOVA, dispersion tests and SIMPER | Figure 3 (nMDS), Table 3, Table S3 |
| `04_brt_total_metrics.R` | BRTs for total density, richness and biomass | Figure 4 |
| `05_brt_species.R` | BRTs for each fish species | Figure 5, Table S5 |
| `05b_redraw_figure5.R` | Redraws Figure 5 from saved results (seconds) | Figure 5 |
| `06_exposure_slope_contrasts.R` | Differences between exposure zones in the slopes of the fish–coral pairs used in the temporal analysis (INLA) | contrast tables |
| `06b_exposure_contrast_figure.R` | Figure and table summarising the exposure contrasts from 06 | Figure S3, Table S6 |
| `07_temporal_slopes_brms.R` | Bayesian (brms) fish–coral slopes in 2011, 2015/2017 and 2020 for each pair and exposure zone (zero-inflated Poisson on density values); classifies each pair × zone as stable, stronger, weaker, change or bounce back; builds the Figure 6 table | relationship-change tables |
| `07b_temporal_slopes_brms_counts.R` | Check on 07: the same slopes from raw counts with a negative binomial model | relationship-change table (counts) |
| `08_figure6_relationship_change.R` | Draws Figure 6 from the relationship-change table | Figure 6 |
| `08b_figureS4_caterpillar.R` | Draws Figure S4a–d (slopes in each period, one plot per exposure zone) from the brms output of 07 | Figure S4a–d, pattern tables |
| `09_brt_out_of_sample.R` | Out-of-sample accuracy of the BRTs | Table S4 |
| `10_condensed_glmm_tables.R` | Condensed contrast tables from 01 and 02 | Table 2, Table S2 |
| `11_depth_sensitivity.R` | Tests whether pooling shallow and deep transects matters | depth comparison tables |
| `13_benthic_covariation.R` | Correlation matrix of coral genera and other benthic categories | Figure S1 and tables |

`00_setup.R` and the `helpers_*.R` files hold settings and functions used by the scripts above.

## Run times

The Bayesian models (01, 02, 07), the INLA models (06) and the species BRTs (05, 09) take hours. Fitted models and per-species results are saved in the `outputs` folder, so an interrupted run restarts where it stopped. These saved model files are large and are not included in the repository; they are recreated when the scripts are run.

## R packages

brms, emmeans, INLA, gbm, dismo, vegan, permute, dplyr, tidyr, purrr, tibble, stringr, ggplot2, ggrepel, patchwork, cowplot, readr, flextable, officer, coda, broom. Each script checks for the packages it needs and says which are missing.
