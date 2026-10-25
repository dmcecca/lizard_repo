# Data

| File | Description | Used by |
|---|---|---|
| `lizard.fish.benthos.latest.csv` | One row per transect: Site, Year, Exposure, Depth, Transect, fish counts (one column per species, plus Total.density, Species.richness, Total.biomass and trophic-group biomass), and benthic cover (coral genera, Coral.genera, other benthic categories) | scripts 01–07, 09 |
| `fish coral relationship change table 2.csv` | Summary of how each fish–coral relationship changed through time, one row per coral genus and fish species | script 08 |

Column names are converted to R-safe names when the data are read (spaces become dots).

## fish_coral_pairs_figure6.csv

The 88 fish species × coral genus pairs used in the INLA exposure contrasts (script 06) and the temporal brms models (script 07, Figures 6 and S3). Columns: `coral`, `fish`. For each fish species whose BRT converged, the coral genus (or genera) with the highest relative influence was selected; pairs for which the temporal model converged in at least one exposure zone are kept.
