# Data

| File | Description | Used by |
|---|---|---|
| `lizard.fish.benthos.latest.csv` | One row per transect: Site, Year, Exposure, Depth, Transect, fish counts (one column per species, plus Total.density, Species.richness, Total.biomass and trophic-group biomass), and benthic cover (coral genera, Coral.genera, other benthic categories) | scripts 01–07, 09 |
| `fish coral relationship change table 2.csv` | Summary of how each fish–coral relationship changed through time, one row per coral genus and fish species | script 08 |

Column names are converted to R-safe names when the data are read (spaces become dots).
