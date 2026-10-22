# =============================================================================
# 00_setup.R
# Settings and helper functions shared by all scripts.
# Every script starts with:  source("R/00_setup.R")
# Open lizard-fish-coral.Rproj in RStudio first, so the working directory is the project folder.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(tibble); library(stringr); library(ggplot2)
})

# ---- 1. Files and folders ---------------------------------------------------
DATA_FILE <- file.path("data", "lizard.fish.benthos.latest.csv")   # see data/README.md
CHANGE_TABLE_FILE <- file.path("data", "fish coral relationship change table 2.csv")   # input to Figure 6
OUT_DIR <- "outputs"

out_path <- function(...) {
  p <- file.path(OUT_DIR, ...)
  dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)
  p
}

# ---- 2. Global analysis settings -------------------------------------------
SEED_MASTER <- 1234

# Exposure level order:
#   alphabetical (Exposed first)  -> brms GLMMs
#   gradient order                -> ordination, PERMANOVA, SIMPER, INLA slopes
EXPOSURE_LEVELS_GRADIENT <- c("Exposed", "Oblique", "Sheltered", "Lagoon")

exposure_cols <- c(Exposed = "#F8766D", Lagoon = "#7CAE00", Oblique = "#00BFC4", Sheltered = "#C77CFF")
year_alpha    <- c("2011" = 0.25, "2015" = 0.50, "2017" = 0.75, "2020" = 1.00)

# Column ranges (names after make.names())
FISH_FIRST  <- "Abudefduf.sexfasciatus"
FISH_LAST   <- "Zebrasoma.veliferum"
CORAL_FIRST <- "Acanthastrea"
CORAL_LAST  <- "Heliopora"
METRIC_FIRST <- "Total.density"; METRIC_LAST <- "Total.biomass"      # density, richness, biomass
TROPHIC_FIRST <- "B_Benthic.invertivore"; TROPHIC_LAST <- "B_Corallivore"

BENTHIC_VARS <- c("Coralline.algae", "Rubble", "sand", "Soft.coral", "Sponge",
                  "Total.scl.coral.cover", "Macroalgae", "Turf.algae", "Coral.genera")
BENTHIC_EXCLUDE_SITE <- "Resort"   # dropped for the benthic ordination only

# Coral predictors used in the BRT and INLA models (21 genera, including fire coral)
CORAL_PREDICTORS <- c(
  "Acropora", "Montipora", "Astreopora", "Dipsastraea", "Pocillopora", "Echinopora",
  "Diploastrea", "Porites", "Stylophora", "Physogyra", "Acanthastrea", "Mycedium",
  "Leptoria", "Coeloseris", "Favites", "Isopora", "Platygyra", "Pavona", "Psammocora",
  "Fire.coral", "Goniopora"
)

# ---- 3. Helpers --------------------------------------------------------------
read_lizard <- function(file = DATA_FILE) {
  if (!file.exists(file)) {
    stop("Data file not found: ", file, "\nSee data/README.md, or edit DATA_FILE in R/00_setup.R.", call. = FALSE)
  }
  x <- read.csv(file, strip.white = TRUE, check.names = FALSE)
  names(x) <- make.names(names(x), unique = TRUE)
  miss <- setdiff(c("Site", "Year", "Exposure", "Depth", "Transect"), names(x))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
  x
}

col_range <- function(df, first, last) {
  i <- match(first, names(df)); j <- match(last, names(df))
  if (is.na(i) || is.na(j) || i > j) stop("Cannot find column range '", first, "' to '", last, "'.", call. = FALSE)
  names(df)[i:j]
}

safe_file <- function(x) x |> str_replace_all("[^A-Za-z0-9]+", "_") |> str_replace_all("^_|_$", "")

require_pkgs <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Please install: ", paste(miss, collapse = ", "),
                         "\n  (INLA is not on CRAN: https://www.r-inla.org/download-install)", call. = FALSE)
  invisible(TRUE)
}

# Save CSV always; .docx too if flextable is installed
save_table <- function(df, name, docx = TRUE, caption = NULL) {
  utils::write.csv(df, out_path(paste0(name, ".csv")), row.names = FALSE)
  if (docx && requireNamespace("flextable", quietly = TRUE)) {
    ft <- flextable::autofit(flextable::flextable(as.data.frame(df)))
    if (!is.null(caption)) ft <- flextable::set_caption(ft, caption)
    flextable::save_as_docx(ft, path = out_path(paste0(name, ".docx")))
  }
  invisible(df)
}

# ---- 4. Community helper (Fig 3, Table 3, Table S3) ---------------------------
# Aggregate transects to Site x Year x Exposure means, transform, optionally drop
# low-variance variables, return Bray-Curtis distances.
#   "4th-root": (x + 1)^0.25 (fish);  "log": log(x + 1) (benthic categories, coral genera);
#   "none": untransformed means (reproduces the Year-only PERMANOVA values quoted in the Results text)
prep_community <- function(data, vars, transformation = c("4th-root", "log", "none"),
                           variance_filter = TRUE, var_cutoff = 0.01) {
  transformation <- match.arg(transformation)
  agg <- data |>
    group_by(Site, Year, Exposure) |>
    summarise(across(all_of(vars), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
  meta    <- agg |> select(Site, Year, Exposure)
  species <- agg |> select(all_of(vars))
  species_t <- switch(transformation,
    "4th-root" = mutate(species, across(everything(), \(x) (x + 1)^0.25)),
    "log"      = mutate(species, across(everything(), \(x) log(x + 1))),
    "none"     = species)
  if (variance_filter) species_t <- species_t |> select(where(\(x) stats::var(x) > var_cutoff))
  list(meta = meta, species = species_t, dist = vegan::vegdist(species_t, method = "bray"))
}
