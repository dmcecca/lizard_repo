# =============================================================================
# 07_temporal_slope_template.R
# Example: change through time in one fish-coral relationship, using Year-specific slopes within one
# exposure zone (INLA negative binomial). Edit FISH, CORAL and EXPOSURE to run another pair.
#
# Output:  outputs/temporal_slopes/slopes_<fish>_<coral>_<exposure>.csv
#          (slope mean, 95% interval and direction [+, -, n = interval spans zero] for each year)
# Needs:   INLA
# =============================================================================

source("R/00_setup.R")
require_pkgs(c("INLA")); library(INLA)

FISH     <- "Chaetodon.aureofasciatus"   # example only - replace
CORAL    <- "Acropora"                   # example only - replace
EXPOSURE <- "Exposed"
NSAMP    <- 400
set.seed(SEED_MASTER)

d <- read_lizard() |>
  mutate(SiteTransect = interaction(Site, Transect, drop = TRUE), Year = factor(Year)) |>
  filter(Exposure == EXPOSURE) |>
  transmute(Year, SiteTransect, fish = as.integer(.data[[FISH]]), coral = as.numeric(.data[[CORAL]])) |>
  filter(!is.na(fish), !is.na(coral)) |> droplevels()
d$coral_sc <- as.numeric(scale(d$coral))

mod <- inla(fish ~ 1 + Year + coral_sc + coral_sc:Year + f(SiteTransect, model = "iid"),
            data = d, family = "nbinomial", control.compute = list(config = TRUE))

# Year-specific slope = coral_sc + coral_sc:Year<level> (reference year = first level)
samp <- inla.posterior.sample(NSAMP, mod)
fe   <- rownames(mod$summary.fixed)
rows <- rownames(samp[[1]]$latent)[grepl(paste(fe, collapse = "|"), rownames(samp[[1]]$latent))]
B <- do.call(rbind, lapply(samp, \(s) s$latent[rows, 1])); colnames(B) <- fe
yrs <- levels(d$Year)
slopes <- sapply(yrs, \(y) if (y == yrs[1]) B[, "coral_sc"] else B[, "coral_sc"] + B[, grep(paste0("Year", y), fe, value = TRUE)[grep("coral_sc", grep(paste0("Year", y), fe, value = TRUE))]])
summ <- tibble(Year = yrs, mean = colMeans(slopes), lo = apply(slopes, 2, quantile, 0.025), hi = apply(slopes, 2, quantile, 0.975),
               direction = case_when(lo > 0 ~ "+", hi < 0 ~ "-", TRUE ~ "n"))
print(summ)
write.csv(summ, out_path("temporal_slopes", paste0("slopes_", safe_file(FISH), "_", CORAL, "_", EXPOSURE, ".csv")), row.names = FALSE)
