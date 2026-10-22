# =============================================================================
# 10_condensed_glmm_tables.R
# Condensed contrast tables (Table S2), built from the models saved by scripts 01 and 02:
#   Exposure contrasts averaged over the four survey years (6 rows)
#   Year contrasts averaged over the four exposure zones (6 rows)
# Cells show the rate ratio [95% CrI] and Pr(RR>1). Bold in the Word file = 95% CrI excludes 1.
#
# Produces:  outputs/condensed_tables/Table_S2_fish_metrics_condensed.(csv|docx)
#            outputs/condensed_tables/Table_S2_trophic_groups_condensed.(csv|docx)
# Run 01 and 02 first.
# =============================================================================

source(file.path("R", "00_setup.R"))
require_pkgs(c("brms", "emmeans", "coda", "broom", "flextable"))

FISH_MODEL_FOLDER    <- file.path(OUT_DIR, "glmm_fish_metrics", "brms_fish_models")
TROPHIC_MODEL_FOLDER <- file.path(OUT_DIR, "glmm_trophic", "brms_trophic_models")
OUT_FOLDER           <- file.path(OUT_DIR, "condensed_tables")

library(brms); library(emmeans); library(coda); library(broom)

CAPTION <- paste(
  "Pairwise contrasts among exposure zones (averaged over survey years) and among survey years (averaged",
  "over exposure zones). Values are rate ratios (RR) of estimated means, first level / second level,",
  "with 95% credible intervals [CrI] and exceedance probability Pr(RR>1). Averages are on the log scale with",
  "equal weights. Bold: 95% CrI excludes 1.")

fmt_cell <- function(rr, lo, hi, ep) sprintf("%.2f [%.2f, %.2f]; Pr(RR>1)=%.2f", rr, lo, hi, ep)

# Posterior summary of each contrast (link scale -> rate ratio)
summarise_contrasts <- function(ct, ci = 0.95) {
  a <- (1 - ci) / 2
  draws <- as.data.frame(do.call(rbind, as.mcmc(ct)))
  lab <- broom::tidy(ct) |> mutate(.col = colnames(draws))
  imap_dfr(draws, \(v, nm) tibble(
    RR = exp(mean(v)), RR_low = exp(quantile(v, a, names = FALSE)), RR_high = exp(quantile(v, 1 - a, names = FALSE)),
    EP = mean(v > 0), .col = nm)) |> left_join(lab, by = ".col") |> select(contrast, RR, RR_low, RR_high, EP)
}

main_effects_one <- function(model, name) {
  clean <- function(x) gsub(" - ", " \u2212 ", gsub("Year|Exposure", "", x), fixed = TRUE)
  one <- function(by, scope) {
    ct <- contrast(emmeans(model, as.formula(paste("~", by)), type = "link", re_formula = NA), method = "pairwise")
    summarise_contrasts(ct) |> transmute(Model = name, Scope = scope, Contrast = clean(contrast), RR, RR_low, RR_high, EP,
                                         Excludes_1 = RR_low > 1 | RR_high < 1)
  }
  bind_rows(one("Exposure", "Exposure (averaged over years)"), one("Year", "Year (averaged over exposures)"))
}

make_table <- function(folder, name, label_map = NULL) {
  files <- list.files(folder, pattern = "\\.rds$", full.names = TRUE)
  if (!length(files)) stop("No .rds files found in: ", folder, "\nRun scripts 01 and 02 first.")
  ids <- sub("\\.rds$", "", basename(files))
  if (!is.null(label_map)) {                       # match file names to nice labels, keep the order of label_map
    key <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))
    m <- match(key(names(label_map)), key(ids)); keep <- !is.na(m)
    files <- files[m[keep]]; ids <- names(label_map)[keep]; labels <- unname(label_map)[keep]
  } else labels <- ids
  message("Models found: ", paste(labels, collapse = ", "))
  long <- map2_dfr(files, labels, \(f, l) { message("Processing ", l); main_effects_one(readRDS(f), l) })
  wide <- long |> mutate(Cell = fmt_cell(RR, RR_low, RR_high, EP), Model = factor(Model, levels = labels)) |>
    select(Scope, Contrast, Model, Cell) |> pivot_wider(names_from = Model, values_from = Cell)
  dir.create(OUT_FOLDER, showWarnings = FALSE, recursive = TRUE)
  write.csv(wide, file.path(OUT_FOLDER, paste0(name, ".csv")), row.names = FALSE)
  write.csv(long, file.path(OUT_FOLDER, paste0(name, "_numeric.csv")), row.names = FALSE)
  if (requireNamespace("flextable", quietly = TRUE)) {
    ft <- flextable::merge_v(flextable::autofit(flextable::flextable(as.data.frame(wide))), j = "Scope")
    for (l in labels) {
      r <- match(paste(long$Scope, long$Contrast)[long$Model == l & long$Excludes_1], paste(wide$Scope, wide$Contrast))
      if (length(r)) ft <- flextable::bold(ft, i = r, j = l, bold = TRUE)
    }
    ft <- flextable::set_caption(ft, CAPTION)
    flextable::save_as_docx(ft, path = file.path(OUT_FOLDER, paste0(name, ".docx")))
  } else message("flextable not installed: wrote CSV only (install.packages('flextable') for Word output)")
  message("Saved ", name, " to ", OUT_FOLDER)
  invisible(wide)
}

fish_labels <- c(Total_density = "Total density", Species_richness = "Species richness", Total_biomass = "Total biomass")
trophic_labels <- c(B_Benthic.invertivore = "Benthic invertivore", B_Carnivore = "Carnivore", B_Omnivore = "Omnivore",
  B_Piscivore = "Piscivore", B_Planktivore = "Planktivore", B_Grazer = "Grazer - All", B_Detritivore = "Grazer - Detritivore",
  B_Cropper = "Grazer - Cropper", B_Browser = "Grazer - Browser", B_Scraper = "Grazer - Scraper", B_Excavator = "Excavator",
  B_Farmer = "Farmer - All", B_Intensive.farmer = "Intensive farmer", B_Extensive.farmer = "Extensive farmer",
  B_Indeterminate.farmer = "Indeterminate farmer", B_Corallivore = "Corallivore - All",
  B_Obligate.corallivore = "Obligate corallivore", B_Facultative.corallivore = "Facultative corallivore")

make_table(FISH_MODEL_FOLDER,    "Table_S2_fish_metrics_condensed",   fish_labels)
make_table(TROPHIC_MODEL_FOLDER, "Table_S2_trophic_groups_condensed", trophic_labels)
