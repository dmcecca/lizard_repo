# =============================================================================
# 03_community_composition.R
# Composition of fish, benthic and coral assemblages in relation to year, exposure and site.
#
# Outputs (outputs/community/):
#   Figure3_nMDS.(png|pdf|tiff)                     Figure 3 (nMDS and fitted vectors)
#   Table3_PERMANOVA                                 Table 3 (Year + Exposure + Site)
#   PERMANOVA_year_only_factor, betadisper_year      Year-only PERMANOVA and dispersion tests
#   PERMANOVA_pairwise_all                           pairwise PERMANOVA (years, exposure zones)
#   TableS3_SIMPER_between_exposure_by_year_top10    SIMPER, every year and pair of exposure zones
#   TableS3_SIMPER_top10_pooled                      SIMPER, top 10 per assemblage averaged over all comparisons
#
# Data: transects are averaged to Site x Year means (56 rows; 52 for benthic, which excludes the Resort site).
#       Fish are (x + 1)^0.25 transformed; benthic categories and coral genera are log(x + 1) transformed.
#       Variables with variance <= 0.01 are dropped for the nMDS and PERMANOVA (not for SIMPER).
#       Bray-Curtis dissimilarity is used throughout.
# Run time: minutes.
# =============================================================================

# Sets the working directory to the project folder (the folder holding R and data) if it isn't already.
# This works when the script is open in RStudio; otherwise run setwd("path to lizard_repo") first.
if (!file.exists("R/00_setup.R")) setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
source("R/00_setup.R")
require_pkgs(c("vegan", "permute", "ggrepel", "patchwork"))
library(vegan); library(ggrepel); library(patchwork)

res_dir <- "community"
set.seed(SEED_MASTER)       # permutation tests, envfit, nMDS random starts
N_PERM  <- 999
YEAR_AS_FACTOR <- FALSE     # FALSE: Year is numeric (1 df) in Table 3; TRUE: Year is a 4-level factor (3 df)
VECTOR_R_CUTOFF <- 0.5      # fitted vectors with correlation above this value are labelled

# ---- 1. Data ----------------------------------------------------------------
lizard <- read_lizard() |>
  mutate(Site = factor(Site), Exposure = factor(Exposure, levels = EXPOSURE_LEVELS_GRADIENT),
         Year = as.numeric(as.character(Year)), Depth = factor(Depth), Transect = factor(Transect))

fish_vars  <- col_range(lizard, FISH_FIRST, FISH_LAST)
coral_vars <- col_range(lizard, CORAL_FIRST, CORAL_LAST)
stopifnot(all(BENTHIC_VARS %in% names(lizard)))
lizard_benthic <- lizard |> filter(Site != BENTHIC_EXCLUDE_SITE)

assemblages <- list(
  Fish    = list(data = lizard,         vars = fish_vars,    transf = "4th-root", titles = c("a) Fish", "b) Fish")),
  Benthic = list(data = lizard_benthic, vars = BENTHIC_VARS, transf = "log",      titles = c("c) Benthic categories", "d) Benthic categories")),
  Coral   = list(data = lizard,         vars = coral_vars,   transf = "log",      titles = c("e) Coral genera", "f) Coral genera")))

prepped <- imap(assemblages, \(a, nm) prep_community(a$data, a$vars, a$transf, variance_filter = TRUE))

# Table 3 benthic block excludes generic coral richness (Coral.genera)
BENTHIC_VARS_PERMANOVA <- setdiff(BENTHIC_VARS, "Coral.genera")
prepped_t3 <- prepped
prepped_t3$Benthic <- prep_community(lizard_benthic, BENTHIC_VARS_PERMANOVA, "log", variance_filter = TRUE)

# Untransformed means, used for the Year-only tests in 3b-3d as well as the transformed data
prepped_raw <- imap(assemblages, \(a, nm) prep_community(a$data, a$vars, "none", variance_filter = FALSE))

# ---- 2. nMDS (Figure 3) -----------------------------------------------------
# Panels b, d and f show only the fitted vectors with a correlation above VECTOR_R_CUTOFF (weaker vectors are not drawn).
# Fish and coral labels are italic (species and genus names); dots in names are shown as spaces.
nmds_analysis <- function(a, p, label, italic_labels = FALSE) {
  nmds <- vegan::metaMDS(p$dist, k = 2, trymax = 500, trace = 0)
  message(label, " nMDS stress: ", round(nmds$stress, 3))
  nmds_data <- cbind(p$meta, as.data.frame(nmds$points)) |> mutate(Year = factor(Year, levels = sort(unique(Year))))
  
  mds_plot <- ggplot(nmds_data, aes(MDS1, MDS2)) +
    geom_point(aes(colour = Exposure, shape = Exposure, alpha = Year), size = 3) +
    stat_ellipse(aes(group = Exposure, linetype = Exposure), colour = "black", linewidth = 0.8) +
    scale_colour_manual(values = exposure_cols, drop = FALSE, name = "Exposure") +
    scale_alpha_manual(name = "Year", values = year_alpha) +
    labs(title = a$titles[1], x = "NMDS Dimension 1", y = "NMDS Dimension 2") + theme_classic()
  
  fit <- vegan::envfit(nmds, p$species, perm = N_PERM)
  vec <- as.data.frame(vegan::scores(fit, display = "vectors"))
  vec$species     <- gsub("\\.", " ", rownames(vec))
  vec$species     <- dplyr::recode(vec$species, "sand" = "Sand", "Coral genera" = "Generic richness",
                                   "Total scl coral cover" = "Hard coral")   # display names for benthic variables
  vec$correlation <- sqrt(vec$NMDS1^2 + vec$NMDS2^2)
  vec <- vec[vec$correlation > VECTOR_R_CUTOFF, ]          # drop the weaker (previously grey) vectors
  
  vector_plot <- ggplot(vec) +
    geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                 arrow = grid::arrow(length = grid::unit(0.2, "cm")), linewidth = 0.7, colour = "black") +
    ggrepel::geom_text_repel(aes(NMDS1, NMDS2, label = species), size = 3, colour = "black",
                             fontface = if (italic_labels) "italic" else "plain",
                             max.overlaps = Inf, box.padding = 0.35, point.padding = 0.1,
                             min.segment.length = Inf, seed = 1) +
    scale_x_continuous(expand = expansion(mult = 0.15)) + scale_y_continuous(expand = expansion(mult = 0.1)) +
    labs(title = a$titles[2], x = "NMDS Dimension 1", y = "NMDS Dimension 2") + theme_classic()
  list(mds = mds_plot, vectors = vector_plot, stress = nmds$stress)
}

nmds <- list(Fish    = nmds_analysis(assemblages$Fish,    prepped$Fish,    "Fish",    italic_labels = TRUE),
             Benthic = nmds_analysis(assemblages$Benthic, prepped$Benthic, "Benthic"),
             Coral   = nmds_analysis(assemblages$Coral,   prepped$Coral,   "Coral",   italic_labels = TRUE))

# Layout: a and b each span the full width (stacked); c-d and e-f sit below as pairs
fig3 <- nmds$Fish$mds / nmds$Fish$vectors /
  (nmds$Benthic$mds | nmds$Benthic$vectors) /
  (nmds$Coral$mds   | nmds$Coral$vectors) +
  plot_layout(heights = c(1, 1.5, 1, 1), guides = "collect")
FIG3_W <- 10; FIG3_H <- 20
ggsave(out_path(res_dir, "Figure3_nMDS.png"),  fig3, width = FIG3_W, height = FIG3_H, dpi = 300)
ggsave(out_path(res_dir, "Figure3_nMDS.pdf"),  fig3, width = FIG3_W, height = FIG3_H)
ggsave(out_path(res_dir, "Figure3_nMDS.tiff"), fig3, width = FIG3_W, height = FIG3_H, dpi = 600, compression = "lzw")
write.csv(tibble(Dataset = names(nmds), Stress = map_dbl(nmds, "stress")),
          out_path(res_dir, "nMDS_stress.csv"), row.names = FALSE)

# ---- 3. PERMANOVA (Table 3) and dispersion -----------------------------------
# 3a. Table 3:  dist ~ Year + Exposure + Site  (sequential terms)
permanova_global <- imap_dfr(prepped_t3, function(p, nm) {
  meta <- p$meta
  if (YEAR_AS_FACTOR) meta$Year <- factor(meta$Year)
  as.data.frame(vegan::adonis2(p$dist ~ Year + Exposure + Site, data = meta, permutations = N_PERM)) |>
    rownames_to_column("Term") |> mutate(Dataset = nm, .before = 1)
})
save_table(permanova_global, file.path(res_dir, "Table3_PERMANOVA"), caption = "PERMANOVA results")

# 3a2. Table 3 with permutations that respect the repeated sampling of sites
#   Year:     permutations are restricted within sites (years are shuffled only among the same site's years)
#   Exposure: sites are the units; Site x Year means are averaged to one profile per site and whole sites are permuted
permanova_restricted <- imap_dfr(prepped_t3, function(p, nm) {
  meta <- p$meta; meta$Year <- factor(meta$Year); meta$Site <- factor(meta$Site)
  yr <- as.data.frame(vegan::adonis2(p$dist ~ Year, data = meta,
                                     permutations = permute::how(blocks = meta$Site, nperm = N_PERM))) |>
    rownames_to_column("Term") |> filter(Term == "Year") |> mutate(Test = "Year (permuted within sites)")
  site_means <- as.data.frame(p$species) |> mutate(Site = meta$Site) |> group_by(Site) |>
    summarise(across(everything(), mean), .groups = "drop")
  site_meta  <- meta |> distinct(Site, Exposure) |> arrange(Site)
  ex <- as.data.frame(vegan::adonis2(vegan::vegdist(select(site_means, -Site), "bray") ~ Exposure, data = site_meta,
                                     permutations = N_PERM)) |>
    rownames_to_column("Term") |> filter(Term == "Exposure") |> mutate(Test = "Exposure (whole sites permuted)")
  bind_rows(yr, ex) |> mutate(Dataset = nm, n_sites = nlevels(meta$Site), .before = 1)
})
save_table(permanova_restricted, file.path(res_dir, "Table3_PERMANOVA_restricted_permutations"),
           caption = "PERMANOVA with permutations restricted to account for repeated sampling of sites.")

# 3b. Year-only model, with Year as a 4-level factor
year_only <- function(p_list, tag) imap_dfr(p_list, function(p, nm) {
  as.data.frame(vegan::adonis2(p$dist ~ Year, data = mutate(p$meta, Year = factor(Year)), permutations = N_PERM)) |>
    rownames_to_column("Term") |> mutate(Dataset = nm, Transformation = tag, .before = 1)
})
permanova_year_only <- bind_rows(year_only(prepped, "as analysed (4th-root / log)"), year_only(prepped_raw, "untransformed"))
save_table(permanova_year_only, file.path(res_dir, "PERMANOVA_year_only_factor"), docx = FALSE)

# 3c. Homogeneity of multivariate dispersion among years (unrestricted permutations)
disp <- function(p_list, tag) imap_dfr(p_list, function(p, nm) {
  pt <- vegan::permutest(vegan::betadisper(p$dist, factor(p$meta$Year)), permutations = N_PERM)
  tibble(Dataset = nm, Transformation = tag, F = pt$tab[1, "F"], p = pt$tab[1, "Pr(>F)"])
})
dispersion <- bind_rows(disp(prepped, "as analysed (4th-root / log)"), disp(prepped_raw, "untransformed"))
save_table(dispersion, file.path(res_dir, "betadisper_year"), docx = FALSE)

# 3d. Pairwise PERMANOVA between years and between exposure zones (raw, Holm and BH-adjusted p)
#   Year pairs: permutations within sites.  Exposure pairs: site profiles, whole sites permuted.
pairwise_permanova <- function(p, dataset, group_var) {
  meta <- p$meta; g <- factor(meta[[group_var]]); full <- as.matrix(p$dist)
  map_dfr(utils::combn(levels(droplevels(g)), 2, simplify = FALSE), function(pr) {
    keep <- which(g %in% pr)
    if (group_var == "Year") {
      fit <- vegan::adonis2(as.dist(full[keep, keep]) ~ group, data = data.frame(group = droplevels(g[keep])),
                            permutations = permute::how(blocks = droplevels(factor(meta$Site[keep])), nperm = N_PERM))
    } else {
      sp <- as.data.frame(p$species)[keep, , drop = FALSE]; sp$Site <- droplevels(factor(meta$Site[keep]))
      site_means <- sp |> group_by(Site) |> summarise(across(everything(), mean), .groups = "drop")
      site_grp   <- meta[keep, ] |> mutate(Site = droplevels(factor(Site))) |> distinct(Site, .data[[group_var]]) |> arrange(Site)
      fit <- vegan::adonis2(vegan::vegdist(select(site_means, -Site), "bray") ~ group,
                            data = data.frame(group = droplevels(factor(site_grp[[group_var]]))), permutations = N_PERM)
    }
    tibble(Dataset = dataset, Factor = group_var, Level1 = pr[1], Level2 = pr[2],
           Df = fit$Df[1], SumOfSqs = fit$SumOfSqs[1], R2 = fit$R2[1], F = fit$F[1], p = fit$`Pr(>F)`[1])
  }) |> mutate(p_holm = p.adjust(p, "holm"), p_BH = p.adjust(p, "BH"))
}
pairwise_all <- bind_rows(
  imap_dfr(prepped_t3,  \(p, nm) pairwise_permanova(p, nm, "Year")),
  imap_dfr(prepped_t3,  \(p, nm) pairwise_permanova(p, nm, "Exposure")),
  imap_dfr(prepped_raw, \(p, nm) pairwise_permanova(p, nm, "Year")) |> mutate(Factor = "Year (untransformed)"))
save_table(pairwise_all, file.path(res_dir, "PERMANOVA_pairwise_all"), caption = "Pairwise PERMANOVA results")

# ---- 4. SIMPER between exposure zones within each year (Table S3) ---------------
# SIMPER is run for every pair of exposure zones in each year, on all variables (no variance filter).
simper_all_variables <- function(data, vars, transf, dataset) {
  map_dfr(sort(unique(data$Year)), function(y) {
    p <- prep_community(filter(data, Year == y), vars, transf, variance_filter = FALSE)
    sim <- vegan::simper(as.matrix(p$species), group = factor(p$meta$Exposure), permutations = 0)
    map_dfr(names(sim), function(cmp) {
      s <- sim[[cmp]]
      tibble(Dataset = dataset, Year = y, Comparison = cmp, Variable = s$species,
             AverageContribution = as.numeric(s$average[s$species]), SD = as.numeric(s$sd[s$species]),
             Ratio = as.numeric(s$ratio[s$species]), Ava = as.numeric(s$ava[s$species]),
             Avb = as.numeric(s$avb[s$species]), CumSum = as.numeric(s$cusum[s$species]))
    })
  }) |>
    group_by(Dataset, Year, Comparison) |>
    mutate(Percent = 100 * AverageContribution / sum(AverageContribution),
           Rank = rank(-AverageContribution, ties.method = "first")) |>
    ungroup()
}
simper_all <- imap_dfr(assemblages, \(a, nm) simper_all_variables(a$data, a$vars, a$transf, nm))

# 4a. Top 10 variables for every year and pair of exposure zones
simper_by_year <- simper_all |> filter(Rank <= 10) |> arrange(Dataset, Year, Comparison, Rank) |>
  select(-Percent, -Rank)
save_table(simper_by_year, file.path(res_dir, "TableS3_SIMPER_between_exposure_by_year_top10"),
           caption = "SIMPER: top contributors to between-exposure dissimilarity by year")

# 4b. Condensed table: top 10 variables per assemblage, averaged over all years and exposure pairs
simper_pooled <- simper_all |>
  group_by(Dataset, Variable) |>
  summarise(MeanContribution = mean(AverageContribution), PercentOfDissimilarity = mean(Percent),
            MeanRatio = mean(Ratio), TimesInTop10 = sum(Rank <= 10), Comparisons = n(), .groups = "drop") |>
  group_by(Dataset) |> arrange(desc(MeanContribution), .by_group = TRUE) |> slice_head(n = 10) |> ungroup() |>
  mutate(MeanContribution = round(MeanContribution, 4), PercentOfDissimilarity = round(PercentOfDissimilarity, 1),
         MeanRatio = round(MeanRatio, 2))
save_table(simper_pooled, file.path(res_dir, "TableS3_SIMPER_top10_pooled"),
           caption = "SIMPER: top 10 contributors to between-exposure dissimilarity, averaged over years and exposure pairs")
message("Done: Figure 3, Table 3, pairwise PERMANOVA and Table S3 written to ", file.path(OUT_DIR, res_dir))
