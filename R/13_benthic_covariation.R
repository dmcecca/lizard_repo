# =============================================================================
# 13_benthic_covariation.R
# Do the coral genera used as fish predictors covary with the other benthic components?
# Spearman correlations (all 336 transect records) among
#   - the 21 coral genera used in the INLA and BRT models, and
#   - the other benthic categories: macroalgae, turf algae, rubble, sand, soft coral, sponge, coralline algae.
#
# Produces (outputs/benthic_covariation/):
#   Benthic_correlation_matrix.pdf / .png   heatmap of all pairwise correlations (Figure S1)
#   Benthic_correlations_all_pairs.csv      every pair, Spearman rho and p-value
#   Coral_vs_other_benthos.csv              each coral genus against each non-coral category
#   Benthic_covariation_summary.csv         numbers for the text
# =============================================================================

# Sets the working directory to the project folder (the folder holding R and data) if it isn't already.
# This works when the script is open in RStudio; otherwise run setwd("path to lizard_repo") first.
if (!file.exists("R/00_setup.R")) setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
source("R/00_setup.R")

res_dir <- "benthic_covariation"
dat <- read_lizard()

OTHER_BENTHOS <- c("Macroalgae", "Turf.algae", "Rubble", "sand", "Soft.coral", "Sponge", "Coralline.algae")
corals <- intersect(CORAL_PREDICTORS, names(dat))
vars   <- c(corals, OTHER_BENTHOS)
x      <- dat |> select(all_of(vars)) |> mutate(across(everything(), as.numeric))
nice   <- \(v) gsub("\\.", " ", ifelse(v == "sand", "Sand", v))

# Spearman correlations (rank-based, suited to skewed cover data) with p-values
rho <- cor(x, method = "spearman", use = "pairwise.complete.obs")
pairs <- purrr::map_dfr(seq_len(ncol(x) - 1), \(i) purrr::map_dfr((i + 1):ncol(x), \(j) {
  ok <- complete.cases(x[[i]], x[[j]])
  tibble(var1 = vars[i], var2 = vars[j], n = sum(ok), rho = rho[i, j],
         p = suppressWarnings(cor.test(x[[i]], x[[j]], method = "spearman", exact = FALSE)$p.value))
})) |> mutate(p_adj = p.adjust(p, method = "holm"),
              type = case_when(var1 %in% corals & var2 %in% corals ~ "coral-coral",
                               var1 %in% corals | var2 %in% corals ~ "coral-other benthos", TRUE ~ "other-other"))
write.csv(pairs, out_path(res_dir, "Benthic_correlations_all_pairs.csv"), row.names = FALSE)

coral_other <- pairs |> filter(type == "coral-other benthos") |>
  mutate(coral = ifelse(var1 %in% corals, var1, var2), benthos = ifelse(var1 %in% corals, var2, var1)) |>
  select(coral, benthos, n, rho, p, p_adj) |> arrange(desc(abs(rho)))
write.csv(coral_other, out_path(res_dir, "Coral_vs_other_benthos.csv"), row.names = FALSE)

summary_tab <- tibble(
  quantity = c("coral genera", "other benthic categories", "coral x other-benthos pairs",
               "pairs with |rho| >= 0.3", "pairs with |rho| >= 0.5", "largest |rho|",
               "median |rho|", "pairs significant after Holm correction"),
  value = c(length(corals), length(OTHER_BENTHOS), nrow(coral_other),
            sum(abs(coral_other$rho) >= 0.3), sum(abs(coral_other$rho) >= 0.5),
            round(max(abs(coral_other$rho)), 2), round(median(abs(coral_other$rho)), 2), sum(coral_other$p_adj < 0.05)))
write.csv(summary_tab, out_path(res_dir, "Benthic_covariation_summary.csv"), row.names = FALSE)
print(summary_tab); print(head(coral_other, 10))

# Heatmap: corals first, then the other benthic categories
ord  <- vars
long <- as.data.frame(as.table(rho)) |> setNames(c("a", "b", "rho")) |>
  mutate(a = factor(a, levels = ord), b = factor(b, levels = rev(ord)))
p <- ggplot(long, aes(a, b, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  geom_text(data = filter(long, as.character(a) != as.character(b), abs(rho) >= 0.3),
            aes(label = sprintf("%.1f", rho)), size = 2.1) +
  geom_vline(xintercept = length(corals) + 0.5, linewidth = 0.5) +
  geom_hline(yintercept = length(OTHER_BENTHOS) + 0.5, linewidth = 0.5) +
  scale_fill_gradient2(low = "orange", mid = "white", high = "blue", midpoint = 0, limits = c(-1, 1),
                       name = "Spearman rho") +
  scale_x_discrete(position = "top", labels = nice) + scale_y_discrete(labels = nice) +
  labs(x = NULL, y = NULL) + theme_minimal(base_size = 10) +
  theme(axis.text.x.top = element_text(angle = 45, hjust = 0, face = "italic"),
        axis.text.y = element_text(face = "italic"), panel.grid = element_blank(), legend.position = "bottom",
        plot.margin = margin(10, 45, 10, 10))
ggsave(out_path(res_dir, "Benthic_correlation_matrix.pdf"), p, width = 10, height = 10)
ggsave(out_path(res_dir, "Benthic_correlation_matrix.png"), p, width = 10, height = 10, dpi = 300)
message("Done: written to ", file.path(OUT_DIR, res_dir))
