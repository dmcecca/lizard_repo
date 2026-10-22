# =============================================================================
# 08_figure6_relationship_change.R
# Figure 6: change in each fish-coral relationship through the disturbance-recovery cycle.
#
# Input:   data/fish coral relationship change table 2.csv
#          A summary table of slope trajectories, compiled by hand from the year-specific slopes (see script 07).
#          Row 1 = exposure zone, row 2 = trajectory (Stable, Stronger, Weaker, Change, Bounce back),
#          columns A-B = coral genus and fish species. Each cell holds + (positive), - (negative),
#          n (neutral: interval spans zero) or a sequence such as "+ to -".
# Outputs (outputs/figure6/):
#   coral_fish_relationship_tile_plot.(pdf|png)       Figure 6
#   relationship_change_counts.csv                    counts per exposure zone and trajectory
# =============================================================================

source("R/00_setup.R")
require_pkgs(c("readr", "cowplot")); library(readr)

res_dir <- "figure6"
raw <- read_csv(CHANGE_TABLE_FILE, col_names = FALSE, show_col_types = FALSE)

# ---- 1. Column names from the two header rows ----------------------------------
fill_right <- function(x) { last <- NA_character_
  for (i in seq_along(x)) { if (!is.na(x[i]) && x[i] != "") last <- x[i] else if (!is.na(last)) x[i] <- last }
  x }
header1 <- fill_right(as.character(unlist(raw[1, ])))
header2 <- as.character(unlist(raw[2, ]))
new_names <- c("coral_genus", "fish",
               vapply(3:22, \(i) paste(str_squish(header1[i]), str_squish(header2[i]), sep = "_") |>
                        str_to_lower() |> str_replace_all("\\s+", "_"), character(1)))

dat_wide <- raw[-c(1, 2), 1:22]; names(dat_wide) <- new_names
dat_wide <- dat_wide |>
  mutate(coral_genus = na_if(as.character(coral_genus), ""), fish = na_if(str_squish(as.character(fish)), "")) |>
  fill(coral_genus) |> filter(!is.na(fish)) |> mutate(across(-c(coral_genus, fish), as.character))

# ---- 2. Long format ---------------------------------------------------------------
dat_long <- dat_wide |>
  pivot_longer(-c(coral_genus, fish), names_to = "exp_traj", values_to = "symbol") |>
  mutate(symbol = na_if(str_squish(symbol), "")) |> filter(!is.na(symbol)) |>
  separate(exp_traj, into = c("exposure", "trajectory"), sep = "_(?=stable|stronger|weaker|change|bounce_back$)") |>
  mutate(exposure = factor(exposure, levels = c("lagoon", "sheltered", "oblique", "exposed"),
                           labels = c("Lagoon", "Sheltered", "Oblique", "Exposed")),
         trajectory = factor(trajectory, levels = c("stable", "stronger", "weaker", "change", "bounce_back"),
                             labels = c("Stable", "Stronger", "Weaker", "Change", "Bounce back")))

write.csv(dat_long |> count(exposure, trajectory) |> pivot_wider(names_from = trajectory, values_from = n, values_fill = 0),
          out_path(res_dir, "relationship_change_counts.csv"), row.names = FALSE)

# Row positions follow the order of the table (top row at the top of the plot)
fish_order <- dat_wide |> select(coral_genus, fish) |> distinct() |> mutate(row_id = row_number())
n_fish <- nrow(fish_order); fish_order$y_plot <- n_fish - fish_order$row_id + 1
dat_long <- dat_long |> left_join(fish_order, by = c("coral_genus", "fish")) |>
  mutate(coral_genus = factor(coral_genus, levels = unique(dat_wide$coral_genus)))

# ---- 3. Colour classes -------------------------------------------------------------
classify_symbol <- function(x) {
  x <- str_squish(x)
  if (is.na(x) || x == "") return(NA_character_)
  if (x == "+") return("Positive")
  if (x == "-") return("Negative")
  if (x == "n") return("Neutral")
  if (str_detect(x, "to")) return("Change")     # "+ to -" or "- to +": a single change of direction
  "Mixed"                                       # sequences such as "+ - +", "n + n"
}
dat_long <- dat_long |>
  mutate(relation_class = factor(map_chr(symbol, classify_symbol),
                                 levels = c("Negative", "Neutral", "Positive", "Change", "Mixed")),
         symbol_plot = str_replace_all(symbol, "to", "  "))

cell_cols <- c(Negative = "#d73027", Neutral = "#d9d9d9", Positive = "#1a9850", Change = "#fdae61", Mixed = "#8073ac")

# ---- 4. Coral labels and separators -------------------------------------------------
coral_bounds <- fish_order |> group_by(coral_genus) |>
  summarise(ymin = min(y_plot), ymax = max(y_plot), ymid = mean(y_plot), .groups = "drop") |>
  arrange(match(coral_genus, unique(dat_wide$coral_genus)))
coral_seps <- coral_bounds |> mutate(ysep = ymin - 0.5) |> filter(ysep >= 0.5)

coral_label_plot <- ggplot() +
  geom_segment(data = coral_seps, aes(x = 0.92, xend = 1.02, y = ysep, yend = ysep), colour = "grey80", linewidth = 0.3) +
  geom_text(data = coral_bounds, aes(x = 0.90, y = ymid, label = coral_genus), hjust = 1, size = 5, fontface = "italic") +
  scale_x_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, n_fish + 0.5), expand = c(0, 0)) +
  coord_cartesian(clip = "off") + theme_void() + theme(plot.margin = margin(t = 45, r = 0, b = 40, l = 2))

# ---- 5. Main tile plot ----------------------------------------------------------------
main_plot <- ggplot(dat_long, aes(trajectory, y_plot)) +
  geom_tile(aes(fill = relation_class), colour = "black", linewidth = 0.35) +
  geom_text(data = filter(dat_long, relation_class %in% c("Change", "Mixed")), aes(label = symbol_plot), size = 5) +
  geom_segment(data = coral_seps, aes(x = 0.5, xend = 5.5, y = ysep, yend = ysep), inherit.aes = FALSE,
               colour = "grey80", linewidth = 0.3) +
  facet_wrap(~ exposure, nrow = 1) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(breaks = fish_order$y_plot, labels = fish_order$fish, limits = c(0.5, n_fish + 0.5), expand = c(0, 0)) +
  scale_fill_manual(values = cell_cols, na.value = "white", drop = FALSE,
                    guide = guide_legend(override.aes = list(colour = "black", linewidth = 0.35))) +
  coord_cartesian(clip = "off") + labs(x = NULL, y = NULL, fill = "Direction") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(size = 11, face = "bold", angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 11, face = "italic"), legend.position = "bottom",
        legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        panel.spacing.x = unit(1, "lines"), plot.margin = margin(t = 45, r = 10, b = 40, l = 0))

fig6 <- cowplot::plot_grid(coral_label_plot, main_plot, ncol = 2, align = "h", axis = "tb", rel_widths = c(0.16, 0.84))
ggsave(out_path(res_dir, "coral_fish_relationship_tile_plot.pdf"), fig6, width = 18, height = 24)
ggsave(out_path(res_dir, "coral_fish_relationship_tile_plot.png"), fig6, width = 18, height = 24, dpi = 150)
message("Done: Figure 6 written to ", file.path(OUT_DIR, res_dir))
