# =============================================================================
# helpers_figure5.R
# Functions for Figure 5 (species x coral-genus heatmap) and the species relative-influence table.
# Used by 05_brt_species.R and 05b_redraw_figure5.R.
# =============================================================================

# Colours show the signed influence of each coral genus (orange = negative, blue = positive).
# Exposure has no direction, so it is shown as a last grey column, shaded by relative influence (%).
draw_figure5 <- function(results, out_dir, row_height = 0.20, species_text = 9, width = 11, split = TRUE) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  mat <- results
  mat$Species   <- trimws(mat$Species)
  mat$Predictor <- trimws(mat$Predictor)
  species_order <- sort(unique(mat$Species))
  coral_order   <- sort(setdiff(unique(mat$Predictor), "Exposure"))
  pred_order    <- c(coral_order, intersect("Exposure", mat$Predictor))

  grid <- expand.grid(Species = species_order, Predictor = pred_order, stringsAsFactors = FALSE) |>
    dplyr::left_join(mat[, c("Species", "Predictor", "SignedInfluence", "BootMean")], by = c("Species", "Predictor"))
  grid$Signed    <- ifelse(grid$Predictor == "Exposure", NA, grid$SignedInfluence)
  grid$ExposureInfluence <- ifelse(grid$Predictor == "Exposure", grid$BootMean, NA)
  grid$Predictor <- factor(grid$Predictor, levels = pred_order)      # keeps Exposure as the last column
  lim  <- max(abs(grid$Signed), na.rm = TRUE)
  nice <- function(x) gsub("\\.", " ", x)

  make_heat <- function(sp) {
    d <- grid |> dplyr::filter(Species %in% sp) |> dplyr::mutate(Species = factor(Species, levels = rev(sp)))
    p <- ggplot2::ggplot(d, ggplot2::aes(Predictor, Species)) +
      ggplot2::geom_tile(ggplot2::aes(fill = Signed), colour = "grey80", linewidth = 0.25) +
      # species names are repeated on the right so rows can be matched from either side
      ggplot2::geom_text(data = dplyr::distinct(d, Species),
                         ggplot2::aes(x = length(pred_order) + 0.6, y = Species, label = nice(Species)),
                         inherit.aes = FALSE, hjust = 0, size = species_text / ggplot2::.pt,
                         colour = "black", fontface = "italic") +
      ggplot2::scale_fill_gradient2(low = "orange", mid = "white", high = "blue", midpoint = 0,
                                    limits = c(-lim, lim), name = "Influence (coral genera)", na.value = "white")
    if ("Exposure" %in% pred_order)
      p <- p + ggplot2::geom_tile(data = dplyr::filter(d, !is.na(ExposureInfluence)),
                                  ggplot2::aes(alpha = ExposureInfluence), fill = "grey20",
                                  colour = "grey80", linewidth = 0.25) +
        ggplot2::scale_alpha_continuous(range = c(0.08, 1), name = "Exposure influence (%)")
    p +
      ggplot2::scale_y_discrete(expand = c(0, 0), labels = nice) +
      ggplot2::scale_x_discrete(expand = c(0, 0), position = "top", labels = nice) +
      ggplot2::coord_cartesian(xlim = c(0.5, length(pred_order) + 0.5), clip = "off") +
      ggplot2::labs(x = NULL, y = NULL) + ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(axis.text.x.top = ggplot2::element_text(face = "italic", angle = 45, hjust = 0, vjust = 0, size = 10, colour = "black"),
                     axis.text.y = ggplot2::element_text(face = "italic", size = species_text, colour = "black"),
                     panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
                     legend.position = "bottom", legend.title = ggplot2::element_text(size = 10, face = "bold"),
                     legend.text = ggplot2::element_text(size = 9), legend.key.width = ggplot2::unit(1.6, "cm"),
                     plot.margin = ggplot2::margin(10, 190, 10, 10))
  }

  save_fig <- function(p, name, n_rows) {
    h <- n_rows * row_height + 3.5                    # extra height for the top labels and legends
    ggplot2::ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = width, height = h, device = "pdf", limitsize = FALSE)
    ggplot2::ggsave(file.path(out_dir, paste0(name, ".png")), p, width = width, height = h, dpi = 300, limitsize = FALSE)
  }
  save_fig(make_heat(species_order), "Figure5_tall", length(species_order))
  if (split) {                                         # the same figure in two halves
    half <- ceiling(length(species_order) / 2)
    save_fig(make_heat(species_order[seq_len(half)]), "Figure5_part1", half)
    save_fig(make_heat(species_order[-seq_len(half)]), "Figure5_part2", length(species_order) - half)
  }
  invisible(NULL)
}

# Relative-influence table for all species: species A-Z, predictors in order of decreasing influence
species_influence_table <- function(results) {
  tab <- results[order(results$Species, -results$Influence), ]
  tab$Species   <- gsub("\\.", " ", tab$Species)
  tab$Predictor <- gsub("\\.", " ", tab$Predictor)
  tab$SignedInfluence[tab$Predictor == "Exposure"] <- NA        # exposure has no direction
  data.frame(Species = tab$Species, Predictor = tab$Predictor,
             `Influence (%)` = round(tab$Influence, 2), `Bootstrap mean (%)` = round(tab$BootMean, 2),
             `Lower 95% CI` = round(tab$Lower95, 2), `Upper 95% CI` = round(tab$Upper95, 2),
             `Successful bootstrap fits` = tab$NSuccess, `Signed influence` = round(tab$SignedInfluence, 2),
             check.names = FALSE)
}
