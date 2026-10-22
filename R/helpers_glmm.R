# =============================================================================
# helpers_glmm.R
# Functions shared by 01_fish_summary_glmms.R and 02_trophic_glmms.R.
# Model, one per response column:   Response ~ Year * Exposure + (1 | SiteDepthTransect)
# Family chosen automatically:  all values > 0 -> Gamma(log);  any zeros -> hurdle_lognormal()
# =============================================================================

require_pkgs(c("brms", "posterior", "emmeans", "coda", "broom", "cowplot"))
suppressPackageStartupMessages({ library(brms); library(posterior); library(emmeans); library(coda); library(cowplot) })
options(mc.cores = min(4, parallel::detectCores()))

priors_gamma <- c(
  prior(normal(0, 1.5), class = "b"), prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), class = "sd"), prior(gamma(1, 1), class = "shape"))
priors_hurdle_logn <- c(
  prior(normal(0, 1.5), class = "b"), prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), class = "sd"), prior(student_t(3, 0, 2.5), class = "sigma"),
  prior(student_t(3, 0, 2.5), class = "Intercept", dpar = "hu"))
ctrl <- list(adapt_delta = 0.97, max_treedepth = 12)

prep_glmm_data <- function(dat) {
  dat |> mutate(Year = as.factor(Year), Exposure = as.factor(Exposure),
                SiteDepthTransect = interaction(Site, Depth, Transect, sep = "-", drop = TRUE))
}

fit_one_glmm <- function(resp_col, data, seed = SEED_MASTER) {
  message("Fitting: ", resp_col)
  df <- data |>
    transmute(Year, Exposure, SiteDepthTransect, Response = as.numeric(.data[[resp_col]])) |>
    filter(!is.na(Response), !is.na(Year), !is.na(Exposure), !is.na(SiteDepthTransect))
  if (nrow(df) == 0)         { warning("Skipping ", resp_col, ": no complete rows."); return(NULL) }
  if (all(df$Response == 0)) { warning("Skipping ", resp_col, ": all values are zero."); return(NULL) }
  n_zero <- sum(df$Response == 0); prop_zero <- mean(df$Response == 0)
  message("  zeros = ", n_zero, " (", round(100 * prop_zero, 1), "%)")

  if (n_zero == 0) {
    fit <- brm(Response ~ Year * Exposure + (1 | SiteDepthTransect), data = df,
               family = Gamma(link = "log"), prior = priors_gamma,
               chains = 4, iter = 3000, warmup = 1000, cores = min(4, parallel::detectCores()),
               control = ctrl, seed = seed)
    fam <- "Gamma(log)"
  } else {
    fit <- brm(bf(Response ~ Year * Exposure + (1 | SiteDepthTransect), hu ~ 1), data = df,
               family = hurdle_lognormal(), prior = priors_hurdle_logn,
               chains = 4, iter = 3000, warmup = 1000, cores = min(4, parallel::detectCores()),
               control = ctrl, seed = seed)
    fam <- "hurdle_lognormal"
  }
  list(model = fit, family = fam, zeros = n_zero, prop_zero = prop_zero)
}

fit_all_glmms <- function(response_cols, data, model_dir, info_file, id_col = "Response") {
  fits <- purrr::imap(response_cols, ~ fit_one_glmm(.x, data, seed = 100 + .y))
  names(fits) <- response_cols
  fits <- fits[!vapply(fits, is.null, logical(1))]
  models <- purrr::map(fits, "model")
  dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)
  purrr::iwalk(models, ~ saveRDS(.x, file.path(model_dir, paste0(safe_file(.y), ".rds"))))
  info <- tibble(Response = names(fits), Family = purrr::map_chr(fits, "family"),
                 Zeros = purrr::map_int(fits, "zeros"), ZeroProp = purrr::map_dbl(fits, "prop_zero"))
  names(info)[1] <- id_col
  write.csv(info, info_file, row.names = FALSE)
  models
}

# ---- Population-level effects table -------------------------------------------
clean_term_label <- function(x) {
  x <- sub("^Intercept$", "Intercept", x)
  x <- sub("^Year([0-9]+)$", "Year: \\1", x)
  x <- sub("^Exposure(.*)$", "Exposure: \\1", x)
  x <- sub("^Year([0-9]+):Exposure(.*)$", "Year \\1 × \\2", x)
  x <- sub("^hu_Intercept$", "Zero hurdle: Intercept", x)
  x
}
row_order_key <- function(lbl) {
  if (lbl == "Intercept") return("A")
  if (grepl("^Year:", lbl)) return("B")
  if (grepl("^Exposure:", lbl)) return("C")
  if (grepl("^Year .* × ", lbl)) return("D")
  if (grepl("^Zero hurdle:", lbl)) return("E")
  "Z"
}

# Posterior mean, 95% CrI and exceedance probability (EP): Pr(b > 0) if mean > 0, else Pr(b < 0)
summarise_fixef_with_ep <- function(model, model_name, ci = 0.95) {
  alpha <- (1 - ci) / 2
  dr <- as_draws_df(model)
  fe_cols <- grep("^b_", names(dr), value = TRUE)
  if (!length(fe_cols)) stop("No fixed-effect draws found for ", model_name)
  dr[, fe_cols, drop = FALSE] |>
    pivot_longer(everything(), names_to = "term_raw", values_to = "value") |>
    group_by(term_raw) |>
    summarise(estimate = mean(value),
              conf.low  = quantile(value, alpha, names = FALSE),
              conf.high = quantile(value, 1 - alpha, names = FALSE),
              p_gt0 = mean(value > 0), p_lt0 = mean(value < 0), .groups = "drop") |>
    mutate(EP = ifelse(estimate >= 0, p_gt0, p_lt0), term = sub("^b_", "", term_raw), Model = model_name) |>
    select(Model, term, estimate, conf.low, conf.high, EP)
}

build_fixef_table <- function(models) {
  ft <- purrr::imap_dfr(models, ~ summarise_fixef_with_ep(.x, .y) |>
    mutate(term_label = clean_term_label(term),
           Result = paste0(sprintf("%.2f [%.2f, %.2f]", estimate, conf.low, conf.high), "; ",
                           sprintf("Pr(%s0)=%.2f", ifelse(estimate >= 0, ">", "<"), EP)),
           order_key = vapply(term_label, row_order_key, character(1))))
  wide <- ft |> select(term_label, Model, Result, order_key) |>
    arrange(order_key, term_label) |> distinct(term_label, Model, .keep_all = TRUE) |>
    select(-order_key) |> pivot_wider(names_from = Model, values_from = Result)
  r2 <- tibble(term_label = "Bayesian R²")
  for (nm in names(models)) r2[[nm]] <- sprintf("%.2f", bayes_R2(models[[nm]])[1])
  bind_rows(wide, r2)
}

# ---- Pairwise contrasts (Tables 2 / S2) ---------------------------------------
fmt_rr_ci_ep <- function(rr, lo, hi, ep) {
  paste0(sprintf("%.2f", rr), " [", sprintf("%.2f", lo), ", ", sprintf("%.2f", hi),
         "]; Pr(RR>1)=", sprintf("%.2f", ep))
}

# Posterior summary of each emmeans contrast on the link scale; RR = exp(difference); EP = Pr(diff > 0)
summarise_contrast_draws <- function(ct_emm, by_var = NULL, ci = 0.95) {
  alpha <- (1 - ci) / 2
  draws_df <- as.data.frame(do.call(rbind, as.mcmc(ct_emm)))
  lab <- broom::tidy(ct_emm) |> mutate(.col = colnames(draws_df))
  res <- purrr::imap_dfr(draws_df, ~ {
    v <- .x
    tibble(diff_mean = mean(v), diff_low = quantile(v, alpha, names = FALSE),
           diff_high = quantile(v, 1 - alpha, names = FALSE), EP = mean(v > 0),
           RR_mean = exp(mean(v)), RR_low = exp(quantile(v, alpha, names = FALSE)),
           RR_high = exp(quantile(v, 1 - alpha, names = FALSE)), .col = .y)
  }) |> left_join(lab, by = ".col")
  keep_by <- if (is.null(by_var)) character(0) else intersect(by_var, names(res))
  res |> select(all_of(keep_by), contrast, diff_mean, diff_low, diff_high, RR_mean, RR_low, RR_high, EP)
}

contrast_tables_for_model <- function(model, model_name, ci = 0.95) {
  minus <- function(x) gsub(" - ", " − ", x, fixed = TRUE)
  ct1 <- contrast(emmeans(model, ~ Exposure | Year, type = "link", re_formula = NA), method = "pairwise")
  t1 <- summarise_contrast_draws(ct1, "Year", ci) |>
    mutate(Scope = "Exposure within Year", Label = paste0(minus(contrast), " | Year ", Year),
           Result = fmt_rr_ci_ep(RR_mean, RR_low, RR_high, EP)) |> select(Scope, Label, Result)
  ct2 <- contrast(emmeans(model, ~ Year | Exposure, type = "link", re_formula = NA), method = "pairwise")
  t2 <- summarise_contrast_draws(ct2, "Exposure", ci) |>
    mutate(Scope = paste0("Year within Exposure (", Exposure, ")"),
           Label = paste0(minus(contrast), " | Exposure ", Exposure),
           Result = fmt_rr_ci_ep(RR_mean, RR_low, RR_high, EP)) |> select(Scope, Label, Result)
  bind_rows(t1, t2) |> mutate(Model = model_name)
}

build_contrast_table <- function(models) {
  purrr::imap_dfr(models, contrast_tables_for_model) |>
    distinct(Scope, Label, Model, Result) |> pivot_wider(names_from = Model, values_from = Result)
}

# ---- Estimated marginal means and plotting -------------------------------------
standardize_emm_df <- function(emm_obj) {
  df <- as.data.frame(emm_obj)
  est <- intersect(c("response", "emmean", "rate", "prob"), names(df))[1]
  lo  <- intersect(c("lower.HPD", "lower.CL", "asymp.LCL", "LCL"), names(df))[1]
  hi  <- intersect(c("upper.HPD", "upper.CL", "asymp.UCL", "UCL"), names(df))[1]
  if (anyNA(c(est, lo, hi))) stop("Could not standardize emmeans output.")
  df |> rename_with(~ "response", all_of(est)) |> rename_with(~ "lower", all_of(lo)) |> rename_with(~ "upper", all_of(hi))
}

# Population-level means (random effects excluded), response scale, 95% CrI.
get_emm_df <- function(model, label, level = 0.95) {
  emm <- emmeans(model, ~ Year * Exposure, type = "response", re_formula = NA, level = level)
  standardize_emm_df(emm) |> mutate(Model = label, Year = as.factor(Year), Exposure = as.factor(Exposure))
}

# compact = TRUE gives the smaller panels used for Fig S1
plot_emm <- function(emm_df, panel_title, y_lab, compact = FALSE) {
  s <- if (compact) 0.85 else 1
  ggplot(emm_df, aes(Year, response, colour = Exposure, group = Exposure)) +
    geom_point(position = position_dodge(width = 0.22), size = 2.4 * s) +
    geom_line(position = position_dodge(width = 0.22), linewidth = 0.7 * s) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.14, linewidth = 0.55 * s,
                  position = position_dodge(width = 0.22)) +
    scale_colour_manual(values = exposure_cols, drop = FALSE) +
    { if (compact) scale_y_continuous(labels = function(x) ifelse(x < 0, "", x)) } +  # Fig S1: hide labels < 0
    labs(title = panel_title, x = "Year", y = y_lab, colour = "Exposure") +
    theme_classic(base_size = if (compact) 10.5 else 12) +
    theme(plot.title = element_text(face = "bold", size = if (compact) 10.5 else 12, hjust = 0),
          axis.title.y = element_text(face = "bold", size = if (compact) 9.5 else 12),
          axis.text = element_text(colour = "black", size = if (compact) 8.5 else 10),
          legend.title = element_text(face = "bold"), plot.margin = margin(8, 8, 8, 8))
}

# ---- Condensed contrasts among main categories --------
# Instead of 48 contrasts per response (Exposure within each Year; Year within each Exposure), report
#   (a) Exposure contrasts averaged over the four survey years   (6 pairwise contrasts)
#   (b) Year contrasts averaged over the four exposure zones     (6 pairwise contrasts)
# Averaging is on the link (log) scale with equal weights, within each posterior draw, so the
# credible intervals and exceedance probabilities propagate the uncertainty correctly.
# RR = ratio of estimated means (first level / second level); EP = Pr(RR > 1).
main_effect_contrasts_one <- function(model, model_name, ci = 0.95) {
  clean <- function(x) gsub(" - ", " − ", gsub("Year|Exposure", "", x), fixed = TRUE)
  one <- function(by, scope) {
    ct <- contrast(emmeans(model, as.formula(paste("~", by)), type = "link", re_formula = NA), method = "pairwise")
    summarise_contrast_draws(ct, NULL, ci) |>
      transmute(Model = model_name, Scope = scope, Contrast = clean(contrast),
                RR = RR_mean, RR_low, RR_high, EP,
                Excludes_1 = RR_low > 1 | RR_high < 1)
  }
  bind_rows(one("Exposure", "Exposure (averaged over years)"), one("Year", "Year (averaged over exposures)"))
}

build_main_effect_long <- function(models, ci = 0.95) {
  purrr::imap_dfr(models, \(m, nm) main_effect_contrasts_one(m, nm, ci))
}

# Wide, manuscript-ready layout: one row per contrast, one column per response.
# `labels` = named vector, e.g. c(Total.density = "Total density"); order of labels sets column order.
widen_main_effect <- function(long, labels = NULL) {
  if (is.null(labels)) labels <- setNames(unique(long$Model), unique(long$Model))
  long |>
    mutate(Cell = fmt_rr_ci_ep(RR, RR_low, RR_high, EP), Model = factor(Model, levels = names(labels), labels = labels)) |>
    select(Scope, Contrast, Model, Cell) |>
    pivot_wider(names_from = Model, values_from = Cell)
}

# Writes the wide table (.csv, and .docx with bold = 95% CrI excludes 1) and the numeric long table (.csv)
save_main_effect_table <- function(long, name, labels = NULL, caption = NULL) {
  wide <- widen_main_effect(long, labels)
  utils::write.csv(wide, out_path(paste0(name, ".csv")), row.names = FALSE)
  utils::write.csv(long, out_path(paste0(name, "_numeric.csv")), row.names = FALSE)
  if (requireNamespace("flextable", quietly = TRUE)) {
    if (is.null(labels)) labels <- setNames(unique(long$Model), unique(long$Model))
    ft <- flextable::autofit(flextable::flextable(as.data.frame(wide)))
    ft <- flextable::merge_v(ft, j = "Scope")
    for (nm in names(labels)) {
      rows <- long |> filter(Model == nm) |> mutate(r = match(paste(Scope, Contrast), paste(wide$Scope, wide$Contrast)))
      rows <- rows$r[rows$Excludes_1]
      if (length(rows)) ft <- flextable::bold(ft, i = rows, j = labels[[nm]], bold = TRUE)
    }
    if (!is.null(caption)) ft <- flextable::set_caption(ft, caption)
    flextable::save_as_docx(ft, path = out_path(paste0(name, ".docx")))
  }
  invisible(wide)
}

MAIN_EFFECT_CAPTION <- paste(
  "Pairwise contrasts among exposure zones (averaged over survey years) and among survey years (averaged",
  "over exposure zones). Values are rate ratios (RR) of estimated means, first level / second level,",
  "with 95% credible intervals [CrI] and exceedance probability Pr(RR>1). Averages are on the log scale with",
  "equal weights. Bold: 95% CrI excludes 1.")
