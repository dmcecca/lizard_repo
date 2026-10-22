# =============================================================================
# helpers_brt.R
# Boosted regression tree (BRT) functions shared by scripts 04, 05 and 09.
#   fit_brt()             one dismo::gbm.step fit
#   extract_rel_inf()     relative influence of every predictor
#   bootstrap_rel_inf()   bootstrap mean and 95% CI of relative influence
# =============================================================================
require_pkgs(c("gbm", "dismo"))
suppressPackageStartupMessages({ library(gbm); library(dismo) })

BRT_PARAMS <- list(tree.complexity = 3, learning.rate = 0.002, bag.fraction = 0.75, max.trees = 2000)
N_BOOT <- 100     # bootstrap refits per model (more refits give more stable CIs but take longer)

# One gbm.step fit; NULL on failure
fit_brt <- function(df, response, predictors, family, verbose = FALSE) {
  tryCatch(
    dismo::gbm.step(data = df, gbm.x = match(predictors, names(df)), gbm.y = match(response, names(df)),
                    family = family, tree.complexity = BRT_PARAMS$tree.complexity,
                    learning.rate = BRT_PARAMS$learning.rate, bag.fraction = BRT_PARAMS$bag.fraction,
                    max.trees = BRT_PARAMS$max.trees, plot.main = FALSE, verbose = verbose),
    error = function(e) { message("  gbm.step failed: ", conditionMessage(e)); NULL })
}

# Relative influence with every predictor listed (0 if never used in a split)
extract_rel_inf <- function(brt_model, predictors) {
  tmp <- summary(brt_model, plotit = FALSE)
  out <- data.frame(var = predictors, rel.inf = 0, stringsAsFactors = FALSE)
  m <- match(out$var, as.character(tmp$var))
  out$rel.inf[!is.na(m)] <- as.numeric(tmp$rel.inf[m[!is.na(m)]])
  out
}

boot_once <- function(b, df, response, predictors, family) {
  df_b <- df[sample.int(nrow(df), replace = TRUE), , drop = FALSE]
  if ("Exposure" %in% names(df_b)) df_b$Exposure <- factor(df_b$Exposure, levels = levels(df$Exposure))
  fit_b <- fit_brt(df_b, response, predictors, family)
  if (is.null(fit_b)) return(NULL)
  out <- extract_rel_inf(fit_b, predictors)
  names(out)[2] <- "rel.inf.boot"
  out
}

# Bootstrap mean and 95% percentile CI of relative influence.
# n_cores > 1: parallel cluster on Windows, forked processes elsewhere.
bootstrap_rel_inf <- function(df, response, predictors, family, n_boot = N_BOOT, n_cores = 1) {
  empty <- data.frame(var = predictors, boot_mean = NA_real_, lower = NA_real_, upper = NA_real_, n_success = 0)
  if (nrow(df) < 2) return(empty)

  if (n_cores > 1 && .Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(n_cores); on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterSetRNGStream(cl, 123)
    parallel::clusterExport(cl, c("fit_brt", "extract_rel_inf", "boot_once", "BRT_PARAMS"), envir = globalenv())
    parallel::clusterEvalQ(cl, { library(dismo); library(gbm); NULL })
    boot_list <- parallel::parLapply(cl, seq_len(n_boot), boot_once, df = df, response = response,
                                     predictors = predictors, family = family)
  } else if (n_cores > 1) {
    RNGkind("L'Ecuyer-CMRG")
    boot_list <- parallel::mclapply(seq_len(n_boot), boot_once, df = df, response = response,
                                    predictors = predictors, family = family, mc.cores = n_cores)
  } else {
    boot_list <- lapply(seq_len(n_boot), boot_once, df = df, response = response,
                        predictors = predictors, family = family)
  }

  boot_tbl <- do.call(rbind, boot_list[!vapply(boot_list, function(x) is.null(x) || inherits(x, "try-error"), logical(1))])
  if (is.null(boot_tbl) || nrow(boot_tbl) == 0) return(empty)

  out <- empty
  for (i in seq_along(predictors)) {
    v <- boot_tbl$rel.inf.boot[boot_tbl$var == predictors[i]]
    if (length(v)) {
      out$boot_mean[i] <- mean(v, na.rm = TRUE)
      out$lower[i]     <- as.numeric(quantile(v, 0.025, na.rm = TRUE))
      out$upper[i]     <- as.numeric(quantile(v, 0.975, na.rm = TRUE))
      out$n_success[i] <- sum(is.finite(v))
    }
  }
  out
}
