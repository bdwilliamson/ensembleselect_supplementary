## average rank stat for SL + knockoffs
#' @param X the original covariates
#' @param X_k the knockoff covariates
#' @param y the outcome
#' @param cores the number of cores
#' @param ... additional arguments to SuperLearner
stat.sl_avg_rank <- function(X, X_k, y, cores = 1, ...) {
  parallel <- TRUE
  if (!requireNamespace("doMC", quietly = TRUE)) {
    warning("doMC is not installed. Without parallelization, the statistics will be slower to compute",
            call. = FALSE, immediate. = TRUE)
    parallel <- FALSE
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    warning("parallel is not installed. Without parallelization, the statistics will be slower to compute.",
            call. = FALSE, immediate. = TRUE)
    parallel <- FALSE
  }
  if (parallel) {
    ncores <- parallel::detectCores(all.tests = TRUE, logical = TRUE)
    if (cores == 2) {
      cores <- min(2, ncores)
    }
    else {
      if (cores > ncores) {
        warning(paste("The requested number of cores is not available. Using instead",
                      ncores, "cores"), immediate. = TRUE)
        cores <- ncores
      }
    }
    if (cores > 1) {
      doMC::registerDoMC(cores = cores)
      parallel <- TRUE
    }
    else {
      parallel <- FALSE
    }
  }
  swap <- rbinom(ncol(X), 1, 0.5)
  swap.M <- matrix(swap, nrow = nrow(X), ncol = length(swap),
                  byrow = TRUE)
  X.swap <- X * (1 - swap.M) + X_k * swap.M
  Xk.swap <- X * swap.M + X_k * (1 - swap.M)
  ## run the Super Learner
  x <- as.data.frame(cbind(X.swap, Xk.swap))
  if (any(duplicated(names(x)))) {
    names(x) <- paste0("V", 1:ncol(x))
  }
  mod <- SuperLearner::SuperLearner(Y = y, X = x, ...)
  p <- ncol(X)
  orig <- 1:p
  ## get the knockoff statistics
  if (sum(stringr::str_count(names(x), "V")) == 0) {
    feature_start <- "X"
  } else {
    feature_start <- "V"
  }
  all_import <- extract_import_sl(mod, x_nms = names(x), import_type = "all") %>%
    mutate(feature_num = as.numeric(stringr::str_remove(feature, feature_start))) %>%
    arrange(feature_num)
  Z_all <- all_import %>%
    mutate(wt_screened_rank = rank / colSums(mod$whichScreen)) %>%
    mutate(desc_rank = rank(-wt_screened_rank)) %>%
    arrange(feature_num)
  Z <- Z_all$desc_rank
  W <- (abs(Z[orig]) - abs(Z[orig + p])) * (1 - 2 * swap)
  return(W)
}
