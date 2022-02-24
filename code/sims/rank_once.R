## run the simulation once
#' @param n the sample size
#' @param p the number of features
#' @param family the type of outcome data to generate (e.g., continuous ["gaussian"] or binary ["binomial"])
#' @param linear the type of linear predictor (linear, e.g., X\beta; or nonlinear)
#' @param x_dist the type of x distribution (i.e., "normal" for MVN or "nonnormal" for complex covariates)
#' @param rho the correlation between features
#' @param beta the beta_0 parameter
#' @param estimator_type the type of estimator to run ("SL")
#' @param extra_layer the type of added layer (i.e., "none" or "SS")
#' @param selection_type the type of variable selection to run (i.e., 'ranks', 'BH', 'KF')
#' @param learners the candidate learners (a list, with screens) to pass to SL
#' @param uni_learners the candidate learners to pass to univariate regressions in SPVIM
#' @param b the number of bootstrap replicates
#' @param thresh the threshold for SL variable selection (without SS or knockoffs)
#' @param rank the rank-based threshold for SL variable selection (without SS or knockoffs)
#' @param q the SS parameter controlling how many features can make it into a model
#' @param fdr the knockoff false discovery rate
#' @param offset the knockoff offset (1 is knockoff+, 0 is not)
#' @param knockoff_method the knockoff method ('asdp', 'equi', 'sdp')
do_one <- function(n, p, family, linear, x_dist, rho, beta,
                   estimator_type, selection_type, extra_layer,
                   learners, uni_learners, b, thresh, rank, q, fdr = 0.5, offset = 1, knockoff_method = 'asdp') {
  # generate data
  if (linear) {
    data_gen <- gen_data_lm
  } else {
    data_gen <- gen_data_nlm
  }
  dat <- data_gen(n, p, rho, matrix(beta), family, x_dist)
  test_dat <- data_gen(10000, p, rho, matrix(beta), family, x_dist)
  dat$delta <- matrix(1, nrow = n, ncol = p)

  # set up return tibble
  selected_sets <- NULL
  # get true performance
  if (linear) {
    linear_predictor <- cbind(1, as.matrix(test_dat$x))%*%beta
  } else {
    linear_predictor <- nl_conditional_mean(test_dat$x, beta)
  }
  if (family == "binomial") {
    true_perf <- vimp::measure_auc(linear_predictor, test_dat$y)$point_est
  } else {
    true_perf <- vimp::measure_r_squared(linear_predictor, test_dat$y)$point_est
  }
  sl_opts <- get_sl_opts(family)
  no_screen_lib <- unlist(lapply(learners, function(x) x[!grepl("screen", x) & !grepl("All", x)]))
  no_screen_no_glmnet_lib <- no_screen_lib[!grepl("glmnet", no_screen_lib)]
  if (selection_type == "sl_rank") {
    if (extra_layer == "none") {
      mod_set_lst <- get_sl_rank_mod_plus_set(dat, learners, V = 5, sl_opts, rank, p)
      mod <- mod_set_lst$mod
      selected_set <- mod_set_lst$set
    } else if (extra_layer == "KF") {
      second_order_knockoffs <- knockoff::create.second_order(as.matrix(dat$x), method = knockoff_method, shrink = TRUE)
      colnames(second_order_knockoffs) <- paste0("V", (p+1):(2*p))
      dat2 <- dat
      dat2$x <- cbind(dat$x, second_order_knockoffs)
      mod_set_lst <- get_sl_rank_mod_plus_set(dat2, learners, V = 5, sl_opts, rank, p)
      mod <- mod_set_lst$mod
      if (sum(stringr::str_count(names(dat2$x), "V")) == 0) {
        feature_start <- "X"
      } else {
        feature_start <- "V"
      }
      all_import <- extract_import_sl(mod, x_nms = names(dat2$x), import_type = "all") %>%
        mutate(feature_num = as.numeric(stringr::str_remove(feature, feature_start))) %>%
        arrange(feature_num)
      Z_all <- all_import %>%
        mutate(wt_screened_rank = rank / colSums(mod$whichScreen)) %>%
        mutate(desc_rank = rank(-wt_screened_rank)) %>%
        arrange(feature_num)
      Z <- Z_all$desc_rank
      kf_W <- (abs(Z[1:p]) - abs(Z[(1:p) + p]))
      kf_thresh <- knockoff::knockoff.threshold(kf_W, fdr = fdr, offset = offset)
      selected_set <- as.numeric(kf_W >= kf_thresh)
    } else if (extra_layer == "SS") {
      # run SL + stability selection
      mod <- stabs::stabsel(x = dat$x, y = dat$y,
                            fitfun = flevr::SL_stabs_fitfun,
                            args.fitfun = list(screens = FALSE,
                                               SL.library = learners, cvControl = list(V = 5),
                                               family = sl_opts$fam, method = sl_opts$method),
                            assumption = "none", verbose = FALSE, eval = TRUE, papply = mclapply,
                            mc.preschedule = FALSE, mc.cores = 1,
                            cutoff = 0.7, B = b/2, sampling.type = "SS", PFER = 4)
      selected_set <- rep(0, p)
      selected_set[mod$selected] <- 1
    }
  } else if (selection_type == "spvim_rank") {
    if (extra_layer == "none") {
      mod_set_lst <- get_spvim_mod_plus_set(dat, learners, uni_learners, V = 2, sl_opts, rank, fdr, p, procedure = "rank")
      mod <- mod_set_lst$mod
      selected_set <- mod_set_lst$set
    } else if (extra_layer == "KF") {
      second_order_knockoffs <- knockoff::create.second_order(as.matrix(dat$x), method = knockoff_method, shrink = TRUE)
      colnames(second_order_knockoffs) <- paste0("V", (p+1):(2*p))
      dat2 <- dat
      dat2$x <- cbind(dat$x, second_order_knockoffs)
      mod_set_lst <- get_spvim_mod_plus_set(dat2, learners, uni_learners, V = 2, sl_opts, rank, fdr, p, procedure = "rank")
      mod <- mod_set_lst$mod
      kf_W <- abs(mod_set_lst$ests$est) - abs(mod_set_lst$ests$est[(1:p) + p])
      kf_thresh <- knockoff::knockoff.threshold(kf_W, fdr = fdr, offset = 1)
      selected_set <- as.numeric(kf_W >= kf_thresh)
    } else if (extra_layer == "SS") {
      stop("Extra layer 'SS' is not currently implemented for SPVIMs.")
    }
  } else if (selection_type == "spvim_BH") {
    if (extra_layer == "none") {
      mod_set_lst <- get_spvim_mod_plus_set(dat, learners, uni_learners, V = 2, sl_opts, rank, fdr, p, procedure = "BH")
      selected_set <- mod_set_lst$set
    } else if (extra_layer == "KF") {
      stop("Extra layer 'KF' is not currently implemented for procedure 'BH'.")
    } else if (extra_layer == "SS") {
      stop("Extra layer 'SS' is not currently implemented for SPVIMs.")
    }
  } else if (selection_type == "sage_rank") {
    if (extra_layer == "none") {
      mod_set_lst <- get_sage_mod_plus_set(dat, learners, uni_learners, V = 2, sl_opts, rank, fdr, p, procedure = "rank")
      mod <- mod_set_lst$mod
      selected_set <- mod_set_lst$set
    } else if (extra_layer == "KF") {
      second_order_knockoffs <- knockoff::create.second_order(as.matrix(dat$x), method = knockoff_method, shrink = TRUE)
      colnames(second_order_knockoffs) <- paste0("V", (p+1):(2*p))
      dat2 <- dat
      dat2$x <- cbind(dat$x, second_order_knockoffs)
      mod_set_lst <- get_sage_mod_plus_set(dat2, learners, uni_learners, V = 2, sl_opts, rank, fdr, p, procedure = "rank")
      mod <- mod_set_lst$mod
      kf_W <- abs(mod_set_lst$ests$est) - abs(mod_set_lst$ests$est[(1:p) + p])
      kf_thresh <- knockoff::knockoff.threshold(kf_W, fdr = fdr, offset = 1)
      selected_set <- as.numeric(kf_W >= kf_thresh)
    } else if (extra_layer == "SS") {
      stop("Extra layer 'SS' is not currently implemented for SAGEs.")
    }
  } else if (selection_type == "sl_none") {
    selected_set <- rep(1, p)
  } else if (selection_type == "sl_screen") {
    if (extra_layer == "none") {
      mod_set_lst <- get_sl_rank_mod_plus_set(dat, learners, V = 5, sl_opts, rank, p)
      mod <- mod_set_lst$mod
      selected_set <- mod_set_lst$set
    } else if (extra_layer == "KF") {
      second_order_knockoffs <- knockoff::create.second_order(as.matrix(dat$x), method = knockoff_method, shrink = TRUE)
      colnames(second_order_knockoffs) <- paste0("V", (p+1):(2*p))
      dat2 <- dat
      dat2$x <- cbind(dat$x, second_order_knockoffs)
      mod_set_lst <- get_sl_rank_mod_plus_set(dat2, learners, V = 5, sl_opts, rank, p)
      mod <- mod_set_lst$mod
      if (sum(stringr::str_count(names(dat2$x), "V")) == 0) {
        feature_start <- "X"
      } else {
        feature_start <- "V"
      }
      all_import <- extract_import_sl(mod, x_nms = names(dat2$x), import_type = "all") %>%
        mutate(feature_num = as.numeric(stringr::str_remove(feature, feature_start))) %>%
        arrange(feature_num)
      Z_all <- all_import %>%
        mutate(wt_screened_rank = rank / colSums(mod$whichScreen)) %>%
        mutate(desc_rank = rank(-wt_screened_rank)) %>%
        arrange(feature_num)
      Z <- Z_all$desc_rank
      kf_W <- (abs(Z[1:p]) - abs(Z[(1:p) + p]))
      kf_thresh <- knockoff::knockoff.threshold(kf_W, fdr = fdr, offset = offset)
      selected_set <- as.numeric(kf_W >= kf_thresh)
    } else if (extra_layer == "SS") {
      # run SL + stability selection
      mod <- stabs::stabsel(x = dat$x, y = dat$y,
                            fitfun = sl_stabs_fitfun,
                            args.fitfun = list(screens = TRUE,
                                               SL.library = learners, cvControl = list(V = 5),
                                               family = sl_opts$fam, method = sl_opts$method),
                            assumption = "none", verbose = FALSE, eval = TRUE, papply = mclapply,
                            mc.preschedule = FALSE, mc.cores = 1,
                            cutoff = 0.7, B = b/2, sampling.type = "SS", PFER = 4)
      selected_set <- rep(0, p)
      selected_set[mod$selected] <- 1
    }
  } else {
    stop("The requested variable selection type is not currently supported.")
  }

  # get performance of selected set (evaluated once)
  if (sum(selected_set) == 1) {
    # add a column of zeros to trick glmnet
    selected_vars <- cbind(0, dat$x[, selected_set == 1, drop = FALSE])
    test_selected_vars <- cbind(0, test_dat$x[, selected_set == 1, drop = FALSE])
  } else {
    selected_vars <- dat$x[, selected_set == 1, drop = FALSE]
    test_selected_vars <- test_dat$x[, selected_set == 1, drop = FALSE]
  }
  if (sum(selected_set) == 0) {
    selected_mod_preds <- rep(mean(test_dat$y), length(test_dat$y))
  } else {
    if (!grepl("screen", args$extra_layer)) {
      selected_mod <- SuperLearner::SuperLearner(Y = dat$y, X = as.data.frame(selected_vars),
                                                 SL.library = no_screen_lib, cvControl = list(V = 5), family = sl_opts$fam, method = sl_opts$method)
    } else {
      selected_mod <- SuperLearner::SuperLearner(Y = dat$y, X = as.data.frame(selected_vars),
                                                 SL.library = learners, cvControl = list(V = 5), family = sl_opts$fam, method = sl_opts$method)
    }
    selected_mod_preds <- predict(selected_mod, newdata = as.data.frame(test_selected_vars))$pred
  }
  if (family == "binomial") {
    selected_mod_perf <- vimp::measure_auc(selected_mod_preds, test_dat$y)$point_est
  } else {
    selected_mod_perf <- vimp::measure_r_squared(selected_mod_preds, test_dat$y)$point_est
  }
  set_tib <- tibble::tibble(estimator_type = estimator_type,
                            extra_layer = extra_layer,
                            selection_type = selection_type,
                            selected_vars = paste(which(selected_set == 1),
                                                  collapse = "; "),
                            perf = selected_mod_perf, true_perf = true_perf)
  selected_sets <- dplyr::bind_rows(selected_sets, set_tib)
  output <- selected_sets %>%
    tibble::add_column(family = family, linear = linear, x_dist = x_dist,
                       b = b, thresh = thresh, rank = rank, n = n, p = p,
                       .before = "estimator_type")
  return(output)
}
