## run the simulation once
#' @param indx the monte-carlo ID
#' @param n the sample size
#' @param p the number of features
#' @param family the type of outcome data to generate (e.g.,
#' continuous ["gaussian"] or binary ["binomial"])
#' @param link the link function
#' @param linear the type of linear predictor
#' (linear, e.g., X\beta; or nonlinear)
#' @param x_dist the type of x distribution
#' (i.e., "normal" for MVN or "nonnormal" for complex covariates)
#' @param rho the correlation between features
#' @param beta the beta_0 parameter
#' @param missing_perc the percentage of missing data
#' (maximum column-wise, on average)
#' @param missing_type if missing_perc > 0, then MAR nested or MAR complex
#' @param estimator_type the type of estimator to run
#' (i.e., "lasso", "SL", "SPVIM", or "SAGE+")
#' @param extra_layer the type of added layer (e.g., "SS" or "KF")
#' @param learners the candidate learners (a list, with screens) to pass to SL
#' @param spvim_learners the candidate learners to pass to SPVIM
#' @param uni_learners the candidate learners to pass to SPVIM
#' for univariate regressions
#' @param b the number of bootstrap replicates
#' @param M the number of multiple imputations
#' @param alpha the per-comparison error rate
#' @param thresh the threshold for SL variable selection
#' (without SS or knockoffs)
#' @param rank the rank-based threshold for SL variable selection
#' (without SS or knockoffs)
#' @param ss_pfer the SS parameter controlling the per-family error rate
#' (expected number of false positives)
#' @param ss_pi the SS parameter controlling the threshold (specify either this or \code{ss_q})
#' @param ss_q the SS parameter controlling the number of variables considered (specify either this
#'   or \code{ss_pi})
#' @param kf_fdr the knockoff false discovery rate
#' @param pfp_q the parameter controlling the proportion of false positives
#' (for SPVIM)
#' @param gfwer_k the parameter controlling the number of type I
#' errors we can make (for SPVIM)
do_one <- function(indx = 1, n = 100, p = 15, family = "binomial", link = NULL,
                   linear = TRUE, x_dist = "normal", rho = 0, beta = rep(0, p),
                   missing_perc = 0, missing_type = "nested",
                   estimator_type = "SL", extra_layer = "none",
                   learners = c("SL.ranger", "SL.glmnet", "SL.mean"),
                   spvim_learners = c("SL.ranger", "SL.glmnet", "SL.mean"),
                   uni_learners = "SL.polymars", b = 100, M = 1,
                   alpha = 0.05, thresh = 0.2, rank = 7,
                   ss_pfer = p * alpha, ss_pi = 0.75, ss_q = 7,
                   kf_fdr = 0.5,
                   pfp_q = 0.2, gfwer_k = 5, data_only = FALSE) {
  # generate data
  if (linear) {
    data_gen <- gen_data_lm
  } else {
    data_gen <- gen_data_nlm
  }
  if (nchar(link) > 0) {
    family <- paste0(family, "-", link)
  }
  dat <- data_gen(n, p, rho, matrix(beta), family, x_dist)
  test_dat <- data_gen(5e4, p, rho, matrix(beta), family, x_dist)
  full_dat <- dat
  # if the percentage of missing data > 0, then make data missing
  if (missing_perc > 0) {
    missing_pattern <- make_missing_pattern(type = missing_type,
                                            p = p)
    # make data missing; may need to run this several times
    ampute_error <- TRUE
    while(ampute_error) {
      dat <- try(expr = mice::ampute(full_dat, prop = missing_perc,
                                     patterns = missing_pattern, mech = "MAR",
                          freq = mice::ampute.default.freq(missing_pattern),
                          weights = mice::ampute.default.weights(missing_pattern,
                                                                 "MAR"),
                          std = FALSE, bycases = TRUE) %>%
                   magrittr::use_series("amp") %>%
                   tibble::as_tibble(),
                 silent = TRUE)
      if (class(dat)[1] != "try-error") {
        ampute_error <- FALSE
      }
    }
    # multiply impute
    mice_preds <- mice::quickpred(data = dat, include = "y", mincor = 0.25)
    mice_impute <- mice::mice(data = dat, m = M, method = "pmm",
                              predictorMatrix = mice_preds,
                              where = is.na(dat), printFlag = FALSE)
    # create a list of lists: each element of imputed data is a list
    # with x and y
    imputed_data <- mice::complete(mice_impute, action = "all")
  } else { # all deltas are zero
    # no need to multiply impute
    full_dat <- dat
    imputed_data <- list(dat)
  }
  # get true performance
  true_perf <- get_true_performance(test_dat, linear, family, beta, x_dist)
  if (data_only) {
    return(true_perf)
  }
  # set up return tibble
  selected_sets <- NULL
  est_vim <- NULL
  fam_for_ests <- switch((grepl("binomial", family)) + 1,
                         "gaussian", "binomial")
  for (m in 1:length(imputed_data)) {
    # fit the estimator
    if (!grepl("base", estimator_type)) {
      if (extra_layer == "none") {
        est_lst <- fit_extra_layer_none(
          est_type = estimator_type,
          dat = imputed_data[[m]],
          family = fam_for_ests, learners = learners,
          spvim_learners = spvim_learners,
          uni_learners = uni_learners,
          rank = rank, pfp_q = pfp_q,
          alpha = alpha,
          gfwer_k = gfwer_k, p = p
        )
      } else if (extra_layer == "SS") {
        est_lst <- fit_extra_layer_SS(
          est_type = estimator_type,
          dat = imputed_data[[m]],
          family = fam_for_ests, learners = learners,
          b = b, p = p, pfer = ss_pfer,
          pi = ss_pi, q = ss_q
        )
      } else if (extra_layer == "KF") {
        est_lst <- fit_extra_layer_KF(
          est_type = estimator_type,
          dat = imputed_data[[m]],
          family = fam_for_ests, kf_fdr = kf_fdr
        )
      } else {
        stop("The requested extra layer is not currently implemented.
             Please use either 'none', 'SS', or 'KF'.")
      }
      this_vim <- est_lst$vim
      selected_set <- est_lst$selected_set
      fit_out <- est_lst$fit_out
    } else {
      this_vim <- tibble::tibble(feature = paste0("V", 1:p), est = NA,
                                 rank = NA, feature_num = 1:p)
      selected_set <- rep(1, p)
      fit_out <- ""
    }
    # get performance of selected set (evaluated once)
    x <- imputed_data[[m]] %>% select(-y)
    y <- imputed_data[[m]] %>% pull(y)
    test_x <- test_dat %>% select(-y)
    test_y <- test_dat %>% pull(y)
    no_screen_lib <- unlist(lapply(learners,
                                   function(x) x[!grepl("screen", x) &
                                                   !grepl("All", x)]))
    if (!is.list(selected_set)) {
      selected_mod_perf <- get_selected_set_perf(
        selected_set = selected_set, x = x, y = y,
        test_x = test_x, test_y = test_y,
        estimator_type = estimator_type,
        learner_lib = no_screen_lib, family = fam_for_ests
      )
      est_perf <- selected_mod_perf$perf
      final_fit_out <- selected_mod_perf$fit_out
      selected_var_str <- paste(which(selected_set == 1), collapse = "; ")
      sens <- get_sensitivity(selected_var_str, 1:6)
      spec <- get_specificity(selected_var_str, 1:6, p)
    } else {
      # it's SPVIM; looking at multiple selected sets
      selected_mod_perf <- sapply(
        1:length(selected_set),
        function(s) {
          get_selected_set_perf(
            selected_set = selected_set[[s]], x = x, y = y, test_x = test_x,
            test_y = test_y, estimator_type = estimator_type,
            learner_lib = no_screen_lib, family = fam_for_ests
          )
          }, simplify = FALSE)
      est_perf <- lapply(selected_mod_perf, function(l) l$perf)
      final_fit_out <- lapply(selected_mod_perf, function(l) l$fit_out)
      names(selected_mod_perf) <- names(selected_set)
      selected_var_str <- lapply(
        selected_set, function(z) paste(which(z == 1), collapse = "; ")
        )
      sens <- unlist(lapply(
        selected_var_str, get_sensitivity, true_active_set = 1:6
      ))
      spec <- unlist(lapply(
        selected_var_str, get_specificity, true_active_set = 1:6, p = p
      ))
      selected_var_str <- unlist(selected_var_str)
    }
    miss_percs <- colMeans(is.na(dat))
    miss_summ <- summary(miss_percs)
    miss_summ_text <- sprintf(
      "Min. = %1$.3f; 1st Qu. = %2$.3f; Median = %3$.3f; Mean = %4$.3f; 3rd Qu. = %5$.3f; Max. = %6$.3f.",
                              miss_summ[1], miss_summ[2], miss_summ[3],
                              miss_summ[4], miss_summ[5], miss_summ[6]
      )
    this_est_type <- switch((is.list(selected_set)) + 1,
                            estimator_type,
                            paste(estimator_type, names(selected_set), sep = "-"))
    set_tib <- tibble::tibble(estimator_type = this_est_type,
                              extra_layer = extra_layer, m = m,
                              selected_vars = selected_var_str,
                              perf = unlist(est_perf),
                              true_perf = true_perf,
                              sensitivity = sens, specificity = spec,
                              overall_miss = mean(miss_percs),
                              miss_summ = miss_summ_text,
                              fit_out = fit_out,
                              final_fit_out = final_fit_out)
    selected_sets <- dplyr::bind_rows(selected_sets, set_tib)
    est_vim <- dplyr::bind_rows(est_vim, this_vim %>%
                                  mutate(estimator_type = estimator_type,
                                         extra_layer = extra_layer, m = m,
                                         .before = "feature"))
  }
  set_output <- selected_sets %>%
    tibble::add_column(mc_id = indx, family = family, linear = linear,
                       x_dist = x_dist, missing_perc = missing_perc,
                       missing_type = missing_type, b = b, m_tot = M,
                       thresh = thresh, rank_thr = rank, n = n, p = p,
                       .before = "estimator_type")
  vim_output <- est_vim %>%
    tibble::add_column(mc_id = indx, family = family, linear = linear,
                       x_dist = x_dist, missing_perc = missing_perc,
                       missing_type = missing_type, b = b, m_tot = M,
                       thresh = thresh, rank_thr = rank, n = n, p = p,
                       .before = "estimator_type")
  list(selected = set_output, vim = vim_output)
}
