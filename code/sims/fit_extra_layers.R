# fit the different extra layer types

# fit extra layer == "none"
# @param est_type the type of estimator
# @param dat the dataset (a list with x and y)
# @param family the family to consider (i.e., "binomial")
# @param learners the learners for SL
# @param uni_learners the univariate learners for SPVIM
# @param pfp_q the desired proportion of false positives
# @param alpha the per-comparison alpha
# @param gfwer_k the desired expected number of type I errors
# @param p the number of covariates
# @return the selected set
fit_extra_layer_none <- function(est_type = "lasso", dat = NULL, family = "binomial",
                                 learners = c("SL.ranger", "SL.glmnet", "SL.mean"),
                                 spvim_learners = c("SL.ranger", "SL.glmnet", "SL.mean"),
                                 uni_learners = "SL.polymars", rank = 7, pfp_q = 0.2,
                                 alpha = 0.1, gfwer_k = 5, p = 500) {
  opts <- get_sl_opts(family)
  if (est_type == "lasso") {
    # run vanilla CV lasso
    mod <- glmnet::cv.glmnet(x = as.matrix(dat %>% select(-y)),
                             y = dat %>% pull(y), nfolds = 10, intercept = TRUE,
                             family = opts$fam)
    vim_select_tib <- flevr::extract_importance_glmnet(fit = mod, feature_names = paste0("V", 1:p),
                                                       coef = 1) %>%
      mutate(feature_num = as.numeric(gsub("V", "", feature))) %>%
      select(feature, importance, rank, feature_num) %>%
      rename(est = importance) %>%
      arrange(feature_num)
    est_vim <- vim_select_tib
    selected_set <- as.numeric(vim_select_tib$est != 0)
  } else if (grepl("SL", est_type)) {
    if (grepl("base", est_type)) {
      # no selection
      mod <- NULL
      est_vim <- tibble::tibble(feature = paste0("V", 1:p),
                                est = NA, rank = NA, feature_num = 1:p)
      selected_set <- rep(1, p)
    } else {
      # run SL
      y <- dat %>% pull(y)
      x <- dat %>% select(-y)
      mod <- SuperLearner::SuperLearner(Y = y, X = x, SL.library = learners,
                                        cvControl = list(V = 5,
                                                         stratifyCV = (opts$fam$family == "binomial")),
                                        family = opts$fam, method = opts$method)
      vim_select_tib <- flevr::extrinsic_selection(fit = mod,
                                                   feature_names = paste0("V", 1:p),
                                                   threshold = rank, import_type = "all",
                                                   x = x, y = y) %>%
        mutate(feature_num = as.numeric(gsub("V", "", feature)),
               est = rank) %>%
        select(feature, est, rank, feature_num, selected) %>%
        arrange(feature_num)
      selected_set <- as.numeric(vim_select_tib$selected)
      est_vim <- vim_select_tib %>% select(-selected)
    }
  } else if (grepl("SPVIM", est_type)) {
    # set up
    q_1 <- pfp_q
    fdr_q <- 1 - sqrt(1 - q_1)
    control_lists <- list(
      list(list(quantity = "gFWER", base_method = "Holm", k = gfwer_k), alpha = alpha),
      list(list(quantity = "PFP", base_method = "Holm", k = gfwer_k, q = pfp_q), alpha = alpha),
      list(list(quantity = "FDR", base_method = "Holm", k = gfwer_k, q = fdr_q), alpha = fdr_q)
      # for FDR control at level 0.2, set alpha = q = 1 - sqrt(1 - 0.2)
    )
    this_v <- 2
    this_v_2 <- 2
    y <- dat %>% pull(y)
    x <- dat %>% select(-y)
    # estimate SPVIM
    gamma <- 1
    nrounds <- 10
    for (i in seq_len(nrounds)) {
      z_w_lst <- vimp::sample_subsets(p = ncol(x), n = nrow(x), gamma = gamma)
      if (nrow(z_w_lst$Z) > ncol(x) | (i == nrounds)) {
        break
      } else {
        gamma <- gamma + 0.5
      }
    }
    est_spvims <- vimp::sp_vim(Y = y, X = x, V = this_v, gamma = gamma,
                               type = ifelse(opts$fam$family == "gaussian", "r_squared", "auc"),
                               stratified = !(opts$fam$family == "gaussian"),
                               SL.library = spvim_learners, univariate_SL.library = uni_learners,
                               cvControl = list(V = this_v_2, stratifyCV = !(opts$fam$family == "gaussian")),
                               family = opts$fam, method = opts$method)
    mod <- est_spvims$preds_lst
    est_vim <- est_spvims$mat %>%
      mutate(feature = paste0("V", s), feature_num = as.numeric(s), rank = rank(-abs(est))) %>%
      select(feature, est, p_value, rank, feature_num)
    selected_lst <- lapply(control_lists, function(control) {
      flevr::intrinsic_selection(spvim_ests = est_spvims, sample_size = nrow(x),
                                 alpha = control[[2]], feature_names = paste0("V", 1:ncol(x)),
                                 control = control[[1]])
    })
    selected_set <- lapply(selected_lst, function(l) as.numeric(l$selected))
    names(selected_set) <- c("gFWER", "PFP", "FDR")
  } else {
    stop("The entered estimator is not currently supported for
         extra layer = 'none'. Please enter one of 'lasso', 'SL', or 'SPVIM'.")
  }
  fit_out <- get_fit_out(mod)
  list(vim = est_vim, selected_set = selected_set, fit_out = fit_out)
}

# fit extra layer == "SS"
# @param est_type the type of estimator
# @param dat the dataset (a list with x and y)
# @param family the family to consider (i.e., "binomial")
# @param learners the learners for SL
# @param b the number of replicates for SS
# @param p the number of covariates
# @param pfer the max per-family error rate
# @param pi the threshold for stable variables (specify either this or q)
# @param q the number of vars to consider (specify either this or pi)
# @return the selected set
fit_extra_layer_SS <- function(est_type = "lasso", dat = NULL,
                               family = "binomial",
                               learners = c("SL.ranger", "SL.glmnet",
                                            "SL.mean"),
                               b = 100, p = 500,
                               pfer = p * 0.05, pi = 0.75, q = 10) {
  x <- as.matrix(dat %>% select(-y))
  y <- dat %>% pull(y)
  opts <- get_sl_opts(family)
  if (opts$fam$family == "binomial") {
    n_0 <- sum(y == 0)
    n_1 <- sum(y == 1)
    k_0 <- floor(n_0 * 0.5)
    k_1 <- ceiling(n_1 * 0.5)
    indices_0 <- rep(c(0, 1), c(n_0 - k_0, k_0))
    indices_1 <- rep(c(0, 1), c(n_1 - k_1, k_1))
    folds_0 <- replicate(b / 2, sample(indices_0))[sample(1:n_0), , drop = FALSE]
    folds_1 <- replicate(b / 2, sample(indices_1))[sample(1:n_1), , drop = FALSE]
    folds <- matrix(0, nrow = nrow(x), ncol = b / 2)
    folds[y == 0, ] <- folds_0
    folds[y == 1, ] <- folds_1
  } else {
    folds <- stabs::subsample(rep(1, nrow(x)), B = b / 2)
  }
  if (est_type == "lasso") {
    # run CV lasso + stability selection
    mod <- stabs::stabsel(x = x, y = y,
                          fitfun = stabs::glmnet.lasso,
                          args.fitfun = list(family = opts$fam, type = "conservative"),
                          assumption = "none", verbose = FALSE,
                          eval = TRUE, papply = mclapply,
                          mc.preschedule = FALSE, mc.cores = 1,
                          # cutoff = pi, PFER = pfer,
                          q = q, PFER = pfer,
                          B = b/2, sampling.type = "SS", folds = folds)
    if ((mod$cutoff > 0.9)) {
      mod <- stabs::stabsel(mod, cutoff = 0.9)
    }
  } else if (est_type == "SL") {
    # run SL + stability selection
    mod <- stabs::stabsel(x = x, y = y,
                          fitfun = flevr::SL_stabs_fitfun,
                          args.fitfun = list(
                            SL.library = learners,
                            cvControl = list(V = 2, stratifyCV = (opts$fam$family == "binomial")),
                            family = opts$fam, method = opts$method
                          ),
                          assumption = "none", verbose = FALSE,
                          eval = TRUE, papply = mclapply,
                          mc.preschedule = FALSE, mc.cores = 1,
                          q = q, PFER = pfer,
                          B = b/2, sampling.type = "SS", folds = folds)
    if (mod$cutoff > 0.8) {
      mod <- stabs::stabsel(mod, cutoff = 0.8)
    }
  } else {
    stop("The entered estimator is not currently supported for
         extra layer = 'SS'. Please enter one of 'lasso', 'SL'.")
  }
  est_vim <- tibble::tibble(
    feature = names(mod$max), est = mod$max, rank = rank(-est),
    feature_num = 1:p
  )
  selected_set <- rep(0, p)
  selected_set[mod$selected] <- 1
  fit_out <- ""
  list(vim = est_vim, selected_set = selected_set, fit_out = fit_out)
}

# fit extra layer == "KF"
# @param est_type the type of estimator
# @param dat the dataset (a list with x and y)
# @param family the family to consider (i.e., "binomial")
# @param kf_fdr the desired fdr
# @return the selected set
fit_extra_layer_KF <- function(est_type = "lasso", dat = NULL,
                               family = "binomial",
                               kf_fdr = 0.5) {
  opts <- get_sl_opts(family)
  if (est_type == "lasso") {
    x <- as.matrix(dat %>% select(-y))
    y <- dat %>% pull(y)
    # run CV lasso + knockoffs; knockoff statistic is
    # difference in absolute value of lasso coefficient
    second_order_knockoffs <- knockoff::create.second_order(x, method = "asdp",
                                                            shrink = FALSE)
    kf_lasso_W <- knockoff::stat.glmnet_coefdiff(X = x,
                                                 X_k = second_order_knockoffs,
                                                 y = y, family = opts$fam,
                                                 cores = 1)
    kf_lasso_thresh <- knockoff::knockoff.threshold(kf_lasso_W,
                                                    fdr = kf_fdr, offset = 1)
    if (is.infinite(kf_lasso_thresh)) {
      kf_lasso_thresh <- knockoff::knockoff.threshold(kf_lasso_W,
                                                      fdr = kf_fdr, offset = 0)
    }
    p <- length(kf_lasso_W)
    est_vim <- tibble::tibble(
      feature = paste0("V", 1:p), est = kf_lasso_W, rank = rank(-kf_lasso_W),
      feature_num = 1:p
    )
    selected_set <- as.numeric(kf_lasso_W >= kf_lasso_thresh)
  } else {
    stop("The entered estimator is not currently supported for
         extra layer = 'KF'. Please enter 'lasso'.")
  }
  fit_out <- ""
  list(vim = est_vim, selected_set = selected_set, fit_out = fit_out)
}
