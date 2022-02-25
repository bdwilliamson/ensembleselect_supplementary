# all functions for missing data imputation

# impute training and testing datasets -----------------------------------------
# @param dataset the dataset with missing data
# @param preimputed_data pre-imputed data, only used if covariates & biomarkers = TRUE
# @param mc_id the current monte-carlo id, only used if covariates & biomarkers = TRUE
# @param cross_fitting_folds the folds for cross-fitting (train/test split)
# @param k the current hold-out fold
# @param M the number of imputations
# @param covariates true/false, whether to use covariates
# @param biomarkers true/false, whether to use biomarkers
# @param no_impute_vars variables to not use in imputations (ID variables)
# @param use_cea true/false, should we use CEA in imputation model?
impute_train_and_test <- function(dataset = NULL, preimputed_data = NULL,
                                  mc_id = 1, cross_fitting = TRUE,
                                  cross_fitting_folds = NULL,
                                  k = 1, M = 1, covariates = FALSE,
                                  biomarkers = TRUE, no_impute_vars = NULL, use_cea = FALSE) {
  # if not the covariate-only or biomarker-only analysis, merge biomarker
  # data with pre-imputed covariate data
  if (covariates & biomarkers) {
    these_preimputed_data <- preimputed_data %>%
      filter(mc_id == !!mc_id, k == !!k)
    covariate_nms <- names(these_preimputed_data %>% select(-mc_id, -k, -imp, -id))
    biomarkers_and_imputed_covars <- lapply(as.list(unique(these_preimputed_data$imp)), function(i) {
      left_join(dataset %>% select(-covariate_nms),
                these_preimputed_data %>% filter(imp == !!i), by = "id")
    })
    # run a single imputation for each set of imputed covariates
    full_imputed_training_data <- as_tibble(rbindlist(lapply(as.list(seq_len(M)), function(m) {
      these_data <- biomarkers_and_imputed_covars[[m]] %>%
        filter(imp == !!m, cross_fitting_folds != !!k)
      imputed_init <- impute_missing_data(dataset = these_data %>% select(-imp), M = 1,
                                          covariates_only = FALSE,
                                          no_impute_vars = no_impute_vars,
                                          use_cea = use_cea)
      imputed_init %>% select(-imp) %>% mutate(imp = these_data$imp)
    })))
    full_imputed_testing_data <- as_tibble(rbindlist(lapply(as.list(seq_len(M)), function(m) {
      these_data <- biomarkers_and_imputed_covars[[m]] %>%
        filter(imp == !!m, cross_fitting_folds == !!k)
      imputed_init <- impute_missing_data(dataset = these_data %>% select(-imp), M = 1,
                                          covariates_only = FALSE,
                                          no_impute_vars = no_impute_vars,
                                          use_cea = use_cea)
      imputed_init %>% select(-imp) %>% mutate(imp = these_data$imp)
    })))
    imputed_training_data <- full_imputed_training_data %>%
      select(-mc_id, -k)
    imputed_testing_data <- full_imputed_testing_data %>%
      select(-mc_id, -k)
    all_imputed_data <- bind_rows(full_imputed_training_data, full_imputed_testing_data) %>%
      arrange(imp, id)
  } else { # run multiple imputation
    # split into test and train
    this_train <- dataset %>%
      filter(cross_fitting_folds != k)
    this_test <- dataset %>%
      filter(cross_fitting_folds == k)
    imputed_training_data <- impute_missing_data(dataset = this_train, M = M,
                                                 covariates_only = covariates & !biomarkers,
                                                 no_impute_vars = no_impute_vars,
                                                 use_cea = use_cea)
    imputed_testing_data <- impute_missing_data(dataset = this_test, M = M,
                                                covariates_only = covariates & !biomarkers,
                                                no_impute_vars = no_impute_vars,
                                                use_cea = use_cea)
    # save the imputed data to a list
    imputed_data <- bind_rows(imputed_training_data, imputed_testing_data) %>%
      arrange(imp, id)
    all_imputed_data <- imputed_data %>% mutate(mc_id = mc_id, k = k, .before = "imp")
  }
  return(list(train = imputed_training_data, test = imputed_testing_data,
              all = all_imputed_data))
}

# impute missing data ----------------------------------------------------------
# @param dataset the dataset with missing data
# @param M the number of imputations
# @param covariates_only is this a covariate-only analysis?
# @param no_impute_vars any variables not to impute
# @param impute_vars the variables to use for imputation
# @param method the imputation method to use
# @param max_iter the maximum number of iterations for imputation
# @param use_cea whether or not to use CEA to predict outcome
impute_missing_data <- function(dataset = NULL, M = 10, covariates_only = FALSE,
                                no_impute_vars = NULL, impute_vars = NULL,
                                method = "pmm", max_iter = 20, use_cea = TRUE) {
  dataset_to_impute <- dataset %>%
    select(-!!no_impute_vars)
  if (any(grepl("race", names(dataset_to_impute)))) {
    race_unknown <- is.na(dataset_to_impute$race_other_or_unknown)
    dataset_to_impute$race_other_or_unknown[race_unknown] <- 1
    dataset_to_impute[race_unknown, c("race_black", "race_american_indian_alaska_native",
                                      "race_asian", "race_native_hawaiian_pacific_islander")] <- 0
  }
  all_nas <- is.na(dataset_to_impute)
  init <- mice::mice(data = dataset_to_impute, m = M, maxit = 0, where = all_nas,
                     printFlag = FALSE, remove_collinear = FALSE, method = method)
  pred_matrix <- init$predictorMatrix
  all_methods <- init$method
  if (!covariates_only) {
    # set up passive imputation of binary calls
    pred_matrix["lab1_molecules_score", "lab1_molecules_neoplasia_call"] <- 0
    pred_matrix["lab1_telomerase_score", "lab1_telomerase_neoplasia_call"] <- 0
    pred_matrix["lab2_fluorescence_score", "lab2_fluorescence_mucinous_call"] <- 0
    pred_matrix["lab4_areg_score", "lab4_areg_mucinous_call"] <- 0
    pred_matrix["lab4_glucose_score", "lab4_glucose_mucinous_call"] <- 0
    pred_matrix[c("lab4_areg_score", "lab4_glucose_score"), "lab4_combined_mucinous_call"] <- 0
    pred_matrix["lab6_m_ab_das_score", "lab6_m_ab_das_neoplasia_call"] <- 0
    pred_matrix["cea", "cea_call"] <- 0
    all_methods["lab1_molecules_neoplasia_call"] <- "~ I(lab1_molecules_score > 25)" # confirmed
    all_methods["lab1_telomerase_neoplasia_call"] <- "~ I(lab1_telomerase_score > 730)" # confirmed
    all_methods["lab2_fluorescence_mucinous_call"] <- "~ I(lab2_fluorescence_score > 1.23)" # confirmed
    all_methods["lab4_areg_mucinous_call"] <- "~ I(lab4_areg_score > 112)" # confirmed
    all_methods["lab4_glucose_mucinous_call"] <- "~ I(lab4_glucose_score < 50)" # confirmed
    all_methods["lab4_combined_mucinous_call"] <- "~ I((lab4_areg_score > 100) * (lab4_glucose_score < 50))" # confirmed
    all_methods["lab6_m_ab_das_neoplasia_call"] <- "~ I(lab6_m_ab_das_score > 0.104)" # confirmed
    all_methods["cea_call"] <- "~ I(cea > 192)" # confirmed
    # if outcome is high malignancy, don't use CEA to predict
    if (!use_cea) {
      pred_matrix["high_malignancy", "cea"] <- 0
    }
  }
  # do imputation
  imputed_data <- tryCatch(mice::mice(data = dataset_to_impute, m = M, method = all_methods,
                                      where = all_nas, predictorMatrix = pred_matrix,
                                      printFlag = FALSE, maxit = max_iter),
                           error = function(e) {
                             pred_matrix_2 <- pred_matrix
                             pred_matrix_2[, colMeans(all_nas) > 0.5] <- 0
                             method2 <- all_methods
                             method2_2 <- all_methods
                             method2["cea"] <- "cart"
                             method2_2[method2_2 == "pmm"] <- "cart"
                             tryCatch(mice::mice(data = dataset_to_impute, m = M, method = method2,
                                                 where = all_nas, predictorMatrix = pred_matrix_2,
                                                 printFlag = FALSE, maxit = max_iter,
                                                 remove_collinear = FALSE),
                                      error = function(e) {
                                        mice::mice(data = dataset_to_impute, m = M, method = method2_2,
                                                   where = all_nas, predictorMatrix = pred_matrix_2,
                                                   printFlag = FALSE, maxit = max_iter,
                                                   remove_collinear = FALSE)
                             })
                           })
  complete_data <- tibble::as_tibble(mice::complete(imputed_data, action = "long")) %>%
    left_join(dataset %>%
                select(!!no_impute_vars) %>%
                mutate(.id = 1:nrow(dataset)), by = ".id") %>%
    select(-.id) %>%
    rename(imp = .imp)
  if (!covariates_only) {
    if (sum(is.na(complete_data$cea_call)) > 0) {
      complete_data <- complete_data %>%
        mutate(cea_call = as.numeric(cea > 192))
    }
  }
  complete_data
}
