# run the simulation for a specified set of parameters

# load required functions and packages
library("mice")
library("boot")
library("tidyverse")
library("data.table")
library("SuperLearner")
library("glmnet")
library("xgboost")
library("ranger")
library("kernlab")
library("caret")
library("argparse")
library("knockoff")
library("parallel")
library("stabs")
library("nloptr")
library("cvAUC")
# if you need to update: remotes::install_github("bdwilliamson/flevr@main", upgrade = "never", lib = .libPaths()[2])
library("flevr")

# pull in job id, set up command-line arguments
if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
  # if you need to update on cluster: remotes::install_github("bdwilliamson/vimp", upgrade = "never", lib = .libPaths()[2])
  library("vimp")
  job_id <- 1001 # for n = 200, p = 500
  code_dir <- "sim_code/"
  prefix <- "sim_code/"
} else {
  # if you need to update on cluster: remotes::install_github("bdwilliamson/vimp", upgrade = "never", lib = .libPaths()[2])
  library("vimp", lib.loc = .libPaths()[2])
  # edit the following line if using a different system than Slurm
  job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  code_dir <- "./"
  # edit the following line to where you wish to save results
  prefix <- paste0("path_to_results")
}
source(paste0(code_dir, "select_once.R"))
source(paste0(code_dir, "gen_data.R"))
source(paste0(code_dir, "utils.R"))
source(paste0(code_dir, "get_true_performance.R"))
source(paste0(code_dir, "fit_extra_layers.R"))

parser <- ArgumentParser()
parser$add_argument("--sim-name", default = "binomial-probit-nonlinear-nonnormal-nested",
                    help = "the name of the simulation")
parser$add_argument("--nreps-total", type = "double", default = 1000,
                    help = "number of replicates in total")
parser$add_argument("--nreps-per-job", type = "double", default = 1,
                    help = "number of replicates for each job")
parser$add_argument("--b", type = "double", default = 100,
                    help = "number of bootstrap replicates")
parser$add_argument("--m", type = "double", default = 1,
                    help = "number of MI replicates")
est_help <- paste0("the estimator to assess. Possible values are: ",
                    "'lasso'; 'SL' or 'base-SL'; and 'SPVIM' or 'SAGE+'.")
parser$add_argument("--est-type", default = "SPVIM",
                    help = est_help)
layer_help <- paste0("extra layer (i.e., 'none' = none, ",
                    "'SS' = stability selection, 'KF' = knockoffs)")
parser$add_argument("--extra-layer", default = "none",
                    help = layer_help)
args <- parser$parse_args()
print(paste0("Running sim ", args$sim_name, " with ",
             as.character(args$nreps_total), " total replicates and ",
      as.character(args$nreps_per_job), " replicates per job."))
print(paste0("Running estimator ", args$est_type, " with extra layer ",
             args$extra_layer, "."))
print(paste0("Additional parameters: number of bootstrap replicates = ",
             args$b, "; number of MI replicates = ", args$m, "."))

# set up all of the possible simulation parameters
ps <- c(30, 500)
ns <- c(200, 500, 1500, 3000)
miss_perc <- c(0, 0.2, 0.4)
num_unique_n_p_miss <- length(unique(ps)) * length(unique(ns)) *
  length(unique(miss_perc))
families <- c("binomial")

# grab the current simulation parameters
nreps_per_combo <- args$nreps_total/args$nreps_per_job
param_grid <- expand.grid(mc_id = 1:nreps_per_combo, p = ps, n = ns,
                          miss = miss_perc)
current_dynamic_args <- param_grid[job_id, ]
print(paste0("N = ", current_dynamic_args$n, "; p = ", current_dynamic_args$p,
             "; missing proportion = ", current_dynamic_args$miss))

# note the intercept
constant <- ifelse(grepl("nonlinear", args$sim_name), 2, 1)
beta_active <- c(0.5, -1, 1, -0.5, 0.5, 1/3, -1/3) * constant
beta_nonactive <- 0
active_set <- 1:7 # really is the first six, since we use an intercept
non_active_set <- (1:current_dynamic_args$p)[-active_set]
beta_0 <- vector("numeric", current_dynamic_args$p + 1)
beta_0[active_set] <- beta_active
beta_0[non_active_set] <- beta_nonactive

# set up fixed parameters
per_comparison_alpha <- 0.05
sl_threshold <- 0.2
# match SL rank selection with SS q's; note q ~ sqrt(5(p / 2))
sl_rank <- ifelse(current_dynamic_args$p == 30, 9 + 1, 35 + 1)
kf_fdr <- 0.2
# stabsel: recommend alpha <= pfer_max <= p * alpha
# so at p = 30, for per-comparison alpha = 0.2, get 6 (i.e., p * PCER)
ss_pcer <- 0.04
ss_pfer <- current_dynamic_args$p * ss_pcer
ss_pi <- 0.9
ss_q <- ceiling(sqrt(ss_pfer * (2 * ss_pi - 1) * current_dynamic_args$p)) + 1
# since we know the true active set, we can specify a specificity of interest
# here, set target specificity of 0.85 at n = 3000 (p = 30), specificity 0.95 for p = 500.
# this implies that
# k = 4 (= 0.15 * 24) when p = 30; and
# k = 50 (= 0.1 * 494) when p = 500
target_specificity <- set_specificity(n = current_dynamic_args$n, ns = ns,
                                      p = current_dynamic_args$p,
                                      max_spec = ifelse(current_dynamic_args$p == 30, 0.85, 0.90))
gfwer_k <- ceiling((1 - target_specificity) * (current_dynamic_args$p - (length(active_set) - 1)))
# we can also specify a PFP of interest; want more control at higher sample size
# note that at n = 200, we want the maximum PFP
# pfp_q <- gfwer_k / ((current_dynamic_args$p - 6) / current_dynamic_args$p * (sqrt(current_dynamic_args$n) / sqrt(200)) + 6)
pfp_q <- gfwer_k / ((current_dynamic_args$p - 6) / current_dynamic_args$p * (sqrt(current_dynamic_args$n) / sqrt(200)) + gfwer_k)

# set up SuperLearner library
xgb_tune_params <- list(max_depth = c(1, 3), shrinkage = c(0.1))
xgb_learners <- create.Learner("SL.xgboost", tune = xgb_tune_params,
                               detailed_names = TRUE, name_prefix = "xgb")

if (grepl("screen", args$est_type)) {
  screen_lib <- c("screen.glmnet", "screen.corRank.p", "All")
  screen_plus_learner_lib <- sapply(
    1:length(learner_lib),
    function(i) c(learner_lib[i], screen_lib),
    simplify = FALSE
  )
} else {
  screen_plus_learner_lib <- learner_lib
}

# ---------------------------------------------
# replicate the simulation nreps_per_job times
# ---------------------------------------------
current_seed <- job_id
print(current_seed)
set.seed(current_seed)
system.time(sim_output <- sapply(1:args$nreps_per_job, function(i)
  do_one(indx = i + args$nreps_per_job * (current_dynamic_args$mc_id - 1),
         n = current_dynamic_args$n, p = current_dynamic_args$p,
         family = ifelse(grepl("binomial", args$sim_name),
                         "binomial", "gaussian"),
         link = ifelse(grepl("probit", args$sim_name), "probit", ""),
         linear = !grepl("nonlinear", args$sim_name),
         x_dist = ifelse(!grepl("nonnormal", args$sim_name),
                         "normal", "nonnormal"),
         rho = 0, beta = beta_0, missing_perc = current_dynamic_args$miss,
         missing_type = ifelse(grepl("nested", args$sim_name),
                               "monotone", "non-monotone"),
         estimator_type = args$est_type, extra_layer = args$extra_layer,
         learners = screen_plus_learner_lib,
         spvim_learners = spvim_lib, uni_learners = uni_learner_lib,
         b = args$b, M = args$m, alpha = per_comparison_alpha,
         thresh = sl_threshold, rank = sl_rank,
         ss_pfer = ss_pfer, ss_q = ss_q,
         kf_fdr = kf_fdr,
         pfp_q = pfp_q, gfwer_k = gfwer_k),
  simplify = FALSE)
)
selected_output_tib <- tibble::as_tibble(data.table::rbindlist(lapply(sim_output, function(x) x$selected)))
selected_output_tib %>%
  mutate(bias_init = sqrt(n) * (perf - true_perf)) %>%
  group_by(n, p, estimator_type) %>%
  summarize(mn_perf = mean(perf), truth = mean(true_perf), bias = mean(bias_init),
            .groups = "drop")
vim_output_tib <- tibble::as_tibble(data.table::rbindlist(lapply(sim_output, function(x) x$vim)))
unique_sim_name <- gsub("probit-", "", gsub("binomial-", "", args$sim_name))
output_dir <- paste0(prefix, gsub("+", "", args$est_type), "_",
                     args$extra_layer, "/",
                     unique_sim_name, "/")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
file_name <- paste0(args$sim_name, "_",
                    args$est_type, "_",
                    args$extra_layer,
                    "_b_", args$b, "_m_", args$m,
                    "_id_", job_id)
saveRDS(selected_output_tib,
        file = paste0(output_dir, file_name, "_select.rds"))
saveRDS(vim_output_tib,
        file = paste0(output_dir, file_name, "_vim.rds"))
print("Simulation complete!")
