#!/usr/local/bin/Rscript
# compare variable selection based on various VIM procedures:
#   (a) selection based on ranks
#   (b) selection based on BH FDR-controlled p-value approach
#   (c) selection based on knockoffs?
# also, test coverage of VIM rank intervals

# load required functions and packages
library("boot")
library("tidyr")
library("dplyr")
library("stringr")
library("purrr")
library("tibble")
library("SuperLearner")
library("glmnet")
library("xgboost")
library("ranger")
library("argparse")
library("knockoff")
library("parallel")
library("methods")
library("stabs")
library("nloptr")

# pull in job id and set up command-line args
if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
  # if you need to update on cluster: remotes::install_github("bdwilliamson/vimp", upgrade = "never", lib = .libPaths()[2])
  library("vimp")
  job_id <- 1
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
source(paste0(code_dir, "rank_once.R"))
source(paste0(code_dir, "gen_data.R"))
source(paste0(code_dir, "utils.R"))
source(paste0(code_dir, "sl_knockoff_stats.R"))

parser <- ArgumentParser()
parser$add_argument("--sim-name", default = "ranks-binomial-linear-normal",
                    help = "the name of the simulation")
parser$add_argument("--nreps-total", type = "double", default = 1000,
                    help = "number of replicates in total")
parser$add_argument("--nreps-per-job", type = "double", default = 10,
                    help = "number of replicates for each job")
parser$add_argument("--b", type = "double", default = 20,
                    help = "number of bootstrap replicates")
parser$add_argument("--selection-type", default = "sl_rank",
                    help = "type of variable selection to use ('sl_rank', 'spvim_rank', 'sage_rank', or 'spvim_BH')")
parser$add_argument("--extra-layer", default = "none",
                    help = "extra layer ('none' or 'KF' or 'SS')")
args <- parser$parse_args()
print(paste0("Running sim ", args$sim_name, " with ", as.character(args$nreps_total), " total replicates and ",
      as.character(args$nreps_per_job), " replicates per job."))
print(paste0("Using ", args$selection_type, " to do feature selection with extra layer ", args$extra_layer, "."))
print(paste0("Additional parameters: number of bootstrap replicates = ", args$b, "."))

# set up all of the possible simulation parameters
ps <- 100
ns <- c(200, 350, 500, 1000, 2000, 3000, 4000)
families <- c("binomial")
nreps_per_combo <- args$nreps_total/args$nreps_per_job
param_grid <- expand.grid(mc_id = 1:nreps_per_combo, p = ps, n = ns)

# grab the current simulation parameters
current_dynamic_args <- param_grid[job_id, ]

# set up the beta vector
if (grepl("weak", args$sim_name)) {
  beta_active <- c(2.5, -1, -0.5, 0.5, -0.25, -0.25, 0.25)
} else {
  beta_active <- c(2.5, -3, -1, 1, -1.5, -0.5, 0.5)
}
beta_nonactive <- 0
active_set <- 1:7
non_active_set <- (1:current_dynamic_args$p)[-active_set]
beta_0 <- vector("numeric", current_dynamic_args$p + 1)
beta_0[active_set] <- beta_active
beta_0[non_active_set] <- beta_nonactive

xgb_tune_params <- list(max_depth = 1, shrinkage = 0.1)
xgb_learners <- create.Learner("SL.xgboost_new", tune = xgb_tune_params, detailed_names = TRUE, name_prefix = "xgb")
# no screens -- testing what SL ranking does on its own
learner_lib <- c(xgb_learners$names, "SL.glm", "SL.ranger.imp", "SL.glmnet", "SL.mean")
# fit polymars for univariate SPVIM regressions
uni_learner_lib <- "my_SL.polymars"
if (grepl("spvim", args$selection_type)) {
  # knock out xgboost, to *hopefully* decrease computational overhead
  learner_lib <- learner_lib[!grepl("xgboost", learner_lib)]
}
if (grepl("screen", args$selection_type) | grepl("screen", args$extra_layer)) {
  # add screens to SL lib
  screens <- c("screen.glmnet", "screen.corRank.25")
  learner_plus_screen_lib <- lapply(as.list(learner_lib), function(x) c(x, screens))
  learner_lib <- learner_plus_screen_lib
}
## ---------------------------------------------
## replicate the simulation nreps_per_job times
## ---------------------------------------------
current_seed <- job_id
print(current_seed)
set.seed(current_seed)
system.time(sim_output <- replicate(args$nreps_per_job, do_one(n = current_dynamic_args$n, p = current_dynamic_args$p,
                family = ifelse(grepl("binomial", args$sim_name), "binomial", "gaussian"),
                linear = !grepl("nonlinear", args$sim_name),
                x_dist = ifelse(!grepl("nonnormal", args$sim_name), "normal", "nonnormal"),
                rho = 0,
                beta = beta_0,
                estimator_type = "SL",
                extra_layer = args$extra_layer,
                selection_type = args$selection_type,
                learners = learner_lib,
                uni_learners = uni_learner_lib,
                b = args$b,
                thresh = 0.2,
                rank = 25,
                q = 25,
                fdr = .75,
                offset = 0,
                knockoff_method = 'asdp'), simplify = FALSE)
)
sim_output_2 <- lapply(as.list(1:length(sim_output)), function(x) tibble::add_column(sim_output[[x]], mc_id = x + args$nreps_per_job * (job_id - 1)))
sim_output_tib <- do.call(rbind.data.frame, sim_output_2)
file_name <- paste0(args$sim_name, "_",
                    args$selection_type, "_",
                    args$extra_layer,
                    "_b_", args$b,
                    "_id_", job_id, ".rds")
output_dir <- paste0(prefix, args$selection_type, "_", args$extra_layer, "/")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
saveRDS(sim_output_tib, file = paste0(output_dir, file_name))
