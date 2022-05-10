# Running the numerical experiments for the `flevr` paper

This repository contains code to reproduce the analyses in ["Flexible variable selection in the presence of missing data"](https://arxiv.org/abs/2202.12989) by Williamson and Huang (2022+). All analyses were implemented in the freely available R programming language; specifically, version 4.0.2. All analyses use the R package `flevr` version 0.0.2.

The numerical experiments consist of two sections. First, we consider the properties of our proposal under four data-generating mechanisms. Second, we investigate the properties of extrinsic selection incorporating both knockoffs and variable screens.

## Properties under various data-generating mechanisms

The following code will replicate the results in Section 4 of the main manuscript (Figures 1 and 2) and Section 2 of the Supplementary Material (Figures S1--S35).

The simulation uses the following files:
* `submit_all_ms_sims.sh`: Submit all simulations for this section.
* `submit_sim_feature_select.sh`: Batch submission of a group of jobs to a Slurm scheduler.
* `run_sim_feature_select.R`: the main R script for this simulation, corresponding to Scenarios 1 and 3--5. Runs the simulation `nreps_per_job` times for a specified set of parameters.
* `investigate_lasso_performance.R`: the second main R script for this simulation, corresponding to Scenarios 2 and 6--8. Runs the simulation `nreps_per_job` times for a specified set of parameters.
* `select_once.R`: Runs the simulation a single time for a specified set of parameters.
* `gen_data.R`: Generate a dataset.
* `fit_extra_layers.R`: Do variable selection based on the specified algorithm (e.g., lasso, lasso + knockoffs, SL).
* `get_true_performance.R`: Get prediction performance of a selected set of variables on independent data.
* `utils.R`: Utility functions.

Running the following code will submit all of the simulations to a Slurm scheduler:
```{bash}
chmod u+x *.sh
./submit_all_ms_sims.sh
```
If you aren't running on a Slurm scheduler, make sure to edit the appropriate lines (flagged in each file). You can run this code locally, but it will take some time.

Once you have run the simulations and copied the results to directory `sim_output`, the following code reproduces all plots and tables:
* `create_all_plots.sh`: creates all plots for the main manuscript and supplement
* `create_plots.sh`: run `create_plots.R` for a given scenario
* `create_plots.R`: load in the results and create figures for variable selection and prediction performance
* `create_roc_curves.sh`: run `create_roc_curves.R` for a given scenario
* `create_roc_curves.R`: create figures with ROC curves for each estimation procedure and scenario

```{bash}
./create_all_plots.sh
```

## Properties of our extrinsic selector with knockoffs and variable screens

The following code will replicate the results in Section 2.8 of the Supplementary Material (Figures S36--S41).

The simulation uses the following files:
* `submit_all_ms_sims_ranks.sh`: Submit all simulations for this section.
* `submit_sim_rank_select.sh`: Batch submission of a group of jobs to a Slurm scheduler.
* `run_sim_rank_select.R`: the main R script for this simulation. Runs the simulation `nreps_per_job` times for a specified set of parameters.
* `rank_once.R`: Runs the simulation a single time for a specified set of parameters.
* `gen_data.R`: Generate a dataset.
* `sl_knockoff_stats.R`: Knockoff statistic for Super Learner-based knockoffs.
* `utils.R`: Utility functions.

Running the following code will submit all of the simulations to a Slurm scheduler:
```{bash}
chmod u+x *.sh
./submit_all_ms_sims_ranks.sh
```
If you aren't running on a Slurm scheduler, make sure to edit the appropriate lines (flagged in each file). You can run this code locally, but it will take some time.

Once you have run the simulations and copied the results to directory `sim_output`, the following code reproduces all plots and tables:
* `create_all_ranks_plots.sh`: creates all plots for the supplement
* `create_ranks_plots.sh`: run `create_ranks_plots.R` for a given scenario
* `create_ranks_plots.R`: load in the results and create figures

```{bash}
./create_all_ranks_plots.sh
```
