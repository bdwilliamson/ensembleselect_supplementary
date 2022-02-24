#!/bin/bash

# -----------------------------------------------------------
# create plots from sim results
# -----------------------------------------------------------
# accepts the following command-line arguments
# 1: sim name (e.g., "binomial-linear-normal-nested")
# 2: number of bootstrap replicates
# 3: number of imputes
# 4: number of total jobs
sim_name=$1
b=$2
nreps_total=$3

Rscript create_ranks_plots.R --sim-name "$sim_name" --b "$b" --nreps-total "$nreps_total"
