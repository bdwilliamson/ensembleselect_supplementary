#!/bin/bash

# -----------------------------------------------------------
# create plots from sim results
# -----------------------------------------------------------
# accepts the following command-line arguments
# 1: sim name (e.g., "binomial-linear-normal-nested")
# 2: number of bootstrap replicates
# 3: number of imputes
# 4: number of total jobs
# 5: should we put all missing data in? 0/1
sim_name=$1
b=$2
m=$3
nreps_total=$4
output_file="compile_output/${1}_roc.out"

Rscript create_roc_curves.R --sim-name "$sim_name" --b "$b" --m "$m" --nreps-total "$nreps_total" --all-miss $5 > $output_file 2>&1
