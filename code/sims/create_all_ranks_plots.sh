#!/bin/bash

# create all plots from sims

# -----------------------------------------------------------
# create plots from linear + normal + nested
# -----------------------------------------------------------
./create_ranks_plots.sh "ranks-binomial-linear-normal" 20 1000

./create_ranks_plots.sh "ranks-binomial-nonlinear-normal" 20 1000

./create_ranks_plots.sh "ranks-binomial-weaklinear-normal" 20 1000
