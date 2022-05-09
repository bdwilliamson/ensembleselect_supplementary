#!/bin/bash

# create all plots from sims

# -----------------------------------------------------------
# create plots from linear + normal + nested (A)
# -----------------------------------------------------------
main_sim_names=("binomial-probit-linear-normal-nested" "nonlinear-normal-correlated")
for sim_name in ${main_sim_names[@]}; do
  ./create_plots.sh $sim_name 100 1 1000 0
  ./create_plots.sh $sim_name 100 1 1000 0
done
supp_sim_names=("linear-normal-correlated" "linear-normal-uncorrelated" \
                "nonlinear-normal-weak-uncorrelated" \
                "binomial-probit-linear-nonnormal-nested" \
                "binomial-probit-nonlinear-normal-nested" \
                "binomial-probit-nonlinear-nonnormal-nested")
for sim_name in ${supp_sim_names[@]}; do
  ./create_plots.sh $sim_name 100 1 1000 1
done

# ROC curves
sim_names=("binomial-probit-linear-normal-nested" "binomial-probit-nonlinear-nonnormal-nested" \
           "binomial-probit-linear-nonnormal-nested" "binomial-probit-nonlinear-normal-nested")
for sim_name in ${sim_names[@]}; do
  ./create_roc_curves.sh $sim_name 100 1 1000 1
done

# for old sims
./create_plots.sh "binomial-linear-normal-nested" 100 1 1000 &
./create_plots.sh "binomial-linear-nonnormal-nested" 100 1 1000 &
./create_plots.sh "binomial-nonlinear-normal-nested" 100 1 1000 &
./create_plots.sh "binomial-nonlinear-nonnormal-nested" 100 1 1000 &
