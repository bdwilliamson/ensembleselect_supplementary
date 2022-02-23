#!/bin/bash

# submit all sims for the ms; submit on an estimator-specific basis
est_layers=("lasso_none" "base-SL_none" "SL_none" "lasso_KF" "lasso_SS" "SL_SS" "SPVIM_none")
# settings: A, D, B, C
sim_prefix="binomial-probit-"
sim_suffix="-nested"
sim_names=("linear-normal" "nonlinear-nonnormal" "linear-nonnormal" "nonlinear-normal")
for est_layer in "${est_layers[@]}"; do
  IFS='_' read -ra est_layer_array <<< "$est_layer"
  est=${est_layer_array[0]}
  layer=${est_layer_array[1]}
  # if not SS or SPVIM, use 5 reps per job
  if [ $layer == "SS" ] || [ $est == "SPVIM" ]; then
    nreps_per_job=2
  else
    nreps_per_job=5
  fi
  # submit
  output_est=$(echo "$est_layer" | tr '[:upper:]' '[:lower:]')
  for sim_name in "${sim_names[@]}"; do
    echo "Sim: $sim_name; est: $est_layer"
    full_sim_name="${sim_prefix}${sim_name}${sim_suffix}"
    io_prefix="mi_predictiveness/output_${output_est}_${sim_name}${sim_suffix}"
    ./submit_sim_feature_select.sh $full_sim_name 1000 $nreps_per_job 100 1 $est $layer $io_prefix 24 ""
  done
done
