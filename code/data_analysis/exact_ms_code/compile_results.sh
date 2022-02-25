#!/bin/bash

# compile all results
est_procs=("lasso_glm" "lasso-SS_glm" "lasso-KF_glm" \
    "SL_SL" "SL-SS_SL" "SPVIM_SL" "baseSL_SL" "SPVIM-RR_SL")
algos=("1" "2")
output_dir="../results/data_analysis/"

for est_proc in ${est_procs[@]}; do
    IFS='_' read -ra select_reg_array <<< "$est_proc"
    select=${select_reg_array[0]}
    regress=${select_reg_array[1]}
    select_dir=${select//;/-}
    for alg in ${algos[@]}; do
        if [ $select == "SPVIM-RR" ] && [ $alg == "2" ]; then
            # do nothing
            :
        else
            Rscript 06_compile_results.R --M 10 --est-type $regress --selection-type $select --algorithm $alg --output-dir $output_dir
        fi
    done
done
