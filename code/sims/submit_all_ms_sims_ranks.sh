#!/bin/bash

# load modules
ml fhR/4.0.2-foss-2019b
ml jbigkit

# submit all sims for the ms
# --------------------------------------------------------------------------------------------------------
# SETTING A: Outcome regression is a linear model, X ~ Normal
# --------------------------------------------------------------------------------------------------------
#   first, SL with all variables
./submit_sim_rank_select.sh "ranks-binomial-linear-normal" 1000 10 20 "sl_none" "none" "sl_none_none"

./submit_sim_rank_select.sh "ranks-binomial-linear-normal" 1000 10 20 "sl_none" "screen" "sl_none_screen"

#   second, SL with screens
./submit_sim_rank_select.sh "ranks-binomial-linear-normal" 1000 10 20 "sl_screen" "none" "sl_screen_none"

./submit_sim_rank_select.sh "ranks-binomial-linear-normal" 1000 10 20 "sl_screen" "KF" "sl_screen_kf"

./submit_sim_rank_select.sh "ranks-binomial-linear-normal" 1000 10 20 "sl_screen" "SS" "sl_screen_ss"

#   third, SL with rank-based selection
./submit_sim_rank_select.sh "ranks-binomial-linear-normal" 1000 10 20 "sl_rank" "none" "sl_rank_none"

./submit_sim_rank_select.sh "ranks-binomial-linear-normal" 1000 10 20 "sl_rank" "KF" "sl_rank_kf"

./submit_sim_rank_select.sh "ranks-binomial-linear-normal" 1000 10 20 "sl_rank" "SS" "sl_rank_ss"

#   fourth, SPVIM ranks
./submit_sim_rank_select.sh "ranks-binomial-linear-normal" 1000 10 20 "spvim_rank" "none" "spvim_rank"

./submit_sim_rank_select.sh "ranks-binomial-linear-normal" 1000 10 20 "spvim_BH" "none" "spvim_BH"

#   last, SAGE ranks
./submit_sim_rank_select.sh "ranks-binomial-linear-normal" 1000 10 20 "sage_rank" "none" "sage_rank"

# --------------------------------------------------------------------------------------------------------
# SETTING B: Outcome regression is a nonlinear model, X ~ Normal
# --------------------------------------------------------------------------------------------------------
#   first, SL with all variables
./submit_sim_rank_select.sh "ranks-binomial-nonlinear-normal" 1000 10 20 "sl_none" "none" "sl_none_none"

./submit_sim_rank_select.sh "ranks-binomial-nonlinear-normal" 1000 10 20 "sl_none" "screen" "sl_none_screen"

#   second, SL with screens
./submit_sim_rank_select.sh "ranks-binomial-nonlinear-normal" 1000 10 20 "sl_screen" "none" "sl_screen_none"

./submit_sim_rank_select.sh "ranks-binomial-nonlinear-normal" 1000 10 20 "sl_screen" "KF" "sl_screen_kf"

./submit_sim_rank_select.sh "ranks-binomial-nonlinear-normal" 1000 10 20 "sl_screen" "SS" "sl_screen_ss"

#   third, SL with rank-based selection
./submit_sim_rank_select.sh "ranks-binomial-nonlinear-normal" 1000 10 20 "sl_rank" "none" "sl_rank_none"

./submit_sim_rank_select.sh "ranks-binomial-nonlinear-normal" 1000 10 20 "sl_rank" "KF" "sl_rank_kf"

./submit_sim_rank_select.sh "ranks-binomial-nonlinear-normal" 1000 10 20 "sl_rank" "SS" "sl_rank_ss"

#   fourth, SPVIM ranks
./submit_sim_rank_select.sh "ranks-binomial-nonlinear-normal" 1000 10 20 "spvim_rank" "none" "spvim_rank"

./submit_sim_rank_select.sh "ranks-binomial-nonlinear-normal" 1000 10 20 "spvim_BH" "none" "spvim_BH"

#   last, SAGE ranks
./submit_sim_rank_select.sh "ranks-binomial-nonlinear-normal" 1000 10 20 "sage_rank" "none" "sage_rank"

# --------------------------------------------------------------------------------------------------------
# SETTING E: Outcome regression is a weak linear model, X ~ Normal
# --------------------------------------------------------------------------------------------------------
#   first, SL with all variables
./submit_sim_rank_select.sh "ranks-binomial-weaklinear-normal" 1000 10 20 "sl_none" "none" "sl_none_none"

./submit_sim_rank_select.sh "ranks-binomial-weaklinear-normal" 1000 10 20 "sl_none" "screen" "sl_none_screen"

#   second, SL with screens
./submit_sim_rank_select.sh "ranks-binomial-weaklinear-normal" 1000 10 20 "sl_screen" "none" "sl_screen_none"

./submit_sim_rank_select.sh "ranks-binomial-weaklinear-normal" 1000 10 20 "sl_screen" "KF" "sl_screen_kf"

./submit_sim_rank_select.sh "ranks-binomial-weaklinear-normal" 1000 10 20 "sl_screen" "SS" "sl_screen_ss"

#   third, SL with rank-based selection
./submit_sim_rank_select.sh "ranks-binomial-weaklinear-normal" 1000 10 20 "sl_rank" "none" "sl_rank_none"

./submit_sim_rank_select.sh "ranks-binomial-weaklinear-normal" 1000 10 20 "sl_rank" "KF" "sl_rank_kf"

./submit_sim_rank_select.sh "ranks-binomial-weaklinear-normal" 1000 10 20 "sl_rank" "SS" "sl_rank_ss"

#   fourth, SPVIM ranks
./submit_sim_rank_select.sh "ranks-binomial-weaklinear-normal" 1000 10 20 "spvim_rank" "none" "spvim_rank"

./submit_sim_rank_select.sh "ranks-binomial-weaklinear-normal" 1000 10 20 "spvim_BH" "none" "spvim_BH"

#   last, SAGE ranks
./submit_sim_rank_select.sh "ranks-binomial-weaklinear-normal" 1000 10 20 "sage_rank" "none" "sage_rank"
