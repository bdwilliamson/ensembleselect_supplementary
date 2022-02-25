#!/bin/bash

# submit all sims for the ms
# first, compute the number of jobs in the job array (if I haven't passed one in)
if [ "${10}" == "" ]; then
  num_n_p_miss=$9
  njobs=`expr $2 / $3 \* $num_n_p_miss`
  arry="1-$njobs"
else
  arry=${10}
fi

# edit the next line to point to your preferred i/o directory
io_prefix="<path to preferred i/o directory>/$8"
fi
mkdir -p $io_prefix
io_file="$io_prefix/slurm-%A_%a.out"

# Takes 7 command-line arguments
# 1: simulation name (e.g., "binomial-linear-normal-nested")
# 2: number of total reps
# 3: number of reps per job
# 4: number of boostrap reps (for stability selection)
# 5: number of MI reps
# 6: estimator ('lasso' or 'SL')
# 7: extra layer ('' or 'SS' or 'knockoffs')
echo -e \
  '#!/bin/bash\n Rscript run_sim_feature_select.R --sim-name $1' \
  '--nreps-total $2 --nreps-per-job $3 --b $4 --m $5 --est-type $6' \
  '--extra-layer $7' > call_sim_feature_select.sh
chmod u+x call_sim_feature_select.sh
# run the sim
sbatch --time=7-0 --array=$arry -e $io_file -o $io_file \
./call_sim_feature_select.sh $1 $2 $3 $4 $5 $6 $7 ${12}
# clean up
rm call_sim_feature_select.sh
