#!/bin/bash

# submit all VIM/SL variable selection sims
# first, compute the number of jobs in array
num_n_p=7
njobs=`expr $2 / $3 \* $num_n_p`
echo $njobs
io_prefix=$7
io_file="$io_prefix/slurm-%A_%a.out"

echo -e \
    '#!/bin/bash\n Rscript run_sim_vim_ranks.R --sim-name $1 --nreps-total $2 --nreps-per-job $3 --b $4 --selection-type $5 --extra-layer $6' > call_sim_rank_select.sh
chmod u+x call_sim_rank_select.sh

sbatch --time=7-0 --array=1-$njobs -e $io_file -o $io_file ./call_sim_rank_select.sh $1 $2 $3 $4 $5 $6
