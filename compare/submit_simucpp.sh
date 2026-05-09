#!/bin/sh
export code_dir="/home/yanjin41/brmplus_simulation/compare"
export result_dir="/scratch/yanjin41/RRRDOR/brmplus_simulation"
#export N=$1

module load StdEnv/2023
module load gcc/12.3   r/4.3.1


n=$1
R=$2
R --vanilla --max-connections=512 CMD BATCH --no-save --no-restore "--args n=$n R=$R result_dir='${result_dir}'" $code_dir/run_simulation.R $result_dir/Rout/RR_simucpp_N_${n}_R_${R}.Rout
#R CMD BATCH --no-save --no-restore "--args n=$n result_dir='${result_dir}'" $code_dir/exact_n50.R $result_dir/Rout/exact_result_N${n}.Rout


#for i in 50;do for j in 2;do sbatch --account=def-liteep -N 1 --ntasks-per-node=40 -t 30:00 -o $SCRATCH/RRRDOR/highd/highd_${i}_${j}.out -J exact_${i}_${j} $HOME/highd/submit_DY_sepjob.sh $i $j; done done
#for i in 50;do sbatch --account=def-liteep -N 1 --ntasks-per-node=40 -t 30:00 -o $SCRATCH/RRRDOR/exact/exact.out -J exact_${i} $HOME/brmplus/submit_DY_sepjob.sh $i ; done
#sbatch --account=def-liteep -N 1 --ntasks-per-node=40 -t 6:00:00 -o $SCRATCH/RRRDOR/exact/exact.out -J ATE_${i} $SCRATCH/project1/simulation/code/submit_DY_sepjob.sh $i
