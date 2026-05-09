#!/bin/sh
export code_dir="/home/yanjin41/brmplus_simulation/compare"
export result_dir="/scratch/yanjin41/RRRDOR/brmplus_simulation"
#export N=$1

module load StdEnv/2023
module load gcc/12.3   r/4.3.1


n=$1
R=$2
R --vanilla --max-connections=512 CMD BATCH --no-save --no-restore "--args n=$n R=$R result_dir='${result_dir}'" $code_dir/run_rarediff.R $result_dir/Rout/RR_simucpp_N_${n}_R_${R}.Rout