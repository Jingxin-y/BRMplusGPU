#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=1
#SBATCH --time=6:00:00
#SBATCH --job-name=brm_gpu

set -euo pipefail

export package_dir="${package_dir:-/home/yanjin41/brmplus_simulation}"
export code_dir="${code_dir:-${package_dir}/compare}"
export result_dir="${result_dir:-/scratch/yanjin41/RRRDOR/brmplus_simulation}"

module load StdEnv/2023
module load gcc/12.3 r/4.3.1 cuda

export BRM_USE_CUDA=1
export BRM_USE_GPU=1
export BRM_GPU_WORKERS="${BRM_GPU_WORKERS:-1}"
export BRM_EXACT_PARALLEL="${BRM_EXACT_PARALLEL:-1}"
export BRM_EXACT_WORKERS="${BRM_EXACT_WORKERS:-${SLURM_CPUS_PER_TASK:-4}}"
export CUDA_LIB="${CUDA_HOME}/lib64"
export CUDA_ARCH="${CUDA_ARCH:-sm_90}"

mkdir -p "${result_dir}/Rout"

if [ "${BRM_INSTALL_ON_JOB:-1}" = "1" ]; then
  R CMD INSTALL "${package_dir}"
fi

n=$1
R=$2

R --vanilla --max-connections=512 CMD BATCH --no-save --no-restore \
  "--args n=${n} R=${R} result_dir='${result_dir}'" \
  "${code_dir}/run_simulation.R" \
  "${result_dir}/Rout/RR_simucpp_GPU_N_${n}_R_${R}.Rout"
