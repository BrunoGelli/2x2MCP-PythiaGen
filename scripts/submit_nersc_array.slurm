#!/bin/bash
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -J mcp_emitters
#SBATCH -o logs/slurm/%x_%A_%a.out
#SBATCH -e logs/slurm/%x_%A_%a.err
# Set account/site options before submitting, e.g. #SBATCH -A <account>
set -euo pipefail
mkdir -p logs/slurm
export JOB_LIST=${JOB_LIST:-job_list.csv}
export EXEC=${EXEC:-./allInOne_multimeson}
export NTHREADS=${NTHREADS:-1}
./scripts/run_one_job.sh "${SLURM_ARRAY_TASK_ID}"
