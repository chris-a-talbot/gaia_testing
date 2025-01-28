#!/bin/bash
#SBATCH --job-name=gaia_testing_simulation
#SBATCH --mail-user=chtalbot@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=24G
#SBATCH --time=6:00:00
#SBATCH --array=1-10
#SBATCH --account=lsa1
#SBATCH --partition=standard
#SBATCH --output=./logs/%x-%A_%a.log

set -e

CWD=$(pwd)

module use /nfs/turbo/lsa-bradburd/shared/Lmod/
module load Bioinformatics
module load SLiM/4.3

echo "Starting simulation REP=${SLURM_ARRAY_TASK_ID} on node $(hostname) at $(date)"
echo "Current working directory: $CWD"
echo "SLURM Job ID: $SLURM_JOB_ID, Array Task ID: $SLURM_ARRAY_TASK_ID"

slim -d "PWD='$CWD'" -d "S=0.3" -d "REP=${SLURM_ARRAY_TASK_ID}" -d "RUNTIME=45000" simulation.slim > "${CWD}/logs/simulation_rep_${SLURM_ARRAY_TASK_ID}.log" 2>&1

echo "Completed simulation REP=${SLURM_ARRAY_TASK_ID} at $(date)"
