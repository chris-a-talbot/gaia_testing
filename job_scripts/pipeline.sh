#!/bin/bash
#SBATCH --job-name=gaia_pipeline
#SBATCH --mail-user=chtalbot@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=1-10
#SBATCH --account=lsa1
#SBATCH --partition=standard
#SBATCH --output=./logs/hpc/%x-%A_%a.log

set -e  # Exit on error

# Setup environment
CWD=$(pwd)
REPLICATE_ID=${SLURM_ARRAY_TASK_ID}
DISPERSAL=0.3
RUNTIME=40000
REPLICATE_NAME="tree-S${DISPERSAL}-R${REPLICATE_ID}"

STARTED_AT=$(date)

# Load required modules
module use /nfs/turbo/lsa-bradburd/shared/Lmod/
module load Bioinformatics
module load SLiM/4.3
module load R/4.4.0
module load gcc/13.2.0
module load python3.11-anaconda

# Create logs directory if it doesn't exist
mkdir -p "${CWD}/logs/hpc"

# Initialize conda and activate environment
source $(conda info --base)/etc/profile.d/conda.sh
if ! conda activate gaia_testing; then
    echo "Failed to activate conda environment gaia_testing"
    exit 1
fi

echo "Starting pipeline for replicate ${REPLICATE_NAME} on node $(hostname) at $(date)"
echo "Current working directory: $CWD"
echo "SLURM Job ID: $SLURM_JOB_ID, Array Task ID: $SLURM_ARRAY_TASK_ID"

# Function to run a validation step that might return 1
run_validation_step() {
    local step_name=$1
    local command=$2

    echo "Starting ${step_name} at $(date)"
    if eval "$command"; then
        echo "Completed ${step_name} successfully (passed validation) at $(date)"
        return 0
    else
        echo "Completed ${step_name} (failed validation) at $(date)"
        return 1
    fi
}

# Function to run a regular command and check its exit status
run_step() {
    local step_name=$1
    local command=$2

    echo "Starting ${step_name} at $(date)"
    if eval "$command"; then
        echo "Completed ${step_name} successfully at $(date)"
        return 0
    else
        echo "Error in ${step_name} at $(date)"
        return 1
    fi
}

# Run SLiM simulation
if ! run_step "SLiM simulation" "slim -d \"PWD='$CWD'\" -d \"S=$DISPERSAL\" -d \"REP=$REPLICATE_ID\" -d \"RUNTIME=$RUNTIME\" ./scripts/run_slim.slim"; then
    echo "Pipeline failed at SLiM simulation"
    exit 1
fi

# Run coalescent check with special handling
if ! run_validation_step "Coalescent check" "python ./scripts/coalescent_check.py \"$REPLICATE_NAME\""; then
    echo "Coalescent check failed - simulation did not reach coalescence. Skipping remaining steps."
    exit 0  # Exit cleanly since this is an expected possibility
fi

# Run generate samples
if ! run_step "Generate samples" "python ./scripts/generate_samples.py \"$REPLICATE_NAME\""; then
    echo "Pipeline failed at generate samples"
    exit 1
fi

# Run subset tree
if ! run_step "Subset tree" "python ./scripts/subset_tree.py \"$REPLICATE_NAME\""; then
    echo "Pipeline failed at subset tree"
    exit 1
fi

# Run generate subset samples
if ! run_step "Generate subset samples" "python ./scripts/generate_samples.py \"$REPLICATE_NAME\" --subsets"; then
    echo "Pipeline failed at generate subset samples"
    exit 1
fi

# Run pre-GAIA validation with special handling
if ! run_validation_step "Validate pre-GAIA" "python ./scripts/validate_pre_gaia.py \"$REPLICATE_NAME\""; then
    echo "Pre-GAIA validation failed - data issues detected. Skipping remaining steps."
    exit 0  # Exit cleanly since this is an expected possibility
fi

# Run GAIA
if ! run_step "Run GAIA" "Rscript ./scripts/run_gaia.R \"$REPLICATE_NAME\""; then
    echo "Pipeline failed at GAIA"
    exit 1
fi

# Run GAIA analysis
if ! run_step "Analyze GAIA" "python ./scripts/analyze_gaia.py \"$REPLICATE_NAME\""; then
    echo "Pipeline failed at GAIA analysis"
    exit 1
fi

# Add at the end before final echo
COMPLETED_AT=$(date)
echo "
=== Pipeline Summary ===
Started at: ${STARTED_AT}
Completed at: ${COMPLETED_AT}
Tree name: ${REPLICATE_NAME}
Job ID: ${SLURM_JOB_ID}
Array task: ${SLURM_ARRAY_TASK_ID}
===================" >> "${CWD}/logs/hpc/${SLURM_JOB_NAME}-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"