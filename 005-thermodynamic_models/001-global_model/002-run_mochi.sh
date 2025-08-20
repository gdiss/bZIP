#!/bin/bash
#SBATCH --account=thermo_bz
#SBATCH --job-name=mochi
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=main
#SBATCH --mem=200G
#SBATCH --gres=gpu:v100:1
#SBATCH --time=2-6

set -eu
 
function display_memory_usage() {
        set +eu
        echo -n "[INFO] [$(date -Iseconds)] [$$] Max memory usage in bytes: "
        cat /sys/fs/cgroup/memory/slurm/uid_$(id -u)/job_${SLURM_JOB_ID}/memory.max_usage_in_bytes
        echo
}
 
trap display_memory_usage EXIT
 
START=$(date +%s)
STARTDATE=$(date -Iseconds)
echo "[INFO] [$STARTDATE] [$$] Starting SLURM job $SLURM_JOB_ID"
echo "[INFO] [$STARTDATE] [$$] Running in $(hostname -s)"
echo "[INFO] [$STARTDATE] [$$] Working directory: $(pwd)"
 
### PUT YOUR CODE IN THIS SECTION
module purge
source /tachyon/scratch/gdiss/software/miniconda3/bin/activate
conda activate pymochi

export NUMEXPR_NUM_THREADS=10
export MKL_NUM_THREADS=10
export MP_NUM_THREADS=10
export OPENBLAS_NUM_THREADS=10
export VECLIB_MAXIMUM_THREADS=10
export GOTO_NUM_THREADS=10
export BLIS_NUM_THREADS=10

output_directory="002-output/"
design_file="002-design.txt"
stdout="002-stdout.txt"
stderr="002-stderr.txt"

run_mochi.py --output_directory "$output_directory" --model_design "$design_file" 1>$stdout 2>$stderr

 
### END OF PUT YOUR CODE IN THIS SECTION
 
END=$(date +%s)
ENDDATE=$(date -Iseconds)
echo "[INFO] [$ENDDATE] [$$] Workflow execution time \(seconds\) : $(( $END-$START ))"