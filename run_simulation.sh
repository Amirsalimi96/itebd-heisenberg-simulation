#!/bin/bash
#SBATCH --partition=regular
#SBATCH --job-name="iTEBD-Heisenberg"
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=04:00:00
#SBATCH --output=simulation_%j.log
#SBATCH --error=simulation_%j.err

# ==============================================================================
# iTEBD Ground State Simulation Batch Runner (SLURM / HPC Cluster)
# ==============================================================================

echo "Job started at: $(date)"

# Load Anaconda environment (uncomment/modify according to your cluster setup)
# module load Anaconda/Anaconda3
# source activate itebd-env

# Target simulation Python script
SIMULATION_FILE="itebd.py"

# Define output directory
OUTPUT_DIR="results"
mkdir -p "$OUTPUT_DIR"

echo "Running iTEBD simulation..."
echo "Python script: $SIMULATION_FILE"

# Execute simulation script
python3 "$SIMULATION_FILE" > "$OUTPUT_DIR/simulation_output.log" 2>&1

echo "Simulation completed successfully."
echo "Job finished at: $(date)"
