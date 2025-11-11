#!/bin/bash
#SBATCH --job-name=pd_network
#SBATCH --output=logs/network_%j.out
#SBATCH --error=logs/network_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=common

# Load required modules
module load Python/3.9.6
module load GCC/11.2.0

# Activate virtual environment (if using)
# source /path/to/venv/bin/activate

# Run analysis
python scripts/python/PD_network_analysis.py \
    --input data/processed/expression_matrix.csv \
    --metadata data/processed/sample_metadata.csv \
    --output data/results/ \
    --threads 8

echo "Network analysis completed!"
