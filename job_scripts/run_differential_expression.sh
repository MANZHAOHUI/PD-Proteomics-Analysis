#!/bin/bash
#SBATCH --job-name=pd_diffexp
#SBATCH --output=logs/diffexp_%j.out
#SBATCH --error=logs/diffexp_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=common

# Load required modules
module load R/4.1.0

# Run R analysis
Rscript scripts/R/PD_differential_expression_analysis.R \
    --input data/raw/proteomics_data.csv \
    --metadata data/raw/sample_info.csv \
    --output data/results/ \
    --cores 4

echo "Differential expression analysis completed!"
