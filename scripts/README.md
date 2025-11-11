# Analysis Scripts

This directory contains the main analysis scripts for the PD Proteomics project.

## Python Scripts (`python/`)

### PD_network_analysis.py
- MEGENA-like co-expression network construction
- Module detection using multiscale clustering
- Key driver identification
- Bayesian causal network inference
- eQTL analysis integration

**Usage:**
```bash
python scripts/python/PD_network_analysis.py \
    --input data/processed/expression_matrix.csv \
    --output data/results/
```

## R Scripts (`R/`)

### PD_differential_expression_analysis.R
- Complete pipeline for differential expression analysis
- Batch correction using ComBat
- limma/edgeR implementation
- Pathway enrichment analysis
- Sex-stratified analysis

**Usage:**
```bash
Rscript scripts/R/PD_differential_expression_analysis.R \
    --input data/raw/ \
    --output data/results/
```

## Running on HPC

For HPC execution, submit jobs using SLURM:

```bash
sbatch job_scripts/run_network_analysis.sh
```

See the `job_scripts/` directory for example SLURM submission scripts.
