# R Package Dependencies for PD Proteomics Analysis
# Run this script to install all required R packages

# CRAN packages
cran_packages <- c(
  # Data manipulation
  "tidyverse",
  "dplyr",
  "tidyr",
  "readr",
  "data.table",
  
  # Statistical analysis
  "limma",
  "edgeR",
  "DESeq2",
  "sva",
  
  # Visualization
  "ggplot2",
  "pheatmap",
  "ComplexHeatmap",
  "EnhancedVolcano",
  "ggrepel",
  "cowplot",
  "patchwork",
  "RColorBrewer",
  "viridis",
  
  # Network analysis
  "igraph",
  "WGCNA",
  "bnlearn",
  
  # Pathway analysis
  "clusterProfiler",
  "enrichplot",
  "DOSE",
  "pathview",
  
  # Utilities
  "openxlsx",
  "optparse",
  "foreach",
  "doParallel"
)

# Install CRAN packages
install.packages(cran_packages, dependencies = TRUE)

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

bioc_packages <- c(
  "limma",
  "edgeR", 
  "DESeq2",
  "sva",
  "impute",
  "preprocessCore",
  "org.Hs.eg.db",
  "GO.db",
  "KEGG.db",
  "clusterProfiler",
  "enrichplot",
  "DOSE",
  "pathview",
  "ComplexHeatmap",
  "EnhancedVolcano",
  "Rgraphviz"
)

BiocManager::install(bioc_packages, update = FALSE)

# Verify installation
cat("\n=== Verifying package installation ===\n")
all_packages <- c(cran_packages, bioc_packages)
installed <- sapply(all_packages, requireNamespace, quietly = TRUE)
if (all(installed)) {
  cat("All packages successfully installed!\n")
} else {
  cat("The following packages failed to install:\n")
  print(names(installed)[!installed])
}
