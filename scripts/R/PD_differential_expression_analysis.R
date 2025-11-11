#!/usr/bin/env Rscript
# ==============================================================================
# Differential Expression Analysis for Parkinson's Disease Proteomics
# ==============================================================================
# This script performs differential expression analysis on PD proteomic datasets
# Following the pipeline described in the main manuscript
# Author: Research Team
# Date: 2025
# ==============================================================================

# Load required libraries
library(limma)
library(edgeR)
library(DESeq2)
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(EnhancedVolcano)

# Set working directory
setwd("/path/to/pd_proteomics/")

# ==============================================================================
# 1. DATA LOADING AND PREPROCESSING
# ==============================================================================

# Load proteomic data from TMT experiments
load_proteomic_data <- function(cohort_name) {
  # Example for loading PPMI cohort data
  if (cohort_name == "PPMI") {
    # Load from publicly available PPMI repository
    data <- read.csv("data/PPMI_proteomics_tmt.csv", row.names = 1)
    metadata <- read.csv("data/PPMI_metadata.csv")
  } else if (cohort_name == "AMP-PD") {
    # Load from AMP-PD portal
    data <- read.csv("data/AMPPD_proteomics_tmt.csv", row.names = 1)
    metadata <- read.csv("data/AMPPD_metadata.csv")
  } else if (cohort_name == "GNPC") {
    # Load from GNPC database
    data <- read.csv("data/GNPC_proteomics_tmt.csv", row.names = 1)
    metadata <- read.csv("data/GNPC_metadata.csv")
  }
  
  return(list(data = data, metadata = metadata))
}

# Quality control function
perform_qc <- function(data, metadata) {
  # Remove samples with > 50% missing values
  missing_per_sample <- colSums(is.na(data)) / nrow(data)
  good_samples <- names(missing_per_sample[missing_per_sample < 0.5])
  
  # Remove proteins detected in < 50% of samples
  detection_rate <- rowSums(!is.na(data)) / ncol(data)
  good_proteins <- rownames(data)[detection_rate > 0.5]
  
  # Filter data
  data_filtered <- data[good_proteins, good_samples]
  metadata_filtered <- metadata[metadata$Sample_ID %in% good_samples, ]
  
  # Impute missing values using KNN
  library(impute)
  data_imputed <- impute.knn(as.matrix(data_filtered), k = 10)$data
  
  return(list(data = data_imputed, metadata = metadata_filtered))
}

# Batch correction using ComBat
correct_batch_effects <- function(data, metadata) {
  library(sva)
  
  # Create model matrix
  mod <- model.matrix(~ Disease_Status + Age + Sex, data = metadata)
  
  # Apply ComBat
  data_corrected <- ComBat(dat = data, 
                           batch = metadata$Batch, 
                           mod = mod,
                           par.prior = TRUE)
  
  return(data_corrected)
}

# ==============================================================================
# 2. DIFFERENTIAL EXPRESSION ANALYSIS
# ==============================================================================

perform_differential_expression <- function(data, metadata, comparison) {
  # Create design matrix
  if (comparison == "PD_vs_Control") {
    design <- model.matrix(~ 0 + Disease_Status + Age + Sex + PMI, data = metadata)
    colnames(design) <- make.names(colnames(design))
  } else if (comparison == "Braak_stages") {
    design <- model.matrix(~ 0 + Braak_Stage + Age + Sex + PMI, data = metadata)
  }
  
  # Fit linear model using limma
  fit <- lmFit(data, design)
  
  # Define contrasts
  if (comparison == "PD_vs_Control") {
    contrast.matrix <- makeContrasts(
      PD_vs_Control = Disease_StatusPD - Disease_StatusControl,
      levels = design
    )
  } else if (comparison == "Braak_stages") {
    contrast.matrix <- makeContrasts(
      Stage6_vs_0 = Braak_Stage6 - Braak_Stage0,
      Stage5_vs_0 = Braak_Stage5 - Braak_Stage0,
      Stage4_vs_0 = Braak_Stage4 - Braak_Stage0,
      levels = design
    )
  }
  
  # Apply contrasts
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # Extract results
  results <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")
  
  # Add gene annotations
  library(AnnotationDbi)
  results$Gene_Symbol <- mapIds(org.Hs.eg.db,
                                keys = rownames(results),
                                column = "SYMBOL",
                                keytype = "UNIPROT",
                                multiVals = "first")
  
  return(results)
}

# ==============================================================================
# 3. PATHWAY ENRICHMENT ANALYSIS
# ==============================================================================

perform_pathway_enrichment <- function(dep_results, fc_cutoff = 0.25, fdr_cutoff = 0.05) {
  # Get significant DEPs
  sig_up <- dep_results %>%
    filter(logFC > fc_cutoff, adj.P.Val < fdr_cutoff) %>%
    pull(Gene_Symbol)
  
  sig_down <- dep_results %>%
    filter(logFC < -fc_cutoff, adj.P.Val < fdr_cutoff) %>%
    pull(Gene_Symbol)
  
  # GO enrichment for upregulated proteins
  ego_up <- enrichGO(gene = sig_up,
                     OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05)
  
  # GO enrichment for downregulated proteins
  ego_down <- enrichGO(gene = sig_down,
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05)
  
  # KEGG pathway enrichment
  library(DOSE)
  kegg_up <- enrichKEGG(gene = sig_up,
                        organism = "hsa",
                        pvalueCutoff = 0.05)
  
  kegg_down <- enrichKEGG(gene = sig_down,
                         organism = "hsa",
                         pvalueCutoff = 0.05)
  
  return(list(
    go_up = ego_up,
    go_down = ego_down,
    kegg_up = kegg_up,
    kegg_down = kegg_down
  ))
}

# ==============================================================================
# 4. VISUALIZATION FUNCTIONS
# ==============================================================================

# Volcano plot
create_volcano_plot <- function(dep_results, title = "PD vs Control") {
  library(ggplot2)
  library(ggrepel)
  
  # Add significance categories
  dep_results$Significance <- "Not Significant"
  dep_results$Significance[dep_results$adj.P.Val < 0.05 & abs(dep_results$logFC) > 0.25] <- "Significant"
  dep_results$Significance[dep_results$adj.P.Val < 0.01 & abs(dep_results$logFC) > 0.5] <- "Highly Significant"
  
  # Label top proteins
  top_proteins <- dep_results %>%
    filter(adj.P.Val < 0.001) %>%
    arrange(desc(abs(logFC))) %>%
    head(20)
  
  p <- ggplot(dep_results, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = Significance), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("gray50", "blue", "red")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray30") +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "gray30") +
    geom_text_repel(data = top_proteins,
                    aes(label = Gene_Symbol),
                    size = 3,
                    max.overlaps = 20) +
    theme_classic() +
    labs(title = title,
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  return(p)
}

# Heatmap of top DEPs
create_heatmap <- function(data, dep_results, metadata, top_n = 50) {
  # Get top differentially expressed proteins
  top_deps <- dep_results %>%
    arrange(adj.P.Val) %>%
    head(top_n) %>%
    rownames()
  
  # Subset data
  heatmap_data <- data[top_deps, ]
  
  # Scale rows
  heatmap_data_scaled <- t(scale(t(heatmap_data)))
  
  # Create annotation dataframe
  annotation_col <- data.frame(
    Disease_Status = metadata$Disease_Status,
    Brain_Region = metadata$Brain_Region,
    Sex = metadata$Sex,
    row.names = colnames(heatmap_data)
  )
  
  # Define colors
  ann_colors <- list(
    Disease_Status = c(Control = "blue", PD = "red"),
    Brain_Region = c(SN = "darkgreen", Putamen = "orange", Cortex = "purple"),
    Sex = c(Male = "lightblue", Female = "pink")
  )
  
  # Create heatmap
  p <- pheatmap(heatmap_data_scaled,
                clustering_distance_rows = "euclidean",
                clustering_distance_cols = "euclidean",
                clustering_method = "ward.D2",
                annotation_col = annotation_col,
                annotation_colors = ann_colors,
                show_rownames = TRUE,
                show_colnames = FALSE,
                scale = "none",
                color = colorRampPalette(c("navy", "white", "red"))(100),
                main = "Top 50 Differentially Expressed Proteins")
  
  return(p)
}

# ==============================================================================
# 5. SEX-STRATIFIED ANALYSIS
# ==============================================================================

perform_sex_stratified_analysis <- function(data, metadata) {
  # Split by sex
  males <- metadata$Sex == "Male"
  females <- metadata$Sex == "Female"
  
  # Male-specific analysis
  male_results <- perform_differential_expression(
    data[, males],
    metadata[males, ],
    "PD_vs_Control"
  )
  
  # Female-specific analysis
  female_results <- perform_differential_expression(
    data[, females],
    metadata[females, ],
    "PD_vs_Control"
  )
  
  # Find sex-specific DEPs
  male_specific <- setdiff(
    rownames(male_results[male_results$adj.P.Val < 0.05, ]),
    rownames(female_results[female_results$adj.P.Val < 0.05, ])
  )
  
  female_specific <- setdiff(
    rownames(female_results[female_results$adj.P.Val < 0.05, ]),
    rownames(male_results[male_results$adj.P.Val < 0.05, ])
  )
  
  shared <- intersect(
    rownames(male_results[male_results$adj.P.Val < 0.05, ]),
    rownames(female_results[female_results$adj.P.Val < 0.05, ])
  )
  
  return(list(
    male_results = male_results,
    female_results = female_results,
    male_specific = male_specific,
    female_specific = female_specific,
    shared = shared
  ))
}

# ==============================================================================
# 6. MAIN ANALYSIS PIPELINE
# ==============================================================================

main <- function() {
  # Load data for each cohort
  cohorts <- c("PPMI", "AMP-PD", "GNPC")
  all_results <- list()
  
  for (cohort in cohorts) {
    cat(paste("Processing", cohort, "cohort...\n"))
    
    # Load data
    cohort_data <- load_proteomic_data(cohort)
    
    # Quality control
    qc_data <- perform_qc(cohort_data$data, cohort_data$metadata)
    
    # Batch correction
    corrected_data <- correct_batch_effects(qc_data$data, qc_data$metadata)
    
    # Differential expression analysis
    dep_results <- perform_differential_expression(
      corrected_data,
      qc_data$metadata,
      "PD_vs_Control"
    )
    
    # Pathway enrichment
    pathways <- perform_pathway_enrichment(dep_results)
    
    # Create visualizations
    volcano <- create_volcano_plot(dep_results, paste(cohort, "- PD vs Control"))
    heatmap <- create_heatmap(corrected_data, dep_results, qc_data$metadata)
    
    # Sex-stratified analysis
    sex_results <- perform_sex_stratified_analysis(corrected_data, qc_data$metadata)
    
    # Store results
    all_results[[cohort]] <- list(
      dep_results = dep_results,
      pathways = pathways,
      plots = list(volcano = volcano, heatmap = heatmap),
      sex_results = sex_results
    )
    
    # Save results
    write.csv(dep_results, 
              paste0("results/", cohort, "_differential_expression.csv"),
              row.names = TRUE)
    
    ggsave(paste0("figures/", cohort, "_volcano.pdf"), 
           volcano, 
           width = 8, 
           height = 6)
  }
  
  # Meta-analysis across cohorts
  cat("Performing meta-analysis across cohorts...\n")
  
  # Combine results
  combined_pvalues <- do.call(cbind, lapply(all_results, function(x) {
    x$dep_results$P.Value
  }))
  
  # Fisher's method for combining p-values
  library(metap)
  combined_p <- apply(combined_pvalues, 1, function(x) {
    sumlog(x)$p
  })
  
  # Adjust for multiple testing
  combined_fdr <- p.adjust(combined_p, method = "BH")
  
  cat("Analysis complete!\n")
  return(all_results)
}

# Run main analysis
results <- main()
