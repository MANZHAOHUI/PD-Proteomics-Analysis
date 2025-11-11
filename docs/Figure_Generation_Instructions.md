# Figure Generation Instructions for PD Proteomics Paper
## Software Requirements and Detailed Instructions

This document provides comprehensive instructions for generating all figures presented in the Parkinson's Disease proteomics manuscript.

---

## Figure 1: Study Overview and Differential Expression Analysis

### Required Software:
- **R/RStudio** (v4.0+) for statistical analysis
- **GraphPad Prism** (v9.0+) for publication-quality plots
- **Adobe Illustrator** or **Inkscape** for figure assembly

### Panel A: Study Design Flowchart
**Software:** Adobe Illustrator / BioRender

**Steps:**
1. Open BioRender (www.biorender.com) or Adobe Illustrator
2. Create workflow diagram showing:
   - Three cohorts (PPMI, AMP-PD, GNPC) with sample numbers
   - Brain regions (SN, Putamen, Frontal Cortex) 
   - Processing pipeline: Sample collection → Protein extraction → TMT labeling → LC-MS/MS → Data analysis
3. Use consistent color scheme: PD samples (red), Controls (blue)
4. Export as vector format (SVG/PDF) at 300 DPI

### Panel B: Volcano Plot
**Software:** R with EnhancedVolcano package

**R Code:**
```r
library(EnhancedVolcano)
EnhancedVolcano(results,
    lab = results$Gene_Symbol,
    x = 'logFC',
    y = 'adj.P.Val',
    title = 'PD vs Control',
    pCutoff = 0.05,
    FCcutoff = 0.25,
    pointSize = 2.0,
    labSize = 3.0,
    colAlpha = 0.5,
    legendPosition = 'right',
    legendLabSize = 12,
    legendIconSize = 4.0)
```

### Panel C: Pathway Enrichment Bar Plot
**Software:** R with ggplot2

**R Code:**
```r
library(ggplot2)
pathway_plot <- ggplot(pathway_data, aes(x = reorder(Description, -log10(p.adjust)), 
                                         y = -log10(p.adjust),
                                         fill = Category)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_classic() +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "", y = "-Log10(FDR)", title = "Enriched Pathways")
```

---

## Figure 2: Co-expression Network Analysis

### Panel A: Network Overview
**Software:** Cytoscape (v3.9+)

**Steps:**
1. Import network file (edges.csv and nodes.csv)
2. Apply force-directed layout (Layout → Prefuse Force Directed Layout)
3. Parameters:
   - Default Spring Coefficient: 1e-4
   - Default Spring Length: 50
   - Default Node Mass: 3
4. Color nodes by module membership
5. Size nodes by degree centrality
6. Export as high-resolution PDF

### Panel B: Module Hierarchy (Sunburst Plot)
**Software:** R with sunburstR package or Python with plotly

**Python Code:**
```python
import plotly.graph_objects as go

fig = go.Figure(go.Sunburst(
    labels=module_labels,
    parents=module_parents,
    values=module_sizes,
    branchvalues="total",
    marker=dict(colors=module_colors)
))
fig.update_layout(margin=dict(t=0, l=0, r=0, b=0))
fig.show()
```

### Panel C-E: Module Subnetworks
**Software:** Gephi (v0.9.2+)

**Steps:**
1. Import module subnetwork
2. Apply ForceAtlas2 layout algorithm:
   - Scaling: 10.0
   - Gravity: 1.0
   - Edge Weight Influence: 1.0
3. Color scheme:
   - Downregulated proteins: Blue
   - Upregulated proteins: Red
   - Unchanged: Gray
4. Node size: Proportional to degree
5. Edge thickness: Proportional to correlation strength
6. Export using Sigma.js exporter for interactive version

---

## Figure 3: Bayesian Causal Network

### Software:** R with bnlearn and Rgraphviz packages

**R Code:**
```r
library(bnlearn)
library(Rgraphviz)

# Construct network
network <- pc.stable(data, alpha = 0.05)

# Customize appearance
graph_attrs <- list(rankdir = "LR", bgcolor = "white")
node_attrs <- list(shape = "ellipse", fontsize = 12)
edge_attrs <- list(color = "gray50")

# Plot
graphviz.plot(network, 
              graph = graph_attrs,
              node = node_attrs, 
              edge = edge_attrs)
```

---

## Figure 4: Experimental Validation

### Panel A: Experimental Design Schematic
**Software:** BioRender

**Template:** Cell Culture & Co-culture templates
**Elements to include:**
- iPSC differentiation timeline
- Dopaminergic neurons and astrocytes
- CRISPR knockdown representation
- Co-culture setup

### Panels B-E: Bar Plots and Statistical Comparisons
**Software:** GraphPad Prism

**Steps in Prism:**
1. Enter data in Column format
2. Choose Column analysis → t-test or One-way ANOVA
3. Graph type: Column bar graph with scatter points
4. Error bars: SEM
5. Statistical annotations: 
   - *p < 0.05
   - **p < 0.01  
   - ***p < 0.001
6. Export as vector format

### Panel F: Heatmap of Rescued Proteins
**Software:** R with ComplexHeatmap

**R Code:**
```r
library(ComplexHeatmap)
Heatmap(expression_matrix,
        name = "Z-score",
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        show_row_names = FALSE,
        show_column_names = TRUE,
        row_split = protein_categories,
        column_split = sample_groups,
        top_annotation = HeatmapAnnotation(
            Treatment = anno_simple(treatment_groups)
        ))
```

---

## Figure 5: Sex-Specific Analysis

### Panel A: Venn Diagram
**Software:** R with VennDiagram or online tool

**R Code:**
```r
library(VennDiagram)
venn.diagram(
    x = list(Male = male_proteins, Female = female_proteins),
    filename = "sex_venn.png",
    col = c("lightblue", "pink"),
    fill = c("lightblue", "pink"),
    alpha = 0.5,
    cex = 1.5,
    fontfamily = "sans",
    cat.fontfamily = "sans"
)
```

### Panels B-C: Pathway Enrichment Comparison
**Software:** R with clusterProfiler

**Visualization approach:**
1. Create dotplot for male-specific pathways
2. Create dotplot for female-specific pathways
3. Combine using patchwork or cowplot
4. Use consistent color scale for p-values

---

## Supplementary Figures

### Multi-panel Assembly
**Software:** Adobe Illustrator or Inkscape

**Best Practices:**
1. Import all panels as linked files
2. Maintain consistent:
   - Font (Arial or Helvetica)
   - Font sizes (Title: 14pt, Axis labels: 12pt, Tick labels: 10pt)
   - Line weights (1pt for axes, 0.5pt for minor elements)
3. Add panel labels (A, B, C) in bold 16pt font
4. Align panels using smart guides
5. Export final figure at 300 DPI for print, 150 DPI for web

---

## Software-Specific Settings

### GraphPad Prism Settings:
- Graph size: 3.5 inches (single column) or 7 inches (double column)
- Font: Arial
- Line thickness: 1 point
- Symbol size: 8 points
- Export: PDF or EPS format

### Cytoscape Settings:
- Canvas size: 2048 x 2048 pixels minimum
- Node border width: 2.0
- Edge opacity: 50%
- Background: White
- Export: PDF with text as font

### R Plot Export Settings:
```r
ggsave("figure.pdf", 
       width = 7, 
       height = 5, 
       units = "in", 
       dpi = 300)
```

### Python Plot Export Settings:
```python
plt.savefig("figure.pdf", 
            dpi=300, 
            bbox_inches='tight',
            format='pdf')
```

---

## Color Schemes

### Disease Status:
- Control: #0066CC (Blue)
- PD: #CC0000 (Red)

### Brain Regions:
- Substantia Nigra: #2E7D32 (Dark Green)
- Putamen: #F57C00 (Orange)  
- Frontal Cortex: #7B1FA2 (Purple)

### Expression Changes:
- Downregulated: #1976D2 (Blue)
- Upregulated: #D32F2F (Red)
- Unchanged: #757575 (Gray)

### Sex:
- Male: #64B5F6 (Light Blue)
- Female: #F48FB1 (Pink)

---

## Data Availability for Figure Recreation

All raw data files required to generate figures are available at:
- Proteomics data: www.pd-proteome-network.org/data
- Network files: www.pd-proteome-network.org/networks
- Analysis scripts: github.com/pd-proteomics/analysis

For questions about figure generation, contact: 
analysis-support@pd-proteome-network.org
