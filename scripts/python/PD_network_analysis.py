#!/usr/bin/env python3
"""
================================================================================
Network Analysis for Parkinson's Disease Proteomics
================================================================================
This script performs co-expression network analysis using MEGENA-like algorithms
and identifies key driver proteins using Bayesian causal network inference
Author: Research Team
Date: 2025
================================================================================
"""

import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy import sparse
from scipy.cluster import hierarchy
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, DBSCAN
import warnings
warnings.filterwarnings('ignore')

# For parallel processing
from joblib import Parallel, delayed
import multiprocessing

# Set style for plots
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)

# ============================================================================
# 1. DATA LOADING AND PREPROCESSING
# ============================================================================

class ProteomicNetworkAnalysis:
    """
    Main class for proteomic network analysis following MEGENA approach
    """
    
    def __init__(self, expression_data, metadata, genetic_data=None):
        """
        Initialize network analysis
        
        Parameters:
        -----------
        expression_data : pd.DataFrame
            Protein expression matrix (proteins x samples)
        metadata : pd.DataFrame
            Sample metadata including disease status, demographics
        genetic_data : pd.DataFrame, optional
            Genetic data for eQTL analysis
        """
        self.expression_data = expression_data
        self.metadata = metadata
        self.genetic_data = genetic_data
        self.correlation_matrix = None
        self.adjacency_matrix = None
        self.modules = {}
        self.key_drivers = {}
        
    def calculate_correlation_matrix(self, method='pearson', min_samples=20):
        """
        Calculate correlation matrix between all proteins
        """
        print("Calculating correlation matrix...")
        
        if method == 'pearson':
            self.correlation_matrix = self.expression_data.T.corr(method='pearson')
        elif method == 'spearman':
            self.correlation_matrix = self.expression_data.T.corr(method='spearman')
        elif method == 'bicor':
            # Biweight midcorrelation (robust to outliers)
            self.correlation_matrix = self._calculate_bicor()
        
        return self.correlation_matrix
    
    def _calculate_bicor(self):
        """
        Calculate biweight midcorrelation (robust correlation measure)
        """
        from astropy.stats import biweight_midcorrelation
        
        n_proteins = self.expression_data.shape[0]
        bicor_matrix = np.zeros((n_proteins, n_proteins))
        
        for i in range(n_proteins):
            for j in range(i, n_proteins):
                if i == j:
                    bicor_matrix[i, j] = 1.0
                else:
                    corr = biweight_midcorrelation(
                        self.expression_data.iloc[i, :],
                        self.expression_data.iloc[j, :]
                    )
                    bicor_matrix[i, j] = corr
                    bicor_matrix[j, i] = corr
        
        return pd.DataFrame(
            bicor_matrix,
            index=self.expression_data.index,
            columns=self.expression_data.index
        )
    
    def construct_network(self, correlation_threshold=0.3, p_value_threshold=0.05):
        """
        Construct co-expression network using correlation threshold and significance
        """
        print("Constructing co-expression network...")
        
        # Calculate p-values for correlations
        n_samples = self.expression_data.shape[1]
        p_values = self._calculate_correlation_pvalues(n_samples)
        
        # Create adjacency matrix
        self.adjacency_matrix = (
            (np.abs(self.correlation_matrix) > correlation_threshold) &
            (p_values < p_value_threshold)
        ).astype(int)
        
        # Remove self-loops
        np.fill_diagonal(self.adjacency_matrix.values, 0)
        
        # Create NetworkX graph
        self.network = nx.from_pandas_adjacency(self.adjacency_matrix)
        
        print(f"Network contains {self.network.number_of_nodes()} nodes "
              f"and {self.network.number_of_edges()} edges")
        
        return self.network
    
    def _calculate_correlation_pvalues(self, n_samples):
        """
        Calculate p-values for correlation coefficients
        """
        # Fisher's z-transformation
        z_scores = np.arctanh(self.correlation_matrix)
        se = 1.0 / np.sqrt(n_samples - 3)
        z_critical = z_scores / se
        
        # Two-tailed p-values
        p_values = 2 * (1 - stats.norm.cdf(np.abs(z_critical)))
        
        return pd.DataFrame(
            p_values,
            index=self.correlation_matrix.index,
            columns=self.correlation_matrix.columns
        )
    
    def identify_modules_megena(self, min_module_size=15):
        """
        Identify co-expression modules using MEGENA-like approach
        (Multiscale Embedded Gene Co-Expression Network Analysis)
        """
        print("Identifying co-expression modules using MEGENA approach...")
        
        # Step 1: Planar Filtered Network (PFN) construction
        pfn = self._construct_planar_filtered_network()
        
        # Step 2: Multiscale clustering
        modules = self._multiscale_clustering(pfn, min_module_size)
        
        # Step 3: Module hierarchy construction
        self.module_hierarchy = self._construct_module_hierarchy(modules)
        
        # Step 4: Calculate module eigengenes
        self.module_eigengenes = self._calculate_module_eigengenes(modules)
        
        self.modules = modules
        print(f"Identified {len(modules)} modules")
        
        return modules
    
    def _construct_planar_filtered_network(self):
        """
        Construct Planar Filtered Network (PFN) following MEGENA approach
        """
        # Convert correlation matrix to distance matrix
        distance_matrix = 1 - np.abs(self.correlation_matrix)
        
        # Apply minimum spanning tree
        from scipy.sparse.csgraph import minimum_spanning_tree
        mst = minimum_spanning_tree(distance_matrix)
        
        # Add additional edges while maintaining planarity
        # This is a simplified version - full MEGENA uses more sophisticated approach
        pfn = nx.from_scipy_sparse_matrix(mst)
        
        return pfn
    
    def _multiscale_clustering(self, network, min_size):
        """
        Perform multiscale clustering on the network
        """
        # Use Louvain community detection at multiple resolutions
        import community as community_louvain
        
        modules = {}
        resolutions = [0.5, 1.0, 1.5, 2.0]  # Multiple resolution parameters
        
        for res in resolutions:
            partition = community_louvain.best_partition(
                network, 
                resolution=res,
                randomize=False
            )
            
            # Group nodes by module
            for node, module_id in partition.items():
                module_key = f"res_{res}_module_{module_id}"
                if module_key not in modules:
                    modules[module_key] = []
                modules[module_key].append(node)
        
        # Filter out small modules
        modules = {k: v for k, v in modules.items() if len(v) >= min_size}
        
        return modules
    
    def _construct_module_hierarchy(self, modules):
        """
        Construct hierarchical relationships between modules
        """
        hierarchy = {}
        
        # Calculate overlap between modules
        for m1 in modules:
            for m2 in modules:
                if m1 != m2:
                    overlap = len(set(modules[m1]) & set(modules[m2]))
                    if overlap > 0:
                        if m1 not in hierarchy:
                            hierarchy[m1] = {}
                        hierarchy[m1][m2] = overlap
        
        return hierarchy
    
    def _calculate_module_eigengenes(self, modules):
        """
        Calculate module eigengenes (first principal component of module expression)
        """
        eigengenes = {}
        
        for module_name, module_genes in modules.items():
            if len(module_genes) > 1:
                # Get expression data for module genes
                module_expr = self.expression_data.loc[module_genes, :]
                
                # Perform PCA
                pca = PCA(n_components=1)
                eigengene = pca.fit_transform(module_expr.T)
                
                eigengenes[module_name] = eigengene.flatten()
        
        return pd.DataFrame(eigengenes, index=self.expression_data.columns)
    
    def identify_key_drivers(self, module_name):
        """
        Identify key driver proteins within a module using network topology
        """
        print(f"Identifying key drivers for module {module_name}...")
        
        module_genes = self.modules[module_name]
        
        # Create subnetwork for module
        subnetwork = self.network.subgraph(module_genes)
        
        # Calculate various centrality measures
        degree_centrality = nx.degree_centrality(subnetwork)
        betweenness_centrality = nx.betweenness_centrality(subnetwork)
        eigenvector_centrality = nx.eigenvector_centrality(subnetwork, max_iter=1000)
        closeness_centrality = nx.closeness_centrality(subnetwork)
        
        # Combine centrality measures to identify key drivers
        driver_scores = {}
        for node in module_genes:
            if node in degree_centrality:
                driver_scores[node] = (
                    0.25 * degree_centrality[node] +
                    0.25 * betweenness_centrality[node] +
                    0.25 * eigenvector_centrality[node] +
                    0.25 * closeness_centrality[node]
                )
        
        # Sort by driver score
        key_drivers = sorted(driver_scores.items(), key=lambda x: x[1], reverse=True)
        
        return key_drivers

# ============================================================================
# 2. BAYESIAN CAUSAL NETWORK ANALYSIS
# ============================================================================

class BayesianCausalNetwork:
    """
    Construct Bayesian causal networks integrating genetic and expression data
    """
    
    def __init__(self, expression_data, genetic_data, metadata):
        """
        Initialize Bayesian network analysis
        """
        self.expression_data = expression_data
        self.genetic_data = genetic_data
        self.metadata = metadata
        self.eqtls = None
        self.causal_network = None
        
    def identify_eqtls(self, p_threshold=1e-5):
        """
        Identify expression quantitative trait loci (eQTLs)
        """
        print("Identifying eQTLs...")
        
        eqtls = []
        
        for protein in self.expression_data.index:
            protein_expr = self.expression_data.loc[protein, :]
            
            # Test association with each SNP
            for snp in self.genetic_data.columns:
                genotype = self.genetic_data[snp]
                
                # Linear regression
                from scipy import stats
                slope, intercept, r_value, p_value, std_err = stats.linregress(
                    genotype, protein_expr
                )
                
                if p_value < p_threshold:
                    eqtls.append({
                        'protein': protein,
                        'snp': snp,
                        'beta': slope,
                        'p_value': p_value,
                        'r_squared': r_value**2
                    })
        
        self.eqtls = pd.DataFrame(eqtls)
        print(f"Identified {len(self.eqtls)} significant eQTLs")
        
        return self.eqtls
    
    def construct_causal_network(self, method='pc'):
        """
        Construct causal network using PC algorithm or other methods
        """
        print("Constructing causal network...")
        
        if method == 'pc':
            # Use PC algorithm (Peter-Clark)
            from pgmpy.estimators import PC
            from pgmpy.models import BayesianNetwork
            
            # Prepare data
            data = pd.concat([
                self.expression_data.T,
                self.genetic_data
            ], axis=1)
            
            # Run PC algorithm
            pc = PC(data)
            model = pc.estimate(
                variant='stable',
                max_cond_vars=3,
                return_type='dag'
            )
            
            self.causal_network = model
            
        elif method == 'ges':
            # Use Greedy Equivalence Search
            from cdt.causality.graph import GES
            
            ges = GES()
            self.causal_network = ges.predict(self.expression_data.T)
            
        return self.causal_network
    
    def identify_causal_key_drivers(self, target_module_genes):
        """
        Identify causal key drivers for a set of genes
        """
        print("Identifying causal key drivers...")
        
        key_drivers = {}
        
        # For each potential driver
        for driver in self.causal_network.nodes():
            # Count number of targets it regulates
            targets_regulated = 0
            
            for target in target_module_genes:
                if self.causal_network.has_edge(driver, target):
                    targets_regulated += 1
            
            # Calculate driver score
            if targets_regulated > 0:
                driver_score = targets_regulated / len(target_module_genes)
                key_drivers[driver] = driver_score
        
        # Sort by driver score
        key_drivers = sorted(key_drivers.items(), key=lambda x: x[1], reverse=True)
        
        return key_drivers

# ============================================================================
# 3. MODULE PRESERVATION ANALYSIS
# ============================================================================

def assess_module_preservation(reference_modules, test_expression, n_permutations=100):
    """
    Assess preservation of modules in independent dataset
    """
    preservation_scores = {}
    
    for module_name, module_genes in reference_modules.items():
        # Filter to genes present in test data
        common_genes = list(set(module_genes) & set(test_expression.index))
        
        if len(common_genes) < 10:
            continue
        
        # Calculate module connectivity in reference
        ref_connectivity = calculate_module_connectivity(
            reference_expression[module_genes]
        )
        
        # Calculate module connectivity in test
        test_connectivity = calculate_module_connectivity(
            test_expression[common_genes]
        )
        
        # Permutation test
        null_distribution = []
        for _ in range(n_permutations):
            random_genes = np.random.choice(
                test_expression.index, 
                size=len(common_genes),
                replace=False
            )
            null_connectivity = calculate_module_connectivity(
                test_expression[random_genes]
            )
            null_distribution.append(null_connectivity)
        
        # Calculate Z-score
        z_score = (test_connectivity - np.mean(null_distribution)) / np.std(null_distribution)
        
        preservation_scores[module_name] = {
            'z_score': z_score,
            'p_value': 2 * (1 - stats.norm.cdf(abs(z_score))),
            'n_genes': len(common_genes)
        }
    
    return preservation_scores

def calculate_module_connectivity(expression_matrix):
    """
    Calculate average connectivity within a module
    """
    corr_matrix = expression_matrix.T.corr()
    # Remove diagonal
    np.fill_diagonal(corr_matrix.values, np.nan)
    # Calculate mean absolute correlation
    connectivity = np.nanmean(np.abs(corr_matrix.values))
    
    return connectivity

# ============================================================================
# 4. VISUALIZATION FUNCTIONS
# ============================================================================

def visualize_network_modules(network, modules, output_file='network_modules.pdf'):
    """
    Visualize network with module coloring
    """
    plt.figure(figsize=(15, 12))
    
    # Create layout
    pos = nx.spring_layout(network, k=0.5, iterations=50)
    
    # Create color map for modules
    colors = plt.cm.Set3(np.linspace(0, 1, len(modules)))
    node_colors = {}
    
    for i, (module_name, module_genes) in enumerate(modules.items()):
        for gene in module_genes:
            if gene in network.nodes():
                node_colors[gene] = colors[i]
    
    # Draw network
    node_color_list = [node_colors.get(node, 'gray') for node in network.nodes()]
    
    nx.draw_networkx_nodes(
        network, pos,
        node_color=node_color_list,
        node_size=50,
        alpha=0.8
    )
    
    nx.draw_networkx_edges(
        network, pos,
        alpha=0.2,
        edge_color='gray'
    )
    
    plt.title('Protein Co-expression Network Modules', fontsize=16, fontweight='bold')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()

def plot_module_trait_correlation(module_eigengenes, traits, output_file='module_trait_correlation.pdf'):
    """
    Plot correlation between module eigengenes and clinical traits
    """
    # Calculate correlations
    correlations = pd.DataFrame()
    p_values = pd.DataFrame()
    
    for module in module_eigengenes.columns:
        for trait in traits.columns:
            corr, pval = stats.pearsonr(
                module_eigengenes[module],
                traits[trait]
            )
            correlations.loc[module, trait] = corr
            p_values.loc[module, trait] = pval
    
    # Create heatmap
    plt.figure(figsize=(10, 12))
    
    # Create annotation with correlation values and significance stars
    annot_data = correlations.copy()
    for i in range(len(correlations)):
        for j in range(len(correlations.columns)):
            if p_values.iloc[i, j] < 0.001:
                annot_data.iloc[i, j] = f"{correlations.iloc[i, j]:.2f}***"
            elif p_values.iloc[i, j] < 0.01:
                annot_data.iloc[i, j] = f"{correlations.iloc[i, j]:.2f}**"
            elif p_values.iloc[i, j] < 0.05:
                annot_data.iloc[i, j] = f"{correlations.iloc[i, j]:.2f}*"
            else:
                annot_data.iloc[i, j] = f"{correlations.iloc[i, j]:.2f}"
    
    sns.heatmap(
        correlations,
        annot=annot_data,
        fmt='',
        cmap='RdBu_r',
        center=0,
        vmin=-1,
        vmax=1,
        cbar_kws={'label': 'Correlation'},
        linewidths=0.5
    )
    
    plt.title('Module-Trait Correlations', fontsize=14, fontweight='bold')
    plt.xlabel('Clinical Traits', fontsize=12)
    plt.ylabel('Modules', fontsize=12)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()

# ============================================================================
# 5. MAIN ANALYSIS PIPELINE
# ============================================================================

def main():
    """
    Main analysis pipeline
    """
    print("="*80)
    print("PARKINSON'S DISEASE PROTEOMIC NETWORK ANALYSIS")
    print("="*80)
    
    # Load data (example paths - replace with actual data)
    print("\n1. Loading data...")
    expression_data = pd.read_csv('data/PD_proteomics_expression.csv', index_col=0)
    metadata = pd.read_csv('data/PD_metadata.csv')
    genetic_data = pd.read_csv('data/PD_genotypes.csv', index_col=0)
    
    # Initialize network analysis
    print("\n2. Initializing network analysis...")
    network_analysis = ProteomicNetworkAnalysis(
        expression_data, 
        metadata, 
        genetic_data
    )
    
    # Calculate correlations
    print("\n3. Calculating correlation matrix...")
    correlation_matrix = network_analysis.calculate_correlation_matrix(method='bicor')
    
    # Construct network
    print("\n4. Constructing co-expression network...")
    network = network_analysis.construct_network(
        correlation_threshold=0.3,
        p_value_threshold=0.05
    )
    
    # Identify modules
    print("\n5. Identifying co-expression modules...")
    modules = network_analysis.identify_modules_megena(min_module_size=15)
    
    # Identify key drivers for top modules
    print("\n6. Identifying key driver proteins...")
    top_modules = list(modules.keys())[:10]  # Analyze top 10 modules
    
    all_key_drivers = {}
    for module in top_modules:
        key_drivers = network_analysis.identify_key_drivers(module)
        all_key_drivers[module] = key_drivers[:20]  # Top 20 drivers per module
    
    # Bayesian causal network analysis
    print("\n7. Constructing Bayesian causal network...")
    causal_analysis = BayesianCausalNetwork(
        expression_data,
        genetic_data,
        metadata
    )
    
    # Identify eQTLs
    eqtls = causal_analysis.identify_eqtls(p_threshold=1e-5)
    
    # Construct causal network
    causal_network = causal_analysis.construct_causal_network(method='pc')
    
    # Identify causal key drivers
    print("\n8. Identifying causal key drivers...")
    causal_drivers = {}
    for module_name, module_genes in modules.items():
        drivers = causal_analysis.identify_causal_key_drivers(module_genes)
        causal_drivers[module_name] = drivers[:10]  # Top 10 causal drivers
    
    # Visualizations
    print("\n9. Creating visualizations...")
    
    # Network visualization
    visualize_network_modules(network, modules, 'figures/network_modules.pdf')
    
    # Module-trait correlations
    clinical_traits = metadata[['UPDRS_motor', 'Braak_stage', 'Cognitive_score']]
    plot_module_trait_correlation(
        network_analysis.module_eigengenes,
        clinical_traits,
        'figures/module_trait_correlation.pdf'
    )
    
    # Save results
    print("\n10. Saving results...")
    
    # Save modules
    with open('results/modules.txt', 'w') as f:
        for module_name, module_genes in modules.items():
            f.write(f">{module_name}\n")
            f.write(','.join(module_genes) + '\n')
    
    # Save key drivers
    pd.DataFrame(all_key_drivers).to_csv('results/key_drivers.csv')
    
    # Save eQTLs
    eqtls.to_csv('results/eqtls.csv', index=False)
    
    print("\nAnalysis complete!")
    print("="*80)

if __name__ == "__main__":
    main()
