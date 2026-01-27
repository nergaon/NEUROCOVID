# -*- coding: utf-8 -*-
"""
the program:
creates qc plots
"""

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import numpy as np
import sys
sys.path.insert(0, "Z:/Analysis/general_scripts")
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from matplotlib.lines import Line2D
from scipy.stats import zscore

def plot_hist(df, result_dict, output_figure):
    fig, ax = plt.subplots(figsize=(20,10))
    plt.xlabel('Log2(expression)', fontsize=32)
    plt.ylabel('Number of Genes', fontsize=32)
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)
    
    legend_added = {}  # Track which labels have been added
    
    for col in df.columns:
        plt.ylim(0, 4000)  
        plt.xlim(0, 18)    
        hist, bins = np.histogram(df[col], bins=30)
        color = result_dict[col]
        
        # Extract label by removing the trailing number (last part after split by '_')
        label_parts = col.rsplit('_', 1)  # Split from the right to remove the last number
        label = label_parts[0]
        
        # Only add label to legend if this prefix hasn't been added yet
        plot_label = label if label not in legend_added else None
        legend_added[label] = True
        
        ax.plot(bins[:-1], hist, linewidth=2, color=color, label=plot_label)
    
    ax.legend(fontsize=20)
    plt.savefig(output_figure, bbox_inches='tight')
    plt.close() 
    return

def create_dictionary(samples):
    result_dict = {}
    for value in samples:
        if "Cortical_neurons_SARS" in value:
            result_dict[value] = 'black'
        if "Cortical_neurons_microglia_mock" in value:
            result_dict[value] = 'red'
        if "Cortical_neurons_microglia_SARS" in value:
            result_dict[value] = 'pink'
        if "Cortical_neurons_mock" in value:
            result_dict[value] = 'lightblue'
        if "Microglia_mock" in value:
            result_dict[value] = 'green'
        if "Microglia_SARS" in value:
            result_dict[value] = 'purple'
        if "organoids_mock" in value:
            result_dict[value] = 'orange'
        if "organoids_SARS" in value:
            result_dict[value] = 'brown'
    return result_dict

def plot_scatter(df, samples, in_dir):
    #creates pairwise scatter plots comparing expression levels between samples
    #repeats should be scattered around the diagonal
    fig = plt.figure(1)    
    for i in range(len(samples)-1):
        print(samples[i])
        for j in range(i+1, len(samples)):
            plt.subplots(figsize=(10,10))
            plt.scatter(df[df.columns[i]], df[df.columns[j]], color='dimgrey', s=2)
            print(str(i) + ' ' +str(j))
            plt.xlabel(df.columns[i],fontsize=40)
            plt.ylabel(df.columns[j],fontsize=40)
            plt.plot([0, 22], [0, 22], 'k-', color = 'black')
            fig_title = 'log2(counts)'
            plt.title(fig_title)
            output = in_dir + '/scatter/' + samples[i] + '_vs_' + samples[j] + '.jpg'
            plt.savefig(output)
            plt.close()       
    plt.tight_layout()
    plt.show()
    return()

def sort_table(df):
    # Define custom order for cell type and treatment combinations
    cell_treatment_order = [
        "Cortical_neurons_mock",
        "Cortical_neurons_SARS-CoV",
        "organoids_mock",
        "organoids_SARS-CoV",
        "Cortical_neurons_microglia_mock",
        "Cortical_neurons_microglia_SARS-CoV",
        "Microglia_mock",
        "Microglia_SARS-CoV"]
    sorted_cols = []
    for sample in cell_treatment_order:
        matching_cols = [col for col in df.columns if sample in col]
        sorted_cols.extend(matching_cols)
    df = df[sorted_cols]
    return df       

def select_rows(df, min_expression, min_samples, output_table, genes_table):
    #get only express genes and add gene data
    df = df.clip(lower=min_expression) #set all values that are less than min_expression to be min_expression
    threshold = int(min_samples * df.shape[1])  # 10% of total columns with min of 3 samples
    if threshold < 3:
        threshold = 3
    selected_rows = df[(df > min_expression).sum(axis=1) >= threshold]
    genes = genes_table[['Gene stable ID', 'Gene description', 'Chromosome/scaffold name', 'Gene name']]
    selected_rows = selected_rows.merge(genes, left_index=True, right_on='Gene stable ID').set_index('Gene stable ID')
    selected_rows.to_csv(output_table, sep="\t", index=True)  # Tab-separated
    selected_rows_only_numbers = selected_rows.iloc[:,:-3]
    return(selected_rows_only_numbers)

def plot_correlation(df, corr_figure):
    # Compute correlation matrix
    corr = df.corr()
    # Extract metadata from column names
    # Parse cell type and treatment from sample names
    cell_types = []
    treatments = []
    for col in df.columns:
        if "Cortical_neurons_microglia" in col:
            cell_type = "Cortical_neurons_microglia"
        elif "Cortical_neurons" in col:
            cell_type = "Cortical_neurons"
        elif "Microglia" in col:
            cell_type = "Microglia"
        elif "organoids" in col:
            cell_type = "organoids"
        else:
            cell_type = "unknown"     
        if "mock" in col:
            treatment = "mock"
        elif "SARS-CoV" in col or "SARS_CoV" in col:
            treatment = "SARS-CoV"
        else:
            treatment = "unknown"       
        cell_types.append(cell_type)
        treatments.append(treatment)
    # Define color mappings
    cell_colors_map = {
        "Cortical_neurons": "#1f77b4",      # blue
        "Cortical_neurons_microglia": "#ff7f0e",  # orange
        "organoids": "#d62728",              # red
        "Microglia": "#2ca02c"             # green
    }    
    treatment_colors_map = {
        "mock": "#ffffff",                  # white (no fill)
        "SARS-CoV": "#000000"               # black (fill)
    }
    # Create color bars
    col_colors = pd.DataFrame({
        "Cell Type": [cell_colors_map.get(ct, "gray") for ct in cell_types],
        "Treatment": [treatment_colors_map.get(tr, "gray") for tr in treatments]
    }, index=df.columns)
    
    # Create side color bars (use the same mapping for row_colors)
    row_colors = col_colors
    # Plot correlation matrix with color bars on top and side
    g = sns.clustermap(
        corr, cmap="coolwarm", center=0, vmin = -1, vmax = 1,
        col_cluster=False, row_cluster=False,  # Keep original column order
        col_colors=col_colors, row_colors=row_colors,
        figsize=(14, 12),
        cbar_pos=(0.02, 0.8, 0.03, 0.18))
    
    # Remove sample names (axis labels) from x and y axes
    g.ax_heatmap.set_xticklabels([])
    g.ax_heatmap.set_yticklabels([])
    
    # Create custom legend
    cell_legend_1 = Line2D([0], [0], color="#1f77b4", lw=6, label="Cortical_neurons")
    cell_legend_2 = Line2D([0], [0], color="#ff7f0e", lw=6, label="Cortical_neurons_microglia")
    cell_legend_3 = Line2D([0], [0], color="#d62728", lw=6, label="organoids")
    cell_legend_4 = Line2D([0], [0], color="#2ca02c", lw=6, label="Microglia")
    
    treatment_legend_1 = Line2D([0], [0], color="white", lw=6, label="mock", markeredgecolor="black", marker='s', markersize=10)
    treatment_legend_2 = Line2D([0], [0], color="black", lw=6, label="SARS-CoV", marker='s', markersize=10)
    
    # Add legend to the plot
    plt.legend(handles=[cell_legend_1, cell_legend_2, cell_legend_3, cell_legend_4, 
                       treatment_legend_1, treatment_legend_2],
               loc='center', bbox_to_anchor=(1.05, -0.8), fontsize=12,
               title="Cell Type / Treatment")
    
    # Save figure
    plt.savefig(corr_figure, bbox_inches='tight')
    plt.close() 
    return

def plot_pca(df, in_dir, condition):
    """
    Perform PCA on the expression table and plot the first two principal components
    colored by cell type and treatment.
    """
    # Step 1: Transpose and standardize the DataFrame
    df_transposed = df.T
    scaler = StandardScaler()
    df_scaled = scaler.fit_transform(df_transposed)
    
    # Step 2: Perform PCA
    pca = PCA(n_components=4)
    pca_result = pca.fit_transform(df_scaled)
    
    # Step 3: Get explained variance ratio
    explained_variance = pca.explained_variance_ratio_ * 100
    
    # Step 4: Create PCA results DataFrame
    pca_df = pd.DataFrame(pca_result, columns=['PCA1', 'PCA2', 'PCA3', 'PCA4'], index=df.columns)
    
    # Step 5: Extract cell type and treatment from sample names
    cell_types = []
    treatments = []
    
    for col in df.columns:
        # Extract cell type (check longest strings first)
        if "Cortical_neurons_microglia" in col:
            cell_type = "Cortical_neurons_microglia"
        elif "Cortical_neurons" in col:
            cell_type = "Cortical_neurons"
        elif "Microglia" in col:
            cell_type = "Microglia"
        elif "organoids" in col:
            cell_type = "organoids"
        else:
            cell_type = "unknown"
        
        # Extract treatment
        if "mock" in col:
            treatment = "mock"
        elif "SARS-CoV" in col or "SARS_CoV" in col:
            treatment = "SARS-CoV"
        else:
            treatment = "unknown"
        
        cell_types.append(cell_type)
        treatments.append(treatment)
    
    pca_df['cell'] = cell_types
    pca_df['treatment'] = treatments
    
    # Step 6: Define color and marker mappings
    cell_colors_map = {
        "Cortical_neurons": "#1f77b4",      # blue
        "Cortical_neurons_microglia": "#ff7f0e",  # orange
        "Microglia": "#2ca02c",             # green
        "organoids": "#d62728"              # red
    }    
    treatment_markers = {
        "mock": "o",        # circle
        "SARS-CoV": "s"     # square
    }
    
    # Step 7: Plot PCA results
    pairs = [(1,2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    for (pc_x, pc_y) in pairs:
        plt.figure(figsize=(10, 8))
        pca_component_1 = "PCA" + str(pc_x)
        pca_component_2 = "PCA" + str(pc_y)
        
        # Plot each combination of cell type and treatment
        for cell_type in sorted(pca_df['cell'].unique()):
            for treatment in sorted(pca_df['treatment'].unique()):
                subset = pca_df[(pca_df['cell'] == cell_type) & (pca_df['treatment'] == treatment)]                
                if len(subset) > 0:
                    color = cell_colors_map.get(cell_type, "gray")
                    marker = treatment_markers.get(treatment, "o")
                    label = f"{cell_type} - {treatment}"                    
                    plt.scatter(subset[pca_component_1], subset[pca_component_2],
                               color=color, marker=marker, label=label, s=100, alpha=0.7, edgecolors='black')
        
        # Add labels and legend
        plt.xlabel(f'{pca_component_1} ({explained_variance[pc_x - 1]:.2f}%)', fontsize=14)
        plt.ylabel(f'{pca_component_2} ({explained_variance[pc_y - 1]:.2f}%)', fontsize=14)
        plt.title(f'PCA{pc_x} vs PCA{pc_y}', fontsize=16)
        plt.legend(title="Cell Type - Treatment", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)
        plt.tight_layout()
        
        # Save figure
        if condition == "all":
            pca_figure = in_dir  + pca_component_1 + '_' + pca_component_2 + '.jpg'
        else:
            pca_figure = in_dir + condition + "_" + pca_component_1 + '_' + pca_component_2 + '.jpg'
        
        plt.savefig(pca_figure, bbox_inches='tight')
        plt.close()
    
    return

def plot_heatmap_genes(df, output_heatmap, output_table):
    df_filtered = df.apply(zscore, axis=1)  # Standardize across genes (rows)
    
    # Extract cell type and treatment from column names
    cell_types = []
    treatments = []
    
    for col in df.columns:
        # Extract cell type (check longest strings first)
        if "Cortical_neurons_microglia" in col:
            cell_type = "Cortical_neurons_microglia"
        elif "Cortical_neurons" in col:
            cell_type = "Cortical_neurons"
        elif "Microglia" in col:
            cell_type = "Microglia"
        elif "organoids" in col:
            cell_type = "organoids"
        else:
            cell_type = "unknown"
        
        # Extract treatment
        if "mock" in col:
            treatment = "mock"
        elif "SARS-CoV" in col or "SARS_CoV" in col:
            treatment = "SARS-CoV"
        else:
            treatment = "unknown"
        
        cell_types.append(cell_type)
        treatments.append(treatment)
    
    # Define color mappings (same as plot_correlation)
    cell_colors_map = {
        "Cortical_neurons": "#1f77b4",      # blue
        "Cortical_neurons_microglia": "#ff7f0e",  # orange
        "Microglia": "#2ca02c",             # green
        "organoids": "#d62728"              # red
    }
    
    treatment_colors_map = {
        "mock": "#ffffff",                  # white (no fill)
        "SARS-CoV": "#000000"               # black (fill)
    }
    
    # Create color bars (aligned with column order)
    col_colors = pd.DataFrame({
        "Cell Type": [cell_colors_map.get(ct, "gray") for ct in cell_types],
        "Treatment": [treatment_colors_map.get(tr, "gray") for tr in treatments]
    }, index=df.columns)
    
    # Plot matrix with color bars on top and side
    g = sns.clustermap(
        df_filtered, cmap="coolwarm", center=0, 
        col_cluster=False, row_cluster=True,  # Keep original column order
        col_colors=col_colors,
        figsize=(12, 10),
        cbar_pos=(0.02, 0.8, 0.03, 0.18))
    
    g.ax_heatmap.set_xticklabels([])  # Remove x-axis labels
    
    # Get the ordered gene list
    ordered_genes = df_filtered.index[g.dendrogram_row.reordered_ind]
    
    # Create custom legend
    cell_legend_1 = Line2D([0], [0], color="#1f77b4", lw=6, label="Cortical_neurons")
    cell_legend_2 = Line2D([0], [0], color="#ff7f0e", lw=6, label="Cortical_neurons_microglia")
    cell_legend_3 = Line2D([0], [0], color="#2ca02c", lw=6, label="Microglia")
    cell_legend_4 = Line2D([0], [0], color="#d62728", lw=6, label="organoids")
    
    treatment_legend_1 = Line2D([0], [0], color="white", lw=6, label="mock", markeredgecolor="black", marker='s', markersize=10)
    treatment_legend_2 = Line2D([0], [0], color="black", lw=6, label="SARS-CoV", marker='s', markersize=10)
    
    # Add legend to the plot
    plt.legend(handles=[cell_legend_1, cell_legend_2, cell_legend_3, cell_legend_4,
                       treatment_legend_1, treatment_legend_2],
               loc='center', bbox_to_anchor=(1.05, -0.8), fontsize=12,
               title="Cell Type / Treatment")
    
    plt.title("Heatmap of Standardized Log2 Expression Data")
    
    # Save figure
    plt.savefig(output_heatmap, bbox_inches='tight')
    plt.close()
    
    df_ordered = df_filtered.loc[ordered_genes]
    df_ordered.to_csv(output_table, sep="\t", index=True)
    return(df_ordered)

def main():
    min_samples = 0.1 #min % samples with less than min expression counts
    min_expression = 6 #min log2 expression value
    input_genes_data = '/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/expression_table/biomart_genes_details.txt'
    genes_table = pd.read_csv(input_genes_data, delimiter =  "\t")    
    
    in_dir = '/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/expression_table/'
    input_expression = in_dir + 'expression_table_HN1.txt'
    print(input_expression)
    counts_table = pd.read_csv(input_expression, delimiter =  "\s|\t", engine='python', skiprows=1)        
    counts_table = counts_table.reset_index().set_index('Geneid', drop=True)
    counts_table.drop(['index','Chr','Start','End','Strand','Length'], axis=1, inplace=True) 
    counts_table.rename(columns=lambda x: x.replace('.sort.bam', ''), inplace=True)
    counts_table_log = counts_table.replace(0, 1)
    counts_table_log = np.log2(counts_table_log)
    samples = counts_table.columns
    result_dict = create_dictionary(samples) #color of each sample
    output_figure = in_dir + 'distribution_samples.jpg'
    plot_hist(counts_table_log,result_dict, output_figure)
    
    #plot_scatter(counts_table_log, samples, in_dir)
    
    #create a table with only express genes
    sort_table_df = sort_table(counts_table_log)
    output_table = in_dir + 'log_express_gene.txt'
    sort_table_df.to_csv(output_table, sep="\t", index=True)
    #get only express genes and set all the values that are smaller than min_expression to be min_expression
    express_genes = select_rows(sort_table_df, min_expression, min_samples, output_table, genes_table)
    # Calculate standard deviation for each gene and select top N most variable genes
    gene_std = express_genes.std(axis=1)
    top_genes = gene_std.sort_values(ascending=False).head(1000).index
    # Subset the expression table to the top N genes
    expression_table_top = express_genes.loc[top_genes]
    #z normelaize the table
    #create correlation matrix of the 1000 genes that have the highest std
    corr_figure = in_dir + 'correlation.jpg'
    plot_correlation(expression_table_top, corr_figure)
    plot_pca(expression_table_top, in_dir, 'all')
    output_heatmap = in_dir + 'heatmap_top1000_all.jpg'
    output_table = in_dir + 'heatmap_top1000_all.txt'
    plot_heatmap_genes(expression_table_top, output_heatmap, output_table)
    
    #choose genes and do pca for each cell type seperatly. 
    #cells_type = ["MON", "pDCs"]
    #cell_dir = in_dir + "QC/cell/"
    #for cell in cells_type:
    #    df_cell = sort_table_df[sort_table_df.columns[sort_table_df.columns.str.contains(cell)]]
    #    output_table = cell_dir + cell + '_log_express_gene.txt'
    #    express_genes = select_rows(df_cell, min_expression, min_samples, output_table, genes_table)
        # Calculate standard deviation for each gene and select top N most variable genes
    #    gene_std = express_genes.std(axis=1)
    #    top_genes = gene_std.sort_values(ascending=False).head(1000).index
        # Subset the expression table to the top N genes
    #    expression_table_top = express_genes.loc[top_genes]
    #    plot_pca(expression_table_top, cell_dir, cell)
    #    output_heatmap = cell_dir + '/heatmap_top10000_' + cell + '.jpg'
    #    output_table = cell_dir + '/heatmap_top10000_' + cell + '.txt'
    #    plot_heatmap_genes(expression_table_top, output_heatmap, output_table)
    return

if __name__ == "__main__":
    main() 