import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from venn import venn
import matplotlib.pyplot as plt
import json
import seaborn as sns
import matplotlib.patches as mpatches
import gseapy as gp
from gprofiler import GProfiler
from matplotlib.lines import Line2D
from plot_config import colors
from plot_config import sort_table

def differential_expression_deseq2(df, conditions, mock_label, virus_label, fc_threshold, pval_threshold, min_expression, out_dir=None):
    """
    Differential expression using DESeq2 (via pydeseq2)
    Returns:
        up_genes: dict of genes up in virus per condition
        down_genes: dict of genes up in mock per condition
        results_df: dataframe of all DE results
    """
    all_results = []
    up_genes = {}
    down_genes = {}
    print("size of df:", df.shape)
    rows = []
    rows.append({"condition": "All", "n_genes": df.shape[0], "SARS-CoV-2_genes": np.nan, "mock_genes": np.nan})
    for cond in conditions:
        cond_cols = [
            c for c in df.columns
            if get_condition_from_col(c, mock_label, virus_label) == cond
        ]
        if len(cond_cols) == 0:
            print(f"No columns for {cond}, skipping")
            continue
        mock_cols = [c for c in cond_cols if mock_label in c]
        virus_cols = [c for c in cond_cols if virus_label in c]
        print(cond, mock_cols, virus_cols)
        min_samples = min(len(mock_cols), len(virus_cols))  #gene is express in all repeats
        if len(mock_cols) == 0 or len(virus_cols) == 0:
            print(f"No mock/virus columns for {cond}, skipping")
            continue
        count_data = df[cond_cols]
        # ---- Expression filter (DESeq2-style prefiltering) ----
        mask = (
            (count_data[mock_cols] > min_expression).sum(axis=1) >= min_samples
        ) | (
            (count_data[virus_cols] > min_expression).sum(axis=1) >= min_samples
        )
        count_data = count_data.loc[mask]
        print(f"{cond}: {count_data.shape[0]} genes pass expression filter")
        background_genes = set(count_data.index)
        # ---- Sample metadata ----
        coldata = pd.DataFrame(index=count_data.columns)
        coldata["condition"] = [
            "virus" if virus_label in c else "mock"
            for c in count_data.columns]
        # ---- DESeq2 ----
        dds = DeseqDataSet(counts=count_data.T, metadata=coldata, design="~ condition", refit_cooks=True)
        dds.deseq2()
        stats = DeseqStats(dds, contrast=("condition", "virus", "mock"))
        stats.summary()
        res = stats.results_df.copy()
        res["Gene"] = res.index
        res["Condition"] = cond
        # Rename to match your downstream expectations
        res = res.rename(columns={
            "log2FoldChange": "Log2FC",
            "padj": "FDR",
            "pvalue": "PValue"})
        all_results.append(res)
        # ---- Significant genes ----
        sig = res[
            (res["FDR"] < pval_threshold) &
            (abs(res["Log2FC"]) > np.log2(fc_threshold))]
        up_genes[cond] = set(sig[sig["Log2FC"] > 0]["Gene"])
        down_genes[cond] = set(sig[sig["Log2FC"] < 0]["Gene"])
        rows.append({"condition": cond, "n_genes": count_data.shape[0], 
                     "SARS-CoV-2_genes": len(up_genes[cond]), "mock_genes": len(down_genes[cond])})
        enrichmen(background_genes, up_genes[cond], out_dir, cond, "up_SARS-CoV")
        enrichmen(background_genes, down_genes[cond], out_dir, cond, "up_mock")
        #plot heatmap
        cell_colors_map, treatment_colors_map, treatment_markers, cell_types, treatments = colors(df.columns)
        plot_heatmap(virus_label, up_genes[cond], down_genes[cond], count_data, cond, out_dir, cell_colors_map, treatment_colors_map, cell_types, treatments, df)

    results_df = pd.concat(all_results, ignore_index=True)
    sum_table = pd.DataFrame(rows)
    print(sum_table)
    return up_genes, down_genes, results_df

def arrangeDS(count_data_df, samples):
    conditions = ["Cortical_neurons_microglia", "Cortical_neurons", "Microglia", "organoids"]
    mock_label = "mock"
    virus_label = "SARS-CoV"
    #make the df to lood nice in heatmap
    heatmap_df = count_data_df.clip(lower=64) #Replaces every value < 64 with 64
    heatmap_df = np.log2(heatmap_df)
    if samples == "one": #z-score on all samples - we can see only differeances between samples and not between virues and mock
        #row-wise z-score (better pattern visibility)
        heatmap_df = heatmap_df.sub(heatmap_df.mean(axis=1), axis=0)
        heatmap_df = heatmap_df.div(heatmap_df.std(axis=1) + 1e-6, axis=0)
    else: #z-score for each condition
        for cond in conditions:
            cond_cols = [c for c in heatmap_df.columns if get_condition_from_col(c, mock_label, virus_label) == cond]
            count_data = heatmap_df[cond_cols]
            count_data = count_data.sub(count_data.mean(axis=1), axis=0)
            count_data = count_data.div(count_data.std(axis=1) + 1e-6, axis=0)
            heatmap_df.update(count_data)
    return(heatmap_df)

def plot_heatmap(virus_label,up_genes, down_genes, count_data, cond, out_dir, cell_colors_map, treatment_colors_map, cell_types, treatments, df):
    output_heatmap = out_dir + f"heatmap_{cond}_deseq2_FC_1.5.png"
    sig_genes = list(set(up_genes) | set(down_genes))
    #sig_genes = set().union(*up_genes.values(), *down_genes.values())
    print(f"Total significant genes: {len(sig_genes)}")
    heatmap_df = count_data.loc[count_data.index.intersection(sig_genes)]
    heatmap_df = arrangeDS(heatmap_df, "one")
    output_heatmap_df = out_dir + f"df_{cond}_oneCond.csv"
    heatmap_df.to_csv(output_heatmap_df, index = True)
    #only the col from one sample
    col_colors = ["black" if virus_label in c else "pink" for c in heatmap_df.columns]
    color_lut = {"SARS-CoV": "black","mock": "pink"}
    g = sns.clustermap(heatmap_df, cmap="coolwarm", center=0, col_colors=col_colors, cbar_pos=(0.02, 0.8, 0.03, 0.18),
        col_cluster=False, xticklabels=False, yticklabels=False, figsize=(12, 14))
    handles = [mpatches.Patch(color=color_lut[name], label=name) for name in color_lut]
    #g.ax_col_dendrogram.legend(handles=handles, title=cond, loc="center", fontsize=20, bbox_to_anchor=(0.5, 0.7), ncol=2, frameon=False)
    leg = g.ax_col_dendrogram.legend(handles=handles, title=cond, loc="center", bbox_to_anchor=(0.5, 0.7), ncol=2,
        frameon=False, fontsize=20)       # legend labels
    leg.get_title().set_fontsize(20)   # legend title
    g.savefig(output_heatmap, dpi=300, bbox_inches="tight")    
    
    # Get the ordered gene list
    ordered_genes = heatmap_df.index[g.dendrogram_row.reordered_ind]
    heatmap_df_ordered = df.loc[ordered_genes]
    heatmap_df_ordered = arrangeDS(heatmap_df_ordered, "all")
    output_heatmap_df = out_dir + f"df_{cond}_allCond.csv"
    heatmap_df_ordered.to_csv(output_heatmap_df, index = True)
    #all the samples
    col_colors = pd.DataFrame({
        "Cell Type": [cell_colors_map.get(ct, "gray") for ct in cell_types],
        "Treatment": [treatment_colors_map.get(tr, "gray") for tr in treatments]
    }, index=df.columns)
    # Create custom legend
    cell_legend_1 = Line2D([0], [0], color="#1f77b4", lw=6, label="Cortical_neurons")
    cell_legend_2 = Line2D([0], [0], color="#ff7f0e", lw=6, label="Cortical_neurons_microglia")
    cell_legend_3 = Line2D([0], [0], color="#2ca02c", lw=6, label="Microglia")
    cell_legend_4 = Line2D([0], [0], color="#d62728", lw=6, label="organoids")
    
    treatment_legend_1 = Line2D([0], [0], color="pink", lw=6, label="mock", markeredgecolor="pink", marker='s', markersize=10)
    treatment_legend_2 = Line2D([0], [0], color="black", lw=6, label="SARS-CoV", marker='s', markersize=10)
    output_heatmap = out_dir + f"heatmap_{cond}_deseq2_FC_1.5_allSamples.png"
    g = sns.clustermap(heatmap_df_ordered, cmap="coolwarm", center=0, col_colors=col_colors, cbar_pos=(0.02, 0.8, 0.03, 0.18),
        col_cluster=False, row_cluster=False, xticklabels=False, yticklabels=False, figsize=(12, 14))
    
    # Add legend to the plot
    plt.legend(handles=[cell_legend_1, cell_legend_2, cell_legend_3, cell_legend_4,
                       treatment_legend_1, treatment_legend_2],
               loc='center', bbox_to_anchor=(1.05, -0.8), fontsize=12,
               title="Cell Type / Treatment")
    
    #plt.title(f"heatmap_{cond}", loc="center")
    titleName = f"heatmap_{cond}_" + str(len(heatmap_df_ordered)) + "_genes_all_samples"
    g.ax_heatmap.set_title(titleName, y=1.25, fontsize=20)
    g.savefig(output_heatmap, dpi=300, bbox_inches="tight")    
    return

def enrichmen(background_genes, sig_genes, out_dir, cond, direction):
    #this is not working well, need to convert ensembl IDs to gene names first
    print(cond, len(sig_genes), len(background_genes))
    enr = gp.enrichr(
        gene_list=list(sig_genes),
        gene_sets="GO_Biological_Process_2021",
        organism="Human",      #
        background=background_genes,
        cutoff=0.05)
    enrich_res = enr.results
    output_file = out_dir + f"GO_BP_{cond}_{direction}_deseq2_FC_1.5.csv"
    enrich_res.to_csv(output_file, index=False)
    return enrich_res

def get_condition_from_col(col_name, mock_label, virus_label):
    """Extract the condition from column name"""
    if mock_label in col_name:
        idx = col_name.index(mock_label)
    elif virus_label in col_name:
        idx = col_name.index(virus_label)
    else:
        idx = len(col_name)
    return col_name[:idx].rstrip("_")

def dict_to_table(genes_by_condition):
    #Convert a dictionary of gene sets to a boolean dataframe
    all_genes = sorted(set().union(*genes_by_condition.values()))
    df_bool = pd.DataFrame(
    {
        cond: [gene in genes_by_condition[cond] for gene in all_genes]
        for cond in genes_by_condition.keys()
    },
    index=all_genes)
    return df_bool

def plot_venn(gene_dict, title, output_file):
    cell_colors_map = {
        "Cortical_neurons": "#1f77b4",
        "Cortical_neurons_microglia": "#ff7f0e",
        "organoids": "#d62728",
        "Microglia": "#2ca02c"
    }
    colors = [cell_colors_map[name] for name in gene_dict.keys()]
    plt.figure(figsize=(8, 8))
    venn(gene_dict,
        #cmap=cell_colors_map,  # this is not working
        cmap=colors,  
        alpha=0.6)

    plt.title(title)
    plt.savefig(output_file, bbox_inches="tight")
    plt.close()
    return

def main():
    #min_samples = 2 #min 2 samples with less than min expression counts
    min_expression = 64 #min expression value
    fc_threshold = 1.5 #log2​(1)=0
    pval_threshold = 0.05
    #pval_threshold = 0.1
    # df: rows=genes, columns=samples
    # columns names contain both condition and treatment info, e.g.
    # 'Cortical_neurons_microglia_mock_01', 'Cortical_neurons_microglia_SARS-CoV_02', etc.
    in_dir = '/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/expression_table/'
    out_dir = '/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/DE/'
    input_table = in_dir + 'expression_table_HN1.txt'
    counts_table = pd.read_csv(input_table, delimiter =  "\s|\t", engine='python', skiprows=1)        
    counts_table = counts_table.reset_index().set_index('Geneid', drop=True)
    counts_table.drop(['index','Chr','Start','End','Strand','Length'], axis=1, inplace=True) 
    counts_table = sort_table(counts_table)
    counts_table.rename(columns=lambda x: x.replace('.sort.bam', ''), inplace=True)
    conditions = ["Cortical_neurons_microglia", "Cortical_neurons", "Microglia", "organoids"]
    mock_label = "mock"
    virus_label = "SARS-CoV"
    genes_details_input = in_dir + "biomart_genes_details.txt"
    genes_details = pd.read_csv(genes_details_input, sep="\t", index_col=0)  
    genes_details = genes_details.iloc[:, [0, -1]] #first and last col only
    # differential expression analysis
    up_genes, down_genes, de_results = differential_expression_deseq2(counts_table, conditions, mock_label, virus_label, fc_threshold, pval_threshold, min_expression, out_dir)
    output_file = out_dir + "Dseq_results_all_conditions_FC_1.5.csv"
    de_results = de_results.set_index("Gene")
    de_results = de_results.join(genes_details)
    de_results.to_csv(output_file, index=True)
    output_file_up = out_dir + "Dseq_results_upregulated_SARS-CoV_FC_1.5.csv"
    de_up = dict_to_table(up_genes)
    de_up = de_up.join(genes_details)
    de_up.to_csv(output_file_up, index=True)
    output_file_down = out_dir + "Dseq_results_upregulated_mock_FC_1.5.csv"
    de_down = dict_to_table(down_genes)
    de_down = de_down.join(genes_details)
    de_down.to_csv(output_file_down, index=True)
    print("DE analysis finished!")
    plot_venn(up_genes, "Genes up in SARS-CoV", out_dir + "Venn_up_SARS-CoV_dseq_FC_1.5.png")
    plot_venn(down_genes, "Genes up in mock", out_dir + "Venn_up_mock_dseq_FC_1.5.png")
    return

if __name__ == "__main__":
    main() 
