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
        plot_heatmap(virus_label, up_genes, down_genes, count_data, cond, out_dir)
        
    results_df = pd.concat(all_results, ignore_index=True)
    sum_table = pd.DataFrame(rows)
    print(sum_table)
    return up_genes, down_genes, results_df

def plot_heatmap(virus_label,up_genes, down_genes, count_data, cond, out_dir):
    output_heatmap = out_dir + f"heatmap_{cond}_deseq2_FC_1.5.png"
    sig_genes = set().union(*up_genes.values(), *down_genes.values())
    print(f"Total significant genes: {len(sig_genes)}")
    heatmap_df = count_data.loc[count_data.index.intersection(sig_genes)]
    heatmap_df = np.log2(heatmap_df + 1)
    #row-wise z-score (better pattern visibility)
    heatmap_df = heatmap_df.sub(heatmap_df.mean(axis=1), axis=0)
    heatmap_df = heatmap_df.div(heatmap_df.std(axis=1) + 1e-6, axis=0)
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

def build_overlap_df(gene_dict, output_file):
    # Build a dataframe showing overlaps between conditions
    conds = list(gene_dict.keys())
    df_overlap = pd.DataFrame(index=conds, columns=conds, dtype=int)
    for c1 in conds:
        for c2 in conds:
            df_overlap.loc[c1, c2] = len(gene_dict[c1].intersection(gene_dict[c2]))
    df_overlap.to_csv(output_file)
    return df_overlap

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
    #cell_colors_map = {
    #    "Cortical_neurons": "#1f77b4",
    #    "Cortical_neurons_microglia": "#ff7f0e",
    #    "organoids": "#d62728",
    #    "Microglia": "#2ca02c"
    #}
    plt.figure(figsize=(8, 8))
    venn(gene_dict,
        #cmap=cell_colors_map,  # this is not working
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
    counts_table.rename(columns=lambda x: x.replace('.sort.bam', ''), inplace=True)
    conditions = ["Cortical_neurons_microglia", "Cortical_neurons", "Microglia", "organoids"]
    mock_label = "mock"
    virus_label = "SARS-CoV"
    genes_details_input = in_dir + "biomart_genes_details.txt"
    genes_details = pd.read_csv(genes_details_input, sep="\t", index_col=0)  
    genes_details = genes_details.iloc[:, [0, -1]] #first and last col only
    # differential expression analysis
    up_genes, down_genes, de_results = differential_expression_deseq2(counts_table, conditions, mock_label, virus_label, fc_threshold, pval_threshold, min_expression, out_dir)
    # SAVE RESULTS
    #my_dict = {
    #    "up_genes": up_genes,
    #    "down_genes": down_genes,
    #}
    #output_json = out_dir + "Dseq_DE_genes_dict_FC2.json"
    #with open(output_json, "w") as f:
   #     json.dump(my_dict, f, indent=2)
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
    #overlap_up = build_overlap_df(de_up.index.tolist(), out_dir + "Dseq_overlap_up_SARS-CoV.csv")
    #overlap_down = build_overlap_df(de_down.index.tolist(), out_dir + "Dseq_overlap_down_mock.csv")
    #plot_venn(up_genes, "Genes up in SARS-CoV", out_dir + "Venn_up_SARS-CoV_dseq_FC_1.5.png")
    #plot_venn(down_genes, "Genes up in mock", out_dir + "Venn_up_mock_dseq_FC_1.5.png")
    return

if __name__ == "__main__":
    main() 
