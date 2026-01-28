import pandas as pd
import numpy as np
import warnings
from scipy.stats import ttest_ind
from scipy.stats._stats_py import _ttest_ind_from_stats
from scipy.stats._warnings_errors import DegenerateDataWarning
from statsmodels.stats.multitest import multipletests
from venn import venn
import matplotlib.pyplot as plt

def differential_expression(df, conditions, mock_label, virus_label, fc_threshold, pval_threshold, min_expression, min_samples):
    """
    Perform DE analysis for multiple conditions.
    Returns:
        up_genes: dict of genes up in virus
        down_genes: dict of genes up in mock
        results_df: dataframe of all DE results
    """
    min_sd = 1e-4   # for log-scaled data. fileter genes with very low variance
    # use 1e-3 if data are raw counts or CPM/TPM
    all_results = []
    up_genes = {}      # genes up in virus per condition
    down_genes = {}    # genes up in mock per condition
    problem_genes = set() # genes causing precision loss warnings
    for cond in conditions:
        cond_cols = [c for c in df.columns if get_condition_from_col(c, mock_label, virus_label) == cond]
        if len(cond_cols) == 0:
            print(f"No columns for {cond}, skipping")
            continue
        mock_cols = [c for c in cond_cols if mock_label in c]
        virus_cols = [c for c in cond_cols if virus_label in c]
        if len(mock_cols) == 0 or len(virus_cols) == 0:
            print(f"No mock/virus columns for {cond}, skipping")
            continue
        mock_data = df[mock_cols]
        virus_data = df[virus_cols]
        # expression filter. at least 2 samples with min expression, from either condition
        expr_mask = ((mock_data >= min_expression).sum(axis=1) >= min_samples) | \
               ((virus_data >= min_expression).sum(axis=1) >= min_samples)
        mock_mean  = mock_data.mean(axis=1)
        virus_mean = virus_data.mean(axis=1)
        between_sd = pd.concat([mock_mean, virus_mean], axis=1).std(axis=1)
        sd_mask = between_sd >= min_sd
        mask = expr_mask & sd_mask
        mock_data = mock_data[mask]
        virus_data = virus_data[mask]
        cond_results = []
        for gene in mock_data.index:   
            mock_vals = mock_data.loc[gene].values
            virus_vals = virus_data.loc[gene].values
            if np.all(mock_vals == mock_vals[0]) and np.all(virus_vals == virus_vals[0]):
                continue
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                t_stat, p_val = ttest_ind(virus_vals, mock_vals, equal_var=False)
                if any("Precision loss occurred" in str(wi.message) for wi in w):
                    problem_genes.add(gene)
            #t_stat, p_val = ttest_ind(virus_vals, mock_vals, equal_var=False)
            mean_mock = np.mean(mock_vals) + 1e-9 # to avoid division by zero
            mean_virus = np.mean(virus_vals) + 1e-9
            log2fc = np.log2(mean_virus / mean_mock)
            cond_results.append([gene, cond, mean_mock, mean_virus, log2fc, p_val])
        cond_df = pd.DataFrame(cond_results, columns=["Gene","Condition","Mean_Mock", "Mean_Virus","Log2FC","PValue"])
        cond_df["FDR"] = multipletests(cond_df["PValue"], method="fdr_bh")[1]
        # Significant genes
        #sig = cond_df[(cond_df["FDR"] < pval_threshold) & (abs(cond_df["Log2FC"]) > np.log2(fc_threshold))]
        sig = cond_df[(cond_df["PValue"] < pval_threshold) & (abs(cond_df["Log2FC"]) > np.log2(fc_threshold))]
        up_genes[cond] = set(sig[sig["Log2FC"] > np.log2(fc_threshold)]["Gene"])
        down_genes[cond] = set(sig[sig["Log2FC"] < -np.log2(fc_threshold)]["Gene"])
        all_results.append(cond_df)
    results_df = pd.concat(all_results, ignore_index=True)
    print(f"{len(problem_genes)} genes cause precision loss")
    print(problem_genes) 
    problem = pd.DataFrame({"gene": list(problem_genes)})
    problem.to_csv("/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/DE/problem_genes.csv", index=False)
    return up_genes, down_genes, results_df

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
    min_samples = 2 #min 2 samples with less than min expression counts
    min_expression = 6 #min log2 expression value
    fc_threshold = 1.1 #log2​(1.1)≈0.138
    pval_threshold = 0.05
    # df: rows=genes, columns=samples
    # columns names contain both condition and treatment info, e.g.
    # 'Cortical_neurons_microglia_mock_01', 'Cortical_neurons_microglia_SARS-CoV_02', etc.
    in_dir = '/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/expression_table/'
    out_dir = '/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/DE/'
    input_table = in_dir + 'log_express_gene.txt'
    #df = pd.read_csv(input_table, sep="\t", index_col=0).iloc[:, :-3]
    df = pd.read_csv(input_table, sep="\t", index_col=0)
    df_num = df.iloc[:, :-3]  # remove last 3 columns with summary stats
    genes_details = df.iloc[:, -1:]  # store the last columns which is gene name
    conditions = ["Cortical_neurons_microglia", "Cortical_neurons", "Microglia", "organoids"]
    mock_label = "mock"
    virus_label = "SARS-CoV"
    # differential expression analysis
    up_genes, down_genes, de_results = differential_expression(df_num, conditions, mock_label, virus_label, fc_threshold, pval_threshold, min_expression, min_samples)
    # -------------------------------
    # SAVE RESULTS
    # -------------------------------
    output_file = out_dir + "DE_results_all_conditions.csv"
    de_results = de_results.set_index("Gene")
    de_results = de_results.join(genes_details)
    de_results.to_csv(output_file, index=True)
    output_file_up = out_dir + "DE_results_upregulated_SARS-CoV.csv"
    all_up_genes = [g for genes in up_genes.values() for g in genes]
    de_up = de_results[(de_results.index.isin(all_up_genes))]
    de_up.to_csv(output_file_up, index=True)
    output_file_down = out_dir + "DE_results_upregulated_mock.csv"
    all_down_genes = [g for genes in down_genes.values() for g in genes]
    de_down = de_results[(de_results.index.isin(all_down_genes))]
    de_down.to_csv(output_file_down, index=True)
    print("DE analysis finished!")
    overlap_up = build_overlap_df(up_genes, out_dir + "DE_overlap_up_SARS-CoV.csv")
    overlap_down = build_overlap_df(down_genes, out_dir + "DE_overlap_down_mock.csv")
    plot_venn(up_genes, "Genes up in SARS-CoV", out_dir + "Venn_up_SARS-CoV.png")
    plot_venn(down_genes, "Genes up in mock", out_dir + "Venn_up_mock.png")
    return

if __name__ == "__main__":
    main() 
