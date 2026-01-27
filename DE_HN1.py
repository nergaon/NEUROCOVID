import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

def differential_expression(df, conditions, mock_label, virus_label, fc_threshold, pval_threshold, min_expression, min_samples):
    """
    Perform DE analysis for multiple conditions.
    Returns:
        up_genes: dict of genes up in virus
        down_genes: dict of genes up in mock
        results_df: dataframe of all DE results
    """
    all_results = []
    up_genes_all = set()
    down_genes_all = set()
    for cond in conditions:
        # select relevant columns. 
        cond_cols = [c for c in df.columns if get_condition_from_col(c, mock_label, virus_label) == cond]
        if len(cond_cols) == 0:
            print(f"No columns found for {cond}, skipping")
            continue
        print(cond, cond_cols)
        mock_cols = [c for c in cond_cols if mock_label in c]
        virus_cols = [c for c in cond_cols if virus_label in c]
        if len(mock_cols) == 0 or len(virus_cols) == 0:
            print(f"No mock or virus columns for {cond}, skipping")
            continue
        mock_data = df[mock_cols]
        virus_data = df[virus_cols]
        # filter genes with low expression in the condition
        mask = ((mock_data >= min_expression).sum(axis=1) >= min_samples) | ((virus_data >= min_expression).sum(axis=1) >= min_samples)
        mock_data = mock_data[mask]
        virus_data = virus_data[mask]
        # iterate over genes
        for gene in df.index:
            mock_vals = mock_data.loc[gene].values
            virus_vals = virus_data.loc[gene].values
             # skip if all values identical
            if np.all(mock_vals == mock_vals[0]) and np.all(virus_vals == virus_vals[0]):
                continue
            # t-test
            t_stat, p_val = ttest_ind(virus_vals, mock_vals, equal_var=False)
            # fold change
            mean_mock = np.mean(mock_vals) + 1e-9 #1e-9 is a small “pseudo-count” added to avoid division by zero
            mean_virus = np.mean(virus_vals) + 1e-9
            fc = mean_virus / mean_mock
            log2fc = np.log2(fc)
            all_results.append({
                "Gene": gene,
                "Condition": cond,
                "Mean_Mock": mean_mock,
                "Mean_Virus": mean_virus,
                "Log2FC": log2fc,
                "PValue": p_val})
    # build dataframe
    results_df = pd.DataFrame(all_results)
    # FDR correction across all genes and conditions
    results_df["FDR"] = multipletests(results_df["PValue"], method="fdr_bh")[1]
    # select significant genes based on thresholds
    #sig_df = results_df[results_df["FDR"] < pval_threshold]
    sig_df = results_df[results_df["PValue"] < pval_threshold]
    up_genes_all = set(sig_df[sig_df["Log2FC"] > np.log2(fc_threshold)]["Gene"])
    down_genes_all = set(sig_df[sig_df["Log2FC"] < -np.log2(fc_threshold)]["Gene"])
    return up_genes_all, down_genes_all, results_df

def get_condition_from_col(col_name, mock_label, virus_label):
    """Extract the condition from column name, assumes <Condition>_<Treatment>_<Number>"""
    if mock_label in col_name:
        idx = col_name.index(mock_label)
    elif virus_label in col_name:
        idx = col_name.index(virus_label)
    else:
        idx = len(col_name)
    return col_name[:idx].rstrip("_")

def main():
    min_samples = 2 #min 2 samples with less than min expression counts
    min_expression = 6 #min log2 expression value
    fc_threshold=1.1 
    pval_threshold=0.05
    # df: rows=genes, columns=samples
    # columns names contain both condition and treatment info, e.g.
    # 'Cortical_neurons_microglia_mock_01', 'Cortical_neurons_microglia_SARS-CoV_02', etc.
    in_dir = '/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/expression_table/'
    out_dir = '/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/DE/'
    input_table = in_dir + 'log_express_gene.txt'
    df = pd.read_csv(input_table, sep="\t", index_col=0).iloc[:, :-3]
    conditions = ["Cortical_neurons_microglia", "Cortical_neurons", "Microglia", "organoids"]
    mock_label = "mock"
    virus_label = "SARS-CoV"
    # differential expression analysis
    up_genes, down_genes, de_results = differential_expression(df, conditions, mock_label, virus_label, fc_threshold, pval_threshold, min_expression, min_samples)
    print("up genes:", up_genes)
    print("down genes:", down_genes)
    # -------------------------------
    # SAVE RESULTS
    # -------------------------------
    output_file = out_dir + "DE_results_all_conditions.csv"
    de_results.to_csv(output_file, index=False)
    output_file_up = out_dir + "DE_results_upregulated_SARS-CoV.csv"
    de_results[(de_results["Gene"].isin(up_genes))].to_csv(output_file_up, index=False)
    output_file_down = out_dir + "DE_results_upregulated_mock.csv"
    de_results[(de_results["Gene"].isin(down_genes))].to_csv(output_file_down, index=False)
    print("DE analysis finished!")
    print(f"Genes up in SARS-CoV: {len(up_genes)}")
    print(f"Genes up in mock: {len(down_genes)}")

    return

if __name__ == "__main__":
    main() 
