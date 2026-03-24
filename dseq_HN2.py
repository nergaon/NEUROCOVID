import pandas as pd
import numpy as np
import sys
from pathlib import Path
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from venn import venn
import matplotlib.pyplot as plt
from goatools.obo_parser import GODag
import seaborn as sns
import gseapy as gp
from gprofiler import GProfiler
from matplotlib.lines import Line2D
from plot_config import colors
from plot_config import sort_table
go_dag = GODag("/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/DE/GO/go-basic.obo")

def differential_expression_deseq2(df, conditions, mock_label, virus_label, fc_threshold, pval_threshold, min_expression, out_dir=None):
    """
    Differential expression using DESeq2 (via pydeseq2)
    Returns:
        up_genes: dict of genes up in virus per condition
        down_genes: dict of genes up in mock per condition
        results_df: dataframe of all DE results
    """
    all_results_GO = []
    all_results = []
    up_genes = {}
    down_genes = {}
    print("size of df:", df.shape)
    rows = []
    rows.append({"condition": "All", "n_genes": df.shape[0], "SARS-CoV-2_genes": np.nan, "mock_genes": np.nan})
    cell_colors_map, treatment_colors_map, treatment_markers, cell_types, treatments = colors(df.columns)
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
        min_samples = min(len(mock_cols), len(virus_cols))  #gene is express in at least min samples
        if len(mock_cols) == 0 or len(virus_cols) == 0:
            print(f"No mock/virus columns for {cond}, skipping")
            continue
        count_data = df[cond_cols]
        # ---- Expression filter (DESeq2-style prefiltering) ----
        mask = (count_data > min_expression).sum(axis=1) >= min_samples
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
        sig = res[(res["FDR"] < pval_threshold) & (abs(res["Log2FC"]) > np.log2(fc_threshold))]
        up_genes[cond] = set(sig[sig["Log2FC"] > 0]["Gene"])
        down_genes[cond] = set(sig[sig["Log2FC"] < 0]["Gene"])
        rows.append({"condition": cond, "n_genes": count_data.shape[0], 
                     "SARS-CoV-2_genes": len(up_genes[cond]), "mock_genes": len(down_genes[cond])})
        res_GO = enrichment(background_genes, up_genes[cond], out_dir, cond, "SARS-CoV")
        all_results_GO.append(res_GO)
        res_GO = enrichment(background_genes, down_genes[cond], out_dir, cond, "mock")
        all_results_GO.append(res_GO)
        # plot heatmaps for up and down genes separately, ordered by DE stat
        plot_heatmap(
            virus_label,
            up_genes[cond],
            down_genes[cond],
            count_data,
            cond,
            out_dir,
            cell_colors_map,
            treatment_colors_map,
            cell_types,
            treatments,
            df,
            res
        )
    #de analysis for all the mock cols vs all SARS-CoV cols
    count_data_all = df.copy()
    # ---- Expression filter (DESeq2-style prefiltering) ----
    mask = (count_data_all > min_expression).sum(axis=1) >= 7 
    count_data_all = count_data_all.loc[mask]
    print(f"All_conditions: {count_data_all.shape[0]} genes pass expression filter")
    coldata_all = pd.DataFrame(index=count_data_all.columns)
    #Build metadata with BOTH factors
    coldata_all["treatment"] = ["virus" if virus_label in c else "mock" for c in count_data_all.columns]
    coldata_all["condition"] = [get_condition_from_col(c, mock_label, virus_label) for c in count_data_all.columns]
    #Run DESeq2 with covariate model
    dds_all = DeseqDataSet(
        counts=count_data_all.T,
        metadata=coldata_all,
        design="~ condition + treatment",
        refit_cooks=True)
    dds_all.deseq2()
    #Test virus vs mock effect
    stats_all = DeseqStats(dds_all, contrast=("treatment", "virus", "mock"))
    stats_all.summary()
    res_all = stats_all.results_df.copy()
    res_all["Gene"] = res_all.index
    res_all["Condition"] = "All_conditions"
    res_all = res_all.rename(columns={
        "log2FoldChange": "Log2FC",
        "padj": "FDR",
        "pvalue": "PValue"})
    all_results.append(res_all)
    #Extract significant genes
    sig_all = res_all[
        (res_all["FDR"] < pval_threshold) &
        (abs(res_all["Log2FC"]) > np.log2(fc_threshold))]
    up_genes["All_conditions"] = set(sig_all[sig_all["Log2FC"] > 0]["Gene"])
    down_genes["All_conditions"] = set(sig_all[sig_all["Log2FC"] < 0]["Gene"])
    rows.append({"condition": "All_conditions", "n_genes": count_data_all.shape[0], 
                     "SARS-CoV-2_genes": len(up_genes["All_conditions"]), "mock_genes": len(down_genes["All_conditions"])})
    enrichment(background_genes, up_genes["All_conditions"], out_dir, "all", "SARS-CoV")
    enrichment(background_genes, down_genes["All_conditions"], out_dir, "all", "mock")
    plot_heatmap(
        virus_label,
        up_genes["All_conditions"],
        down_genes["All_conditions"],
        df,
        "all",
        out_dir,
        cell_colors_map,
        treatment_colors_map,
        cell_types,
        treatments,
        df,
        res_all
    )
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

def get_ordered_genes_by_stat(sig_genes, de_results, descending=True):
    if len(sig_genes) == 0:
        return []
    stat_col = "stat" if "stat" in de_results.columns else "statistic"
    if stat_col not in de_results.columns:
        return sorted(sig_genes)
    ordered = (
        de_results[de_results["Gene"].isin(sig_genes)]
        .dropna(subset=[stat_col])
        .sort_values(stat_col, ascending=not descending)
    )
    # Preserve stat order but keep each gene name once.
    return list(dict.fromkeys(ordered["Gene"].tolist()))


def order_frame_by_gene_list(frame, ordered_genes):
    """Order rows by gene list without using reindex (safe with duplicate index labels)."""
    if frame.empty or len(ordered_genes) == 0:
        return frame
    rank_map = {gene: i for i, gene in enumerate(ordered_genes)}
    ordered_frame = frame.loc[frame.index.isin(rank_map)].copy()
    if ordered_frame.empty:
        return ordered_frame
    ordered_frame["_gene_rank"] = ordered_frame.index.map(rank_map)
    ordered_frame = ordered_frame.sort_values("_gene_rank", kind="stable")
    ordered_frame = ordered_frame.drop(columns=["_gene_rank"])
    return ordered_frame


def plot_heatmap_direction(
    virus_label,
    direction,
    count_data,
    cond,
    out_dir,
    cell_colors_map,
    treatment_colors_map,
    cell_types,
    treatments,
    df,
    ordered_genes
):
    if len(ordered_genes) == 0:
        print(f"No {direction} genes for {cond}, skipping heatmap")
        return

    heatmap_df = order_frame_by_gene_list(count_data, ordered_genes).dropna(how="all")
    heatmap_df = arrangeDS(heatmap_df, "one")

    output_heatmap_df = out_dir + f"df_{cond}_{direction}_oneCond.csv"
    heatmap_df.to_csv(output_heatmap_df, index=True)

    # Keep sample annotations as top bars (Cell Type + Treatment); keep numeric colorbar on the side.
    col_colors = pd.DataFrame({
        "Cell Type": [
            cell_colors_map.get(get_condition_from_col(c, "mock", virus_label), "gray")
            for c in heatmap_df.columns
        ],
        "Treatment": [
            treatment_colors_map.get("SARS-CoV" if virus_label in c else "mock", "gray")
            for c in heatmap_df.columns
        ]
    }, index=heatmap_df.columns)
    g = sns.clustermap(
        heatmap_df,
        cmap="coolwarm",
        center=0,
        col_colors=col_colors,
        cbar_pos=(0.01, 0.22, 0.012, 0.22),
        row_cluster=False,
        col_cluster=False,
        xticklabels=False,
        yticklabels=True,
        figsize=(12, 14)
    )
    g.fig.suptitle(f"heatmap_{cond}_{direction}_{len(ordered_genes)}_genes", y=0.92, fontsize=20)
    output_heatmap = out_dir + f"heatmap_{cond}_{direction}_deseq2_FC_1.5.png"
    g.savefig(output_heatmap, dpi=300, bbox_inches="tight")
    plt.close(g.fig)

    heatmap_df_ordered = order_frame_by_gene_list(df, ordered_genes).dropna(how="all")
    heatmap_df_ordered = arrangeDS(heatmap_df_ordered, "all")
    output_heatmap_df = out_dir + f"df_{cond}_{direction}_allCond.csv"
    heatmap_df_ordered.to_csv(output_heatmap_df, index=True)

    col_colors_all = pd.DataFrame({
        "Cell Type": [cell_colors_map.get(ct, "gray") for ct in cell_types],
        "Treatment": [treatment_colors_map.get(tr, "gray") for tr in treatments]
    }, index=df.columns).reindex(heatmap_df_ordered.columns)

    cell_legend_1 = Line2D([0], [0], color="#1f77b4", lw=6, label="Cortical_neurons")
    cell_legend_2 = Line2D([0], [0], color="#ff7f0e", lw=6, label="Cortical_neurons_microglia")
    cell_legend_3 = Line2D([0], [0], color="#2ca02c", lw=6, label="Microglia")
    cell_legend_4 = Line2D([0], [0], color="#d62728", lw=6, label="organoids")
    treatment_legend_1 = Line2D([0], [0], color="pink", lw=6, label="mock", markeredgecolor="pink", marker="s", markersize=10)
    treatment_legend_2 = Line2D([0], [0], color="black", lw=6, label="SARS-CoV", marker="s", markersize=10)

    g = sns.clustermap(
        heatmap_df_ordered,
        cmap="coolwarm",
        center=0,
        col_colors=col_colors_all,
        cbar_pos=(0.01, 0.22, 0.012, 0.22),
        row_cluster=False,
        col_cluster=False,
        xticklabels=False,
        yticklabels=True,
        figsize=(12, 14)
    )

    g.ax_heatmap.legend(
        handles=[
            cell_legend_1,
            cell_legend_2,
            cell_legend_3,
            cell_legend_4,
            treatment_legend_1,
            treatment_legend_2,
        ],
        loc="upper left",
        bbox_to_anchor=(-0.42, 1.0),
        fontsize=10,
        title="Cell Type / Treatment"
    )

    title_name = f"heatmap_{cond}_{direction}_{len(heatmap_df_ordered)}_genes_all_samples"
    g.fig.suptitle(title_name, y=0.97, fontsize=20)
    output_heatmap = out_dir + f"heatmap_{cond}_{direction}_deseq2_FC_1.5_allSamples.png"
    g.savefig(output_heatmap, dpi=300, bbox_inches="tight")
    plt.close(g.fig)


def plot_heatmap(virus_label, up_genes, down_genes, count_data, cond, out_dir, cell_colors_map, treatment_colors_map, cell_types, treatments, df, de_results):
    ordered_up_genes = get_ordered_genes_by_stat(up_genes, de_results, descending=True)
    ordered_down_genes = get_ordered_genes_by_stat(down_genes, de_results, descending=False)

    print(f"Total significant up genes for {cond}: {len(ordered_up_genes)}")
    print(f"Total significant down genes for {cond}: {len(ordered_down_genes)}")

    plot_heatmap_direction(
        virus_label,
        "up",
        count_data,
        cond,
        out_dir,
        cell_colors_map,
        treatment_colors_map,
        cell_types,
        treatments,
        df,
        ordered_up_genes
    )

    plot_heatmap_direction(
        virus_label,
        "down",
        count_data,
        cond,
        out_dir,
        cell_colors_map,
        treatment_colors_map,
        cell_types,
        treatments,
        df,
        ordered_down_genes
    )
    return

def enrichment(background_genes, sig_genes, out_dir, cond, direction):
    print(cond, len(sig_genes), len(background_genes))
    enr = gp.enrichr(
        gene_list=list(sig_genes),
        gene_sets="GO_Biological_Process_2025",
        organism="Human",      #
        background=background_genes,
        cutoff=0.05)
    enrich_res = enr.results.copy()
    enrich_res["cond"] = cond
    enrich_res["direction"] = direction
    enrich_res["cond_dir"] = cond + "_" + direction
    output_file = out_dir + f"GO/GO_BP_{cond}_{direction}_deseq2_FC_1.5.csv"
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


def filter_overlapping_sets(gene_dict):
    filtered = {}
    for key, genes in gene_dict.items():
        other_genes = set().union(*(other for other_key, other in gene_dict.items() if other_key != key))
        if len(set(genes) & other_genes) > 0:
            filtered[key] = set(genes)
    return filtered

def plot_venn(gene_dict, title, output_file):
    gene_dict = {key: set(value) for key, value in gene_dict.items()}
    filtered_gene_dict = filter_overlapping_sets(gene_dict)
    if len(filtered_gene_dict) >= 2:
        gene_dict = filtered_gene_dict

    cell_colors_map = {
        "Cortical_neurons": "#1f77b4",
        "Cortical_neurons_microglia": "#ff7f0e",
        "organoids": "#d62728",
        "Microglia": "#2ca02c"
    }

    labels = list(gene_dict.keys())
    colors = [cell_colors_map[name] for name in labels]
    plt.figure(figsize=(8, 8))

    if 2 <= len(gene_dict) <= 3:
        vendor_path = Path(__file__).resolve().parent / "_vendor"
        if str(vendor_path) not in sys.path:
            sys.path.append(str(vendor_path))
        from matplotlib_venn import venn2, venn3

        set_values = [gene_dict[label] for label in labels]
        if len(gene_dict) == 2:
            venn2(subsets=set_values, set_labels=labels, set_colors=colors, alpha=0.6)
        else:
            venn3(subsets=set_values, set_labels=labels, set_colors=colors, alpha=0.6)
    else:
        venn(
            gene_dict,
            cmap=colors,
            alpha=0.6
        )

    plt.title(title)
    plt.savefig(output_file, bbox_inches="tight")
    plt.close()
    return

def main():
    #min_samples = 2 #min 2 samples with less than min expression counts
    min_expression = 64 #min expression value
    fc_threshold = 1.5 #log2​(1)=0
    pval_threshold = 0.05
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
    #replace the ensembl index with gene name
    new_index = (genes_details["Gene name"].reindex(counts_table.index).fillna(pd.Series(counts_table.index, index=counts_table.index)))
    counts_table.index = new_index
    # differential expression analysis
    up_genes, down_genes, de_results = differential_expression_deseq2(counts_table, conditions, mock_label, virus_label, fc_threshold, pval_threshold, min_expression, out_dir)
    output_file = out_dir + "/Dseq_results_all_conditions_FC_1.5.csv"
    de_results = de_results.set_index("Gene")
    de_results.to_csv(output_file, index=True)
    output_file_up = out_dir + "/Dseq_results_upregulated_SARS-CoV_FC_1.5.csv"
    de_up = dict_to_table(up_genes)
    de_up = de_up.merge(genes_details, left_index=True, right_on="Gene name", how="left")
    de_up.to_csv(output_file_up, index=True)
    output_file_down = out_dir + "/Dseq_results_upregulated_mock_FC_1.5.csv"
    de_down = dict_to_table(down_genes)
    de_down = de_down.merge(genes_details, left_index=True, right_on="Gene name", how="left")
    de_down.to_csv(output_file_down, index=True)
    print("DE analysis finished!")
    up_genes.pop("All_conditions", None)
    down_genes.pop("All_conditions", None)
    plot_venn(up_genes, "Genes up in SARS-CoV", out_dir + "Venn_up_SARS-CoV_dseq_FC_1.5.png")
    plot_venn(down_genes, "Genes up in mock", out_dir + "Venn_up_mock_dseq_FC_1.5.png")
    return

if __name__ == "__main__":
    main() 
