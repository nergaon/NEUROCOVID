import pandas as pd

def count_unique_genes(row):
    genes = set()
    for val in row:
        if isinstance(val, str):
            genes.update(g.strip() for g in val.split(","))
    return len(genes)

def functional_category(go_name):
    n = go_name.lower()
    if any(x in n for x in [
        "virus", "interferon", "antiviral"
    ]):
        return "Antiviral / Interferon response"
    if any(x in n for x in [
        "inflammatory", "cytokine", "immune", "lipopolysaccharide"
    ]):
        return "Inflammation & immune signaling"
    if "apopt" in n or "cell death" in n:
        return "Apoptosis & cell survival"
    if any(x in n for x in [
        "migration", "adhesion", "chemotaxis"
    ]):
        return "Cell migration & tissue remodeling"
    if any(x in n for x in [
        "extracellular matrix", "collagen"
    ]):
        return "Extracellular matrix organization"
    if any(x in n for x in [
        "synap", "neuro", "dopamine", "nervous"
    ]):
        return "Neuronal & synaptic processes"
    if any(x in n for x in [
        "ion transport", "membrane depolarization"
    ]):
        return "Ion transport & excitability"
    if any(x in n for x in [
        "gene expression", "biosynthetic", "transcription"
    ]):
        return "Gene regulation"
    return "Other"

def build_go_matrix(all_results_GO_df):
    all_results_GO_df["Genes_list"] = all_results_GO_df["Genes"].str.split(";")
    go_dict = {}
    for _, row in all_results_GO_df.iterrows():
        term = row["Term"]
        cond_dir = row["cond_dir"]
        genes = set(row["Genes_list"])
        if term not in go_dict:
            go_dict[term] = {}
        go_dict[term][cond_dir] = genes
    go_matrix = pd.DataFrame.from_dict(go_dict, orient="index")
    # convert sets to comma-separated strings
    go_matrix = go_matrix.applymap(lambda x: ",".join(sorted(x)) if isinstance(x, set) else "")
    return go_matrix

def remove_redundant_go(all_results_GO_df):
    term_to_genes = (
        all_results_GO_df
        .groupby("Term")["Genes"]
        .apply(lambda x: set(";".join(x).split(";")))
        .to_dict())
    terms = list(term_to_genes.keys())
    to_remove = set()
    for i in range(len(terms)):
        for j in range(len(terms)):
            if i == j:
                continue
            genes_i = term_to_genes[terms[i]]
            genes_j = term_to_genes[terms[j]]
            if genes_i.issubset(genes_j) and len(genes_i) < len(genes_j):
                to_remove.add(terms[i])
    filtered_df = all_results_GO_df[~all_results_GO_df["Term"].isin(to_remove)]
    filtered_df["Parent_Process"] = filtered_df["Term"].apply(get_parent_process)
    return filtered_df

def get_parent_process(term_name):
    for go_id, go_obj in go_dag.items():
        if go_obj.name == term_name:
            parents = go_obj.get_all_parents()
            parent_names = [go_dag[p].name for p in parents]
            return ";".join(parent_names[:3])  # first 3 parents
    return ""

def collapse_identical_gene_sets(df):
    """
    df must contain:
        Term
        Genes  (semicolon-separated string)
        Adjusted P-value
    """
    # Convert gene string to sorted canonical representation
    df["GeneSetKey"] = df["Genes"].apply(
        lambda x: ";".join(sorted(set(x.split(";")))))
    # For identical gene sets → keep row with smallest adjusted p-value
    collapsed = (
        df.sort_values("Adjusted P-value")
          .drop_duplicates(subset="GeneSetKey", keep="first")
          .drop(columns="GeneSetKey"))
    return collapsed

def main():
    filtered_GO_df = remove_redundant_go(all_results_GO_df)
    go_matrix = build_go_matrix(filtered_GO_df)
    padj_summary = (filtered_GO_df.groupby("Term")["Adjusted P-value"].min())
    go_matrix["Min_Adjusted_Pvalue"] = padj_summary
    output_file = out_dir + "GO/Combined_GO_matrix.xlsx"
    go_matrix.to_excel(output_file, index=True, engine="openpyxl")
    
    df = pd.read_excel("/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/DE/GO/Combined_GO_matrix.xlsx", index_col=0)
    df["Min_Adjusted_Pvalue"] = pd.to_numeric(df["Min_Adjusted_Pvalue"], errors="coerce")
    df_filtered = df[df["Min_Adjusted_Pvalue"] <= 0.05].copy()
    gene_cols = df_filtered.columns.drop("Min_Adjusted_Pvalue")
    df_filtered["Unique_gene_count"] = df_filtered[gene_cols].apply(count_unique_genes, axis=1)
    df_filtered["Functional_category"] = [functional_category(name) for name in df_filtered.index]
    df_sorted = df_filtered.sort_values(["Functional_category", "Min_Adjusted_Pvalue"],ascending=[True, True])
    with pd.ExcelWriter("/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/DE/GO/GO_functional_summary.xlsx",
        engine="openpyxl") as writer:
        row = 0
        for category, group in df_sorted.groupby("Functional_category"):
            # Title row
            pd.DataFrame([[category.upper()]]).to_excel(
                writer,
                startrow=row,
                header=False,
                index=False)
            row += 1
            group.to_excel(writer, startrow=row)
            row += len(group) + 2
    return

if __name__ == "__main__":
    main() 