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

def main():
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