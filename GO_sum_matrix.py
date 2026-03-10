import pandas as pd
import re
from glob import glob
from itertools import combinations

# ---------------------------------------------------------
# 3️⃣ Extract significant genes from cell
# ---------------------------------------------------------
def extract_sig_genes(cell, alpha=0.05):
    if pd.isna(cell):
        return []

    cell = str(cell)

    if "| p=" not in cell:
        return []

    genes_part, p_part = cell.split("| p=")

    try:
        pval = float(p_part.strip())
    except:
        return []

    if pval <= alpha:
        genes = [g.strip() for g in genes_part.split(",") if g.strip()]
        return genes
    else:
        return []

# ---------------------------------------------------------
# 4️⃣ Count unique genes across all significant conditions
# ---------------------------------------------------------
def count_unique_genes(row, condition_cols):
    genes = set()
    for col in condition_cols:
        genes.update(extract_sig_genes(row[col]))
    return pd.Series({
        "Unique_gene_count": len(genes),
        "Unique_genes": ", ".join(sorted(genes))})

# ---------------------------------------------------------
# 5️⃣ Clean GO term (remove + and GO ID)
# ---------------------------------------------------------
def clean_term(term):
    term = term.replace("+", "").strip()
    term = re.sub(r"\(GO:\d+\)", "", term)
    return term.strip()

# ---------------------------------------------------------
# 7️⃣ Functional category assignment
# ---------------------------------------------------------
def classify_term(term):
    t = term.lower()
    if any(k in t for k in ["antiviral", "interferon", "mhc", "antigen"]):
        return "Antiviral / Antigen Presentation"
    if any(k in t for k in ["inflammatory", "immune", "b cell", "cytokine"]):
        return "Inflammatory / Immune"
    if any(k in t for k in ["synaptic", "neuron", "learning"]):
        return "Neuronal / Synaptic"
    if any(k in t for k in ["actin", "junction", "cytoskeleton"]):
        return "Cytoskeleton / Adhesion"
    if any(k in t for k in ["apoptotic"]):
        return "Cell Death"
    return "Other"

def main():
    files = glob("/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/DE/GO/*.csv") 
    dfs = [pd.read_csv(f) for f in files]
    df = pd.concat(dfs, ignore_index=True)
    #print(df.columns)
    #Combine genes + p-value into one cell
    df["value"] = (df["Genes"].str.replace(";", ",", regex=False) + " | p=" + df["Adjusted P-value"].astype(str))
    #Pivot to matrix
    matrix = df.pivot_table(index="Term", columns="cond_dir", values="value", aggfunc="first")
    #Add minimum adjusted p-value
    matrix["Min_Adjusted_Pvalue"] = (df.groupby("Term")["Adjusted P-value"].min())
     # 1️⃣ Filter significant GO terms
    df_sig =  matrix[matrix["Min_Adjusted_Pvalue"] <= 0.05].copy()
    df_sig.to_excel("/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/DE/GO/Combined_GO_matrix.xlsx", engine="openpyxl")
    # 2️⃣ Identify condition columns
    # (everything except Term and Min_Adjusted_Pvalue)
    df_sig = df_sig.reset_index().rename(columns={"index": "Term"})
    condition_cols = [c for c in df_sig.columns if c not in ["Term", "Min_Adjusted_Pvalue"]]
    df_sig[["Unique_gene_count", "Unique_genes"]] = df_sig.apply(lambda r: count_unique_genes(r, condition_cols), axis=1)
    #print(df_sig.columns)
    df_sig["Clean_term"] = df_sig["Term"].apply(clean_term)

    # ---------------------------------------------------------
    # 6️⃣ Remove redundant GO terms based on gene sets
    #    - drop any term whose gene set is identical to or a strict
    #      subset of another term's gene set
    #    - keep the term with the largest gene set (and lowest
    #      p-value when counts tie)
    # ---------------------------------------------------------
    def _filter_by_gene_sets(df_input):
        df_work = df_input.copy()
        # convert the comma-separated gene list into a set
        df_work["Gene_set"] = df_work["Unique_genes"].apply(
            lambda s: set(g.strip() for g in s.split(",") if g.strip())
        )
        # sort so that largest sets (and best p-values) come first
        df_work = df_work.sort_values(
            ["Unique_gene_count", "Min_Adjusted_Pvalue"],
            ascending=[False, True],
        )
        keep_idx = []
        kept_sets = []
        for idx, row in df_work.iterrows():
            gs = row["Gene_set"]
            # check if this set is subset/equal of any already kept set
            if any(gs <= ks for ks in kept_sets):
                continue
            keep_idx.append(idx)
            kept_sets.append(gs)
        return df_work.loc[keep_idx].drop(columns=["Gene_set"])

    df_sig = _filter_by_gene_sets(df_sig)

    # 7️⃣ Collapse redundant GO terms (after filtering by gene sets)
    # (keep lowest p-value entry per clean term)
    collapsed = (df_sig.sort_values("Min_Adjusted_Pvalue").groupby("Clean_term", as_index=False).first())
    collapsed["Functional_Category"] = collapsed["Clean_term"].apply(classify_term)
    # 8️⃣ Final table (keep ALL original columns)
    final_cols = (["Term", "Functional_Category"] + list(matrix.columns) + ["Unique_gene_count", "Unique_genes"])
    final_table = collapsed[final_cols].sort_values(["Functional_Category", "Min_Adjusted_Pvalue"])
    # 9️⃣ Save to Excel
    final_table.to_excel("/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/DE/GO/GO_functional_summary.xlsx", index=False)
    return

if __name__ == "__main__":
    main() 