import pandas as pd
import re

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

# ---------------------------------------------------------
# 1️⃣ Filter significant GO terms
# ---------------------------------------------------------
df_sig = matrix[matrix["Min_Adjusted_Pvalue"] <= 0.05].copy()

# ---------------------------------------------------------
# 2️⃣ Identify condition columns
# (everything except Term and Min_Adjusted_Pvalue)
# ---------------------------------------------------------
condition_cols = [c for c in df_sig.columns if c not in ["Term", "Min_Adjusted_Pvalue"]]
df_sig[["Unique_gene_count", "Unique_genes"]] = df_sig.apply(lambda r: count_unique_genes(r, condition_cols), axis=1)
df_sig["Clean_term"] = df_sig["Term"].apply(clean_term)
# ---------------------------------------------------------
# 6️⃣ Collapse redundant GO terms
# (keep lowest p-value entry per clean term)
# ---------------------------------------------------------
collapsed = (df_sig.sort_values("Min_Adjusted_Pvalue").groupby("Clean_term", as_index=False).first())
collapsed["Functional_Category"] = collapsed["Clean_term"].apply(classify_term)
# ---------------------------------------------------------
# 8️⃣ Final table (keep ALL original columns)
# ---------------------------------------------------------
final_cols = (
    ["Functional_Category"] +
    list(matrix.columns) +
    ["Unique_gene_count", "Unique_genes"])
final_table = collapsed[final_cols].sort_values(
    ["Functional_Category", "Min_Adjusted_Pvalue"])
# ---------------------------------------------------------
# 9️⃣ Save to Excel
# ---------------------------------------------------------
final_table.to_excel("GO_functional_summary.xlsx", index=False)