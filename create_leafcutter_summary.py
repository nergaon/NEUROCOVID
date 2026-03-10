#!/usr/bin/env python3
"""
Create a summary table from leafcutter results across different tissue types.
Keeps ALL Success clusters from ANY tissue folder.
"""
import pandas as pd
from pathlib import Path

# Define base directory
base_dir = Path("/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/leafcutter_0.2.9")
counts_file = base_dir / "NEUROCOVID_perind_numers.counts"

# Define tissue folders with results
tissue_folders = {
    "Cortical_neurons": "Cortical_neurons",
    "Cortical_neurons_microglia": "Cortical_neurons_microglia",
    "Microglia": "Microglia",
    "organoids": "organoids"
}

# Load counts file once and build cluster-to-intron mapping
print("Loading counts file and building cluster index...")
# File is space-separated, not tab-separated
counts_df = pd.read_csv(counts_file, sep=" ", index_col=0)

# Build mapping of cluster to introns in counts file
cluster_to_introns = {}
for intron_id in counts_df.index:
    # Extract cluster from intron (format: chr:start:end:cluster_id)
    parts = intron_id.split(':')
    if len(parts) >= 4:
        cluster_id = parts[-1]  # Last part is the cluster ID like "clu_1_-"
        
        if cluster_id not in cluster_to_introns:
            cluster_to_introns[cluster_id] = []
        cluster_to_introns[cluster_id].append(intron_id)

print(f"Found {len(cluster_to_introns)} unique clusters in counts file")

# Get all unique sample names
samples = list(counts_df.columns)

# Store results per cluster per tissue
cluster_data = {}  # cluster -> tissue -> {'p.adjust': val, 'max_abs_deltapsi': val, 'avg_expr_mock': val, 'avg_expr_sarscov': val, 'genes': val}

# Process each tissue type
for tissue_name, folder_name in tissue_folders.items():
    print(f"\nProcessing {tissue_name}...")
    
    tissue_dir = base_dir / folder_name
    significance_file = tissue_dir / "leafcutter_ds_cluster_significance.txt"
    effect_sizes_file = tissue_dir / "leafcutter_ds_effect_sizes.txt"
    
    # Load significance data
    sig_df = pd.read_csv(significance_file, sep="\t")
    
    # Load effect sizes data
    effect_df = pd.read_csv(effect_sizes_file, sep="\t")
    
    # Extract cluster ID from intron column in effect_df
    # Format: chr:left:right:cluster_id
    effect_df['cluster'] = effect_df['intron'].apply(lambda x: x.rsplit(':', 1)[-1])
    
    # Group by cluster and get max abs deltapsi for each
    cluster_max_deltapsi = effect_df.groupby('cluster')['deltapsi'].apply(lambda x: x.abs().max()).to_dict()
    
    # Filter to Success clusters only
    success_clusters = sig_df[sig_df['status'] == 'Success'].copy()
    
    # Extract clean cluster ID from significance file (remove chr prefix)
    success_clusters['cluster_clean'] = success_clusters['cluster'].apply(lambda x: x.split(':', 1)[1] if ':' in x else x)
    
    print(f"  Found {len(success_clusters)} Success clusters")
    
    # Identify tissue samples - precise matching to avoid overlaps
    if tissue_name == "Cortical_neurons":
        tissue_samples_mock = [s for s in samples if s.startswith("Cortical_neurons_") and not s.startswith("Cortical_neurons_microglia_") and "mock" in s.lower()]
        tissue_samples_sarscov = [s for s in samples if s.startswith("Cortical_neurons_") and not s.startswith("Cortical_neurons_microglia_") and ("SARS-CoV" in s or "SARS" in s)]
    elif tissue_name == "Cortical_neurons_microglia":
        tissue_samples_mock = [s for s in samples if s.startswith("Cortical_neurons_microglia_") and "mock" in s.lower()]
        tissue_samples_sarscov = [s for s in samples if s.startswith("Cortical_neurons_microglia_") and ("SARS-CoV" in s or "SARS" in s)]
    elif tissue_name == "Microglia":
        tissue_samples_mock = [s for s in samples if s.startswith("Microglia_") and "mock" in s.lower()]
        tissue_samples_sarscov = [s for s in samples if s.startswith("Microglia_") and ("SARS-CoV" in s or "SARS" in s)]
    elif tissue_name == "organoids":
        tissue_samples_mock = [s for s in samples if s.startswith("organoids_mock_")]
        tissue_samples_sarscov = [s for s in samples if s.startswith("organoids_SARS-CoV_")]
    else:
        tissue_samples_mock = []
        tissue_samples_sarscov = []
    
    print(f"  Mock samples: {len(tissue_samples_mock)}")
    print(f"  SARS-CoV samples: {len(tissue_samples_sarscov)}")
    
    # Process each Success cluster
    for idx, cluster_row in success_clusters.iterrows():
        cluster_id = cluster_row['cluster_clean']
        
        # Get max deltapsi
        max_abs_deltapsi = cluster_max_deltapsi.get(cluster_id, pd.NA)
        if not pd.isna(max_abs_deltapsi):
            max_abs_deltapsi = round(max_abs_deltapsi, 2)
        
        # Get introns for this cluster from pre-built mapping
        cluster_introns = cluster_to_introns.get(cluster_id, [])
        
        # Calculate average expression for mock and SARS-CoV
        if len(cluster_introns) > 0:
            subset_df = counts_df.loc[cluster_introns, :]
            
            if len(tissue_samples_mock) > 0:
                # Check if all samples exist in the dataframe
                missing_samples = [s for s in tissue_samples_mock if s not in subset_df.columns]
                if missing_samples:
                    print(f"    Warning: Missing mock samples for cluster {cluster_id}: {missing_samples}")
                    avg_expr_mock = pd.NA
                else:
                    avg_expr_mock = round(subset_df[tissue_samples_mock].values.flatten().mean(), 2)
            else:
                avg_expr_mock = pd.NA
            
            if len(tissue_samples_sarscov) > 0:
                # Check if all samples exist in the dataframe
                missing_samples = [s for s in tissue_samples_sarscov if s not in subset_df.columns]
                if missing_samples:
                    print(f"    Warning: Missing SARS-CoV samples for cluster {cluster_id}: {missing_samples}")
                    avg_expr_sarscov = pd.NA
                else:
                    avg_expr_sarscov = round(subset_df[tissue_samples_sarscov].values.flatten().mean(), 2)
            else:
                avg_expr_sarscov = pd.NA
        else:
            avg_expr_mock = pd.NA
            avg_expr_sarscov = pd.NA
        
        # Store data for this cluster and tissue
        if cluster_id not in cluster_data:
            cluster_data[cluster_id] = {}
        cluster_data[cluster_id][tissue_name] = {
            'p.adjust': round(cluster_row['p.adjust'], 2) if not pd.isna(cluster_row['p.adjust']) else pd.NA,
            'max_abs_deltapsi': max_abs_deltapsi,
            'avg_expr_mock': avg_expr_mock,
            'avg_expr_sarscov': avg_expr_sarscov,
            'genes': cluster_row['genes']
        }

# Now create the wide dataframe
print("\nCreating wide summary table...")

# Get all unique clusters
all_clusters = sorted(cluster_data.keys())

# Prepare data for dataframe
wide_data = []
for cluster_id in all_clusters:
    row = {'cluster': cluster_id}
    
    # Get genes from the first tissue that has this cluster
    genes = None
    for tissue in tissue_folders.keys():
        if tissue in cluster_data[cluster_id]:
            genes = cluster_data[cluster_id][tissue]['genes']
            break
    row['genes'] = genes
    
    # Add data for each tissue
    for tissue_name in tissue_folders.keys():
        if tissue_name in cluster_data[cluster_id]:
            data = cluster_data[cluster_id][tissue_name]
            row[f'{tissue_name}_p.adjust'] = data['p.adjust']
            row[f'{tissue_name}_deltapsi'] = data['max_abs_deltapsi']
            row[f'{tissue_name}_mock_expr'] = data['avg_expr_mock']
            row[f'{tissue_name}_sarscov_expr'] = data['avg_expr_sarscov']
        else:
            row[f'{tissue_name}_p.adjust'] = pd.NA
            row[f'{tissue_name}_deltapsi'] = pd.NA
            row[f'{tissue_name}_mock_expr'] = pd.NA
            row[f'{tissue_name}_sarscov_expr'] = pd.NA
    
    wide_data.append(row)

# Create results dataframe
results_df = pd.DataFrame(wide_data)

# Save results to file
output_file = base_dir / "leafcutter_summary_table_wide.txt"
results_df.to_csv(output_file, sep="\t", index=False)

print(f"\n✓ Wide summary table created: {output_file}")
print(f"\nTable shape: {results_df.shape[0]} rows x {results_df.shape[1]} columns")
print(f"Number of unique clusters: {len(all_clusters)}")
print("\nColumn names:")
print(results_df.columns.tolist())
print("\nFirst few rows:")
print(results_df.head(10))
print("\nColumn names:")
print(results_df.columns.tolist())
print("\nFirst few rows:")
print(results_df.head(10))


