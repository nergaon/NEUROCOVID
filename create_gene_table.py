#!/usr/bin/env python3
import pandas as pd
import numpy as np

# Gene information
genes_data = {
    'ACE2': 'ENSG00000130234',
    'TMPRSS2': 'ENSG00000184012',
    'Furin': 'ENSG00000140564',
    'Integrin subunit beta1 - ITGB1': 'ENSG00000150093',
    'Integrin subunit beta3 - ITGB3': 'ENSG00000259207',
    'GRP78 (HSPA5)': 'ENSG00000044574',
    'DPP4': 'ENSG00000197635',
    'AXL': 'ENSG00000167601',
    'CD147 (BSG)': 'ENSG00000172270',
    'NRP1': 'ENSG00000099250',
    'CD209': 'ENSG00000090659',
    'CD209L (CLEC4M)': 'ENSG00000104938',
    'Heparan sulfate proteoglycan (HSPG) (SDC2)': 'ENSG00000169439',
    'Vimentin (VIM)': 'ENSG00000026025'
}

# Read the expression table (skip the first comment line)
df = pd.read_csv('/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/expression_table/expression_table_HN1.txt', 
                  sep='\t', skiprows=1, nrows=None)

# Get sample column names (exclude metadata columns)
sample_cols = [col for col in df.columns if col.endswith('.bam')]

# Define tissue types
tissue_types = ['Cortical_neurons', 'organoids', 'Cortical_neurons_microglia', 'Microglia']

# Classify samples by tissue type and condition
sample_groups = {}
for tissue in tissue_types:
    if tissue == 'Cortical_neurons':
        # Exclude Cortical_neurons_microglia when looking for plain Cortical_neurons
        sample_groups[tissue] = {
            'mock': [col for col in sample_cols if 'Cortical_neurons' in col and 'Cortical_neurons_microglia' not in col and 'mock' in col],
            'SARS-CoV': [col for col in sample_cols if 'Cortical_neurons' in col and 'Cortical_neurons_microglia' not in col and 'SARS-CoV' in col]
        }
    else:
        sample_groups[tissue] = {
            'mock': [col for col in sample_cols if tissue in col and 'mock' in col],
            'SARS-CoV': [col for col in sample_cols if tissue in col and 'SARS-CoV' in col]
        }

# Print summary
for tissue, groups in sample_groups.items():
    print(f"{tissue}: {len(groups['mock'])} mock, {len(groups['SARS-CoV'])} SARS-CoV")

# Create output table
output_data = []

for gene_name, gene_id in genes_data.items():
    row = {'Gene name': gene_name, 'Ensembl ID': gene_id}
    
    # Find the row for this gene
    gene_row = df[df['Geneid'] == gene_id]
    
    if len(gene_row) > 0:
        gene_row = gene_row.iloc[0]
        
        # Calculate average expression for each tissue type and condition
        for tissue in tissue_types:
            mock_cols = sample_groups[tissue]['mock']
            sars_cols = sample_groups[tissue]['SARS-CoV']
            
            if mock_cols:
                mock_values = gene_row[mock_cols].astype(float)
                row[f'{tissue}_mock_avg'] = int(round(mock_values.mean(), 0))
            else:
                row[f'{tissue}_mock_avg'] = np.nan
            
            if sars_cols:
                sars_values = gene_row[sars_cols].astype(float)
                row[f'{tissue}_SARS-CoV_avg'] = int(round(sars_values.mean(), 0))
            else:
                row[f'{tissue}_SARS-CoV_avg'] = np.nan
    else:
        print(f"Warning: {gene_id} ({gene_name}) not found in expression table")
        for tissue in tissue_types:
            row[f'{tissue}_mock_avg'] = np.nan
            row[f'{tissue}_SARS-CoV_avg'] = np.nan
    
    output_data.append(row)

# Create output dataframe
output_df = pd.DataFrame(output_data)

# Save to file
output_file = '/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/expression_table/gene_expression_summary_HN1.txt'
output_df.to_csv(output_file, sep='\t', index=False)

print(f"\nOutput saved to: {output_file}")
print("\nSummary table:")
print(output_df.to_string(index=False))
