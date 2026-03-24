# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 11:11:13 2020

@author: hadas
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
'''
plot count and PSI bar from all the clusters in the JSR counts table
'''

def plot_bar_2(persent_cluster, value_cluster, fig_title, main_dir, legend_color_map):
    '''
    input: datafaram with the clusters values in percent
    output: bar plot of the introns
    '''
    output_1 = main_dir + "genes_figs/" + fig_title + '.svg' #png, pdf, jpeg, , tiff, raw, bmp, eps or svg
    output_2 = main_dir + "genes_figs/" + fig_title + '.png'
    #make sure the junctions in the value and % tables are in the same order
    persent_cluster = persent_cluster.sort_index()
    value_cluster = value_cluster.sort_index()
    #get introns position as the legends
    intron_legend = persent_cluster.index
    intron_legend  = intron_legend.tolist() 
    legend_color_map = create_color_map(intron_legend, legend_color_map)

    #drop only the trailing replicate id from sample names (e.g. ..._14)
    persent_cluster.columns = ['_'.join(col.split('_')[:-1]) for col in persent_cluster.columns]
    value_cluster.columns = ['_'.join(col.split('_')[:-1]) for col in value_cluster.columns]
    index = []
    for i in range(0, len(persent_cluster.columns)-1):      
        if persent_cluster.columns[i] != persent_cluster.columns[i+1]:
            index.append(i+1)
    #add col with 0 to seperate different samples
    zeroCol = []
    runNum = 0
    for i in index:
        col_name = str(i)
        i = i + runNum
        persent_cluster.insert(i, col_name, 0)
        value_cluster.insert(i, col_name, 0)
        runNum += 1   
        zeroCol.append(i)
    #the names in x axis
    labels = persent_cluster.columns
    labels = labels.tolist()
    labels = [
        'mock' if 'mock' in str(name).lower() else ('SARS-CoV' if 'sars-cov' in str(name).lower() else name)
        for name in labels
    ]
    unique_names = []
    for i in range(len(labels)):
        if i in zeroCol:
            labels[i] = ""
        #elif i%5 != 1: #second sample will have the name
        elif labels[i] in unique_names: #first sample will have the name
            labels[i] = ' '
        else:
            unique_names.append(labels[i])
    N = int(len(labels)) 
    r = np.zeros(shape=(N))  # the x locations for the groups. 
    for i in range(0, N-1):
        r[i+1] = r[i] + 0.25 
    barWidth = 0.2  # the width of the bars 
    plt.close('all')
    #my_title = gene
    introns_persent = persent_cluster.values
    introns_counts = value_cluster.values
    bars_persent = np.zeros(shape=(N))
    bars_counts = np.zeros(shape=(N))
    # Create subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4, 4), gridspec_kw={'height_ratios': [1, 1]})

    for i in range(0, len(introns_persent)):
        color = legend_color_map[intron_legend[i]]
        ax1.bar(r, introns_persent[i], width = barWidth, bottom=bars_persent, edgecolor = "none", color=color)
        bars_persent = np.add(bars_persent, introns_persent[i])
        ax2.bar(r, introns_counts[i], width = barWidth, bottom=bars_counts, edgecolor = "none", color=color)
        bars_counts = np.add(bars_counts, introns_counts[i])
    
    ax1.set_xticks(r) 
    ax1.set_xticklabels(labels, fontsize=0)  #no cell names in the upper fig
    ax1.tick_params(axis="x", which='both', length=0) 
    ax1.set_ylim(0, 1)
    ax1.set_yticks([0, 1]) 
    ax1.tick_params(axis="y", labelsize=10) 
    ax1.set_ylabel('PSI')
    
    ax2.set_xticks(r) 
    ax2.tick_params(axis="x", which='both', length=0) 
    #ax2.set_xticklabels(labels, fontsize=10, rotation=30) 
    ax2.set_xticklabels(labels, fontsize=10, rotation=0)
    ax2.tick_params(axis="y", labelsize=10)
    #ax2.legend(intron_legend, fontsize=8, markerscale=0.15)
    ax2.set_ylabel('Counts')
    # Add the legend below both plots
    fig.legend(intron_legend, fontsize=8, markerscale=0.15, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=1)
    fig.suptitle(fig_title, fontsize=12)
    # Adjust layout to make space for the legend
    plt.subplots_adjust(hspace=0.4)
    #plt.show()
    plt.savefig(output_1, bbox_inches='tight')
    plt.savefig(output_2, bbox_inches='tight')
    plt.close('all')
    return (legend_color_map)

def create_color_map(intron_legend, legend_color_map):
    """
    Create or update a color map for intron legends.
    - Keeps existing colors.
    - Assigns new colors only to new introns.
    """
    #legend_color_map = {}
    cmap = plt.get_cmap("tab20")  # or any other colormap
    color_index = 0
    used_color = [] #in the same cluster, don't use the same colors
    for intron in intron_legend:
        if intron in legend_color_map:
           used_color.append(legend_color_map[intron]) 
    for intron in intron_legend:
        if intron not in legend_color_map:
            while cmap(color_index % cmap.N) in used_color:
               color_index += 1 
            # Find the next unused color
            legend_color_map[intron] = cmap(color_index % cmap.N)
            color_index += 1
    return legend_color_map

def main():
    cell_type = ['Cortical_neurons_microglia', 'Cortical_neurons', 'Microglia', 'organoids']
    main_dir = '/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/leafcutter_0.2.9/'
    output_dir = os.path.join(main_dir, 'genes_figs')
    if os.path.isdir(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    # Load the full counts file once
    # Format: space-separated; first row = sample names; first col = junction (chr:start:end:clu_N_strand)
    counts_df = pd.read_csv(main_dir + 'NEUROCOVID_perind_numers.counts', sep=' ', index_col=0)

    # Parse junction index into h_junction (chr:start:end) and cluster id (chr:clu_N_strand)
    junc_index = counts_df.index.tolist()
    h_junctions = [j.rsplit(':', 1)[0] for j in junc_index]           # e.g. 'chr1:15947:16607'
    clu_parts   = [j.rsplit(':', 1)[1] for j in junc_index]           # e.g. 'clu_1_-'
    chroms      = [j.split(':')[0]     for j in junc_index]           # e.g. 'chr1'
    clusters    = [c + ':' + p for c, p in zip(chroms, clu_parts)]    # e.g. 'chr1:clu_1_-'
    counts_df.index.name = 'junctions'  # rename index to 'cluster'
    counts_df['cluster'] = clusters
    legend_color_map = {}

    for one_cell in cell_type:
        print(one_cell)
        cell_dir = main_dir + one_cell + '/'

        # Get sample columns for this cell type
        group_df = pd.read_csv(cell_dir + 'group_file.txt', sep=' ', header=None, names=['sample', 'condition'])
        sample_cols = group_df['sample'].tolist()

        # value_table: counts for this cell type's samples + cluster column
        value_table = counts_df[sample_cols + ['cluster']].copy()

        # percent_table: PSI = junction count / total cluster count per sample
        percent_table = value_table.copy()
        for sample in sample_cols:
            cluster_totals = value_table.groupby('cluster')[sample].transform('sum')
            percent_table[sample] = value_table[sample].div(cluster_totals).fillna(0)

        # Success clusters from significance file
        sig_df = pd.read_csv(cell_dir + 'leafcutter_ds_cluster_significance.txt', sep='\t')
        sucess_clusters = sig_df.loc[sig_df['status'] == 'Success', 'cluster'].tolist()

        count_fig = 0
        for cluster in sucess_clusters:
            print(cluster)
            count_fig += 1
            persent_cluster = percent_table.loc[percent_table['cluster'] == cluster].copy()
            #persent_cluster = persent_cluster.sort_index()
            row = sig_df.loc[sig_df['cluster'] == cluster].iloc[0]
            gene_name = row['genes']
            p_adjust = row['p.adjust']
            formatted_p = '{:.3f}'.format(p_adjust) if pd.notna(p_adjust) else 'NA'
            persent_cluster.drop(columns=['cluster'], inplace=True)
            value_cluster = value_table.loc[value_table['cluster'] == cluster].copy()
            value_cluster.drop(columns=['cluster'], inplace=True)
            try:
                fig_title = one_cell + '_' + str(gene_name) + '_' + cluster + '_p=' + formatted_p
            except Exception:
                fig_title = one_cell + '_' + cluster + '_p=' + formatted_p
            legend_color_map = plot_bar_2(persent_cluster, value_cluster, fig_title, main_dir, legend_color_map)
        print(count_fig, 'success clusters in', one_cell)

    return
    
if __name__ == "__main__":
    main()   
    
