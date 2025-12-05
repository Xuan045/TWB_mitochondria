import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Define sets and their corresponding colors
set1 = ['m.489C', 'm.8701G', 'm.9540C', 'm.10398G', 'm.10400T', 'm.10873C', 'm.14783C', 'm.15043A']
set2 = ['m.827G', 'm.6023A', 'm.6216C', 'm.6413C', 'm.16217C', 'm.15535T']

color_palette = {
    variant: '#E68F74' for variant in set1
}
color_palette.update({variant: '#6A92CD' for variant in set2})

# Load SNP label data
bim_data = pd.read_csv('info0.7_extract.bim', sep="\t", header=None, usecols=[3, 5])
bim_data.columns = ['Position', 'Allele']
bim_data['SNP_Label'] = 'm.' + bim_data['Position'].astype(str) + bim_data['Allele']
snp_labels = bim_data['SNP_Label'].tolist()
# save the snp_labels to a file
with open('snp_labels.txt', 'w') as f:
    for item in snp_labels:
        f.write("%s\n" % item)

# Load LD matrix data
ld_matrix = np.loadtxt('info0.7_extract.ld', delimiter='\t')

# Create the heatmap without dendrogram, using the custom color palette
clustermap_fig = sns.clustermap(ld_matrix, cmap='Greys', vmin=0, vmax=1, figsize=(10, 10),
                                xticklabels=snp_labels, yticklabels=snp_labels,
                                row_colors=[color_palette.get(x, '#ffffff') for x in snp_labels],
                                col_colors=[color_palette.get(x, '#ffffff') for x in snp_labels])

# Set a title for the colorbar
clustermap_fig.cax.set_title('r²')

# Rotate the labels for better readability
plt.setp(clustermap_fig.ax_heatmap.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
plt.setp(clustermap_fig.ax_heatmap.get_yticklabels(), rotation=0)

# Set text color and make it bold based on variant label
for tick_label in clustermap_fig.ax_heatmap.get_xticklabels():
    variant_label = tick_label.get_text()
    tick_label.set_color(color_palette.get(variant_label, 'black'))
    tick_label.set_fontweight('bold')

for tick_label in clustermap_fig.ax_heatmap.get_yticklabels():
    variant_label = tick_label.get_text()
    tick_label.set_color(color_palette.get(variant_label, 'black'))
    tick_label.set_fontweight('bold')

# Increase font size for better visibility
clustermap_fig.ax_heatmap.tick_params(labelsize=12)
clustermap_fig.ax_heatmap.set_xlabel(clustermap_fig.ax_heatmap.get_xlabel(), fontsize=18)
clustermap_fig.ax_heatmap.set_ylabel(clustermap_fig.ax_heatmap.get_ylabel(), fontsize=18)

# Save the plot with larger text as "ld_heatmap_poster.png"
plt.savefig('ld_heatmap_poster.png', dpi=300, bbox_inches='tight')
