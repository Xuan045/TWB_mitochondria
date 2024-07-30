import pandas as pd
import re

# The script reads a file containing gene information and extracts the gene name and positions
# It then determines the non-coding regions in the sequence and saves the information to a tab-separated file

gene_positions = []
seq_file='/Users/xuanchou/Documents/TWB_1492/Analysis/seq_info/NC_012920.1_sequence.txt'

with open(seq_file, 'r') as file:
    for line in file:
        if line.startswith(">"):
            # Extract gene name
            gene_name = None
            gene_match = re.search(r'\[gene=([^]]+)\]', line)
            if gene_match:
                gene_name = gene_match.group(1)
            
            # Extract positions
            start_pos = None
            end_pos = None
            location_match = re.search(r'\[location=(complement\()*(\d+)\.\.(\d+)\]*', line)
            if location_match:
                complement = location_match.group(1)
                start_pos = location_match.group(2)
                end_pos = location_match.group(3)
            
            if gene_name and start_pos and end_pos:
                gene_positions.append([gene_name, start_pos, end_pos, 'Y' if complement else 'N'])

# Convert the list to a DataFrame
gene_positions_df = pd.DataFrame(gene_positions, columns=['Gene Name', 'Start Position', 'End Position', 'Complement'])

# Define the full range
full_range = set(range(1, 16570))

# Extract the gene regions from the DataFrame
gene_ranges = set()
for _, row in gene_positions_df.iterrows():
    gene_ranges.update(range(int(row['Start Position']), int(row['End Position']) + 1))

# Determine the non-coding regions
non_coding_ranges = full_range - gene_ranges

# Convert the non-coding ranges to start-end format
non_coding_positions = []
start = None
end = None

for i, pos in enumerate(sorted(non_coding_ranges)):
    # If start is None, this is the start of a new interval
    if start is None:
        start = pos
        end = pos

    # Check if the current position is contiguous with the previous position
    elif pos - 1 == end:
        end = pos

    # If not contiguous or if we have reached the last element in the list, save the current interval
    if pos + 1 not in non_coding_ranges or i == len(non_coding_ranges) - 1:
        non_coding_positions.append(['non-coding', start, end, pd.NA])
        start = None
        end = None

# Convert the non-coding positions list to a DataFrame and append it to the original DataFrame
non_coding_df = pd.DataFrame(non_coding_positions, columns=['Gene Name', 'Start Position', 'End Position', 'Complement'])
combined_df = pd.concat([gene_positions_df, non_coding_df])

# Convert the 'Start Position' and 'End Position' columns to integers for proper sorting
combined_df['Start Position'] = combined_df['Start Position'].astype(int)
combined_df['End Position'] = combined_df['End Position'].astype(int)

# Sort the combined DataFrame by 'Start Position'
combined_df = combined_df.sort_values(by='Start Position').reset_index(drop=True)

# Add gene function
cr_range = sorted(full_range - set(range(577, 16024)))

# Create an empty list to store the function information
functions = []

for i, row in combined_df.iterrows():
    if row['Start Position'] in cr_range and row['End Position'] in cr_range:
        functions.append('Control region')
    elif row['Gene Name'].startswith('TRN'):
        functions.append('tRNA')
    elif row['Gene Name'].startswith('RNR'):
        functions.append('rRNA')
    elif row['Gene Name'] == 'non-coding':
        functions.append('Intergenic region')
    else:
        functions.append('Protein coding')  # If no condition is met, add None to indicate no function
combined_df['Function'] = pd.DataFrame(functions)

col = combined_df.pop('Complement')
combined_df.insert(len(combined_df.columns), 'Complement', col)

# Save the combined DataFrame to a tab-separated file
output_file_path_combined = '/Users/xuanchou/Documents/TWB_1492/Analysis/seq_info/combined_positions.tsv'
combined_df.to_csv(output_file_path_combined, sep='\t', index=False)
