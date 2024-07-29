import pandas as pd

# Define the thresholds for filtering
INFO_THRESHOLD = 0.7
PROB_THRESHOLD = 0

FILE_DIR = '/Users/xuanchou/Documents/TWB_1492/microarray/replicate'
PARA = 'twb2.EAS_kingMT.shapeit2.impute2'
IMPUTED_FILE = f'{FILE_DIR}/{PARA}.diploid'
SAMPLE_FILE = f'{FILE_DIR}/twb2.shapeit2.sample'
INFO_FILE = f'{FILE_DIR}/{PARA}.diploid_info'
OUTPUT_FILE = f'{FILE_DIR}/wgs1465.twb2.imputed.filtered.info{INFO_THRESHOLD}.prob{PROB_THRESHOLD}'


def load_samples(file_path):
    """Load sample information from .sample file."""
    with open(file_path, 'r') as f:
        # Skip the first two lines
        next(f)
        next(f)
        samples = [line.split()[0] for line in f]
    return samples

def process_line(line, prob_threshold):
    """
    Function to process a single line of the imputed file
    Input: a line of the imputed file
    Output: a list of the processed line, with the first five components as variant identifiers
    """
    parts = line.strip().split(' ')
    identifiers = parts[:5]  # Variant identifiers
    
    processed_probs = []
    for i in range(5, len(parts), 3):
        prob_set = [float(parts[j]) for j in range(i, i + 3)]
        max_prob = max(prob_set)
        if max_prob > prob_threshold:
            # Determine the genotype based on the index of the maximum probability
            genotype = '1/1' if prob_set.index(max_prob) == 0 else ('0/0' if prob_set.index(max_prob) == 2 else None)
            processed_probs.append(genotype)
        else:
            processed_probs.append(None)

    return identifiers + processed_probs

def process_imputed_file(file_path, prob_threshold):
    """Process the entire imputed file."""
    with open(file_path, 'r') as file:
        processed_data = [process_line(line, prob_threshold) for line in file]
    return processed_data

def main():
    samples = load_samples(SAMPLE_FILE)
    processed_data = process_imputed_file(IMPUTED_FILE, PROB_THRESHOLD)
    processed_df = pd.DataFrame(processed_data, columns=['Chromosome', 'RS_ID', 'Position', 'Allele_0', 'Allele_1'] + samples)
    # processed_df.to_csv(f'{FILE_DIR}/wgs1465.twb2.imputed', sep='\t', index=False)
    print('Finished processing imputed file.')
    
    info_df = pd.read_csv(INFO_FILE, sep=' ')
    
    # Convert columns to string for merging
    for col in ['Chromosome', 'RS_ID', 'Position', 'Allele_0', 'Allele_1']:
        processed_df[col] = processed_df[col].astype(str)
    for col in ['snp_id', 'rs_id', 'position', 'a0', 'a1']:
        info_df[col] = info_df[col].astype(str)
    
    merged_df = processed_df.merge(info_df, left_on=['Chromosome', 'RS_ID', 'Position', 'Allele_0', 'Allele_1'], right_on=['snp_id', 'rs_id', 'position', 'a0', 'a1'], how='outer')
    # merged_df = merged_df[[col for col in merged_df.columns if col not in ['snp_id', 'rs_id', 'position', 'a0', 'a1']] + samples]

    final_df = merged_df[merged_df['info'] > INFO_THRESHOLD]
    print(f'Finished filtering variants. {len(final_df)} variants remain.')
    
    # Calculate AC, AN, and AF
    final_df['AN'] = final_df[samples].apply(lambda x: x.notnull().sum(), axis=1)
    final_df['AC'] = final_df[samples].apply(lambda x: x[x == '1/1'].count(), axis=1)
    final_df['AF'] = final_df['AC'] / final_df['AN']
    print('Finished calculating AC, AN, and AF.')
    
    # Reorder the columns
    non_sample_cols = [col for col in final_df.columns if col not in samples]
    final_df = final_df[non_sample_cols + samples]
    
    final_df.to_csv(OUTPUT_FILE, sep='\t', index=False)
    print(f'Output written to {OUTPUT_FILE}.')

if __name__ == '__main__':
    main()
