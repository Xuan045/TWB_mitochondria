import pandas as pd
import gzip
import matplotlib.pyplot as plt
import seaborn as sns

INFO_THRESHOLD = 0.7
PROB_THRESHOLD = 0
AF = 0.01

WGS_VCF = '/Users/xuanchou/Documents/TWB_1492/microarray/recode_wgs/TWB1465_wgs_recode.vcf.gz'
FILE_DIR = '/Users/xuanchou/Documents/TWB_1492/microarray/replicate'
IMPUTE_FILE = f'{FILE_DIR}/wgs1465.twb2.imputed.filtered.info{INFO_THRESHOLD}.prob{PROB_THRESHOLD}'
OUTPUT_PLOT = f'{FILE_DIR}/info{INFO_THRESHOLD}_prob{PROB_THRESHOLD}.af_comparison_af{AF}.png'
OUTPUT_TSV = f'{FILE_DIR}/info{INFO_THRESHOLD}_prob{PROB_THRESHOLD}_af_comparison_w_id.tsv'

def parse_vcf(file_path):
    with gzip.open(file_path, 'rb') as file:
        lines = [line.decode('utf-8').strip() for line in file if not line.startswith(b'##')]
    return pd.DataFrame([line.split('\t') for line in lines[1:]], columns=lines[0].split('\t'))

def process_gt(df):
    '''Extract the genotype information from the VCF file and calculate the allele frequency'''
    sample_cols = df.columns[df.columns.str.startswith('NGS')]
    gt_data = {}
    for sample in sample_cols:
        gt_data[sample] = df.apply(lambda x: x[sample].split(':')[0], axis=1)
    
    # Convert the dictionary to a DataFrame
    gt_df = pd.DataFrame(gt_data)

    # Calculate the allele frequency
    gt_df['AF'] = gt_df.apply(cal_af, axis=1)
    
    # Concatenate the variant information
    gt_df = pd.concat([df[['#CHROM', 'POS', 'ID', 'REF', 'ALT']], gt_df], axis=1)
    return gt_df

def cal_af(row):
    '''Calculate the allele frequency for a single row'''
    non_missing_gt = row[5:].apply(lambda gt: gt != './.')  # Exclude missing GT
    ac = (row[5:][non_missing_gt] == '1/1').sum()  # Count homozygous variants
    an = non_missing_gt.sum()  # Count non-missing GT
    return ac / an if an != 0 else 0

def plot_af_comparison(df, x_col, y_col, title, output_file, xlim=(0,1), ylim=(0,1), figsize=(8,6), palette='Set2'):
    """
    Create a scatter plot comparing allele frequencies.

    Parameters:
    data (DataFrame): DataFrame containing the data to be plotted.
    x_col (str): Column name for WGS AF
    y_col (str): Column name for array AF
    title (str): Title of the plot.
    output_file (str): File path to save the plot.
    xlim (tuple): x-axis limits.
    ylim (tuple): y-axis limits.
    figsize (tuple): Figure size.
    palette (str): Color palette for the plot.
    """
    # Calculate the correlation coefficient
    df_clean = df[[x_col, y_col]].dropna()
    df_clean = df_clean[(df_clean[x_col] != 0) & (df_clean[y_col] != 0) & (df_clean[y_col] > AF)]
    n_var = len(df_clean)
    # Pearson correlation coefficient
    pearson_corr = df_clean[x_col].corr(df_clean[y_col])
    
    plt.style.use('default')
    plt.figure(figsize=figsize)
    sns.scatterplot(x=x_col, y=y_col, data=df_clean, palette=palette)

    # Drawing a diagonal line with slope 1
    plt.plot([0, 1], [0, 1], color='red', linestyle='--')

    # Setting labels and limits
    plt.xlabel('WGS_AF')
    plt.ylabel('Array_AF')
    plt.xlim(*xlim)
    plt.ylim(*ylim)

    # Adding a title with the number of points
    plt.title(title)

    # Adding text box with correlation coefficient
    plt.text(0.05, 0.9, f'N = {n_var}\nPearson r = {pearson_corr:.3f}')
    
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    

def main():
    wgs_df = parse_vcf(WGS_VCF)
    print("WGS DataFrame Loaded:", wgs_df.shape)

    wgs_gt_df = process_gt(wgs_df)
    print("Processed WGS DataFrame:", wgs_gt_df.shape)

    # Load imputed data
    impute_df = pd.read_csv(IMPUTE_FILE, sep='\t')
    print("Imputed DataFrame Loaded:", impute_df.shape)

    # Merge two files
    impute_df['merge_cond'] = impute_df['Position'].astype(str) + impute_df['Allele_1'] + impute_df['Allele_0']
    wgs_gt_df['merge_cond'] = wgs_gt_df['POS'].astype(str) + wgs_gt_df['REF'] + wgs_gt_df['ALT']
    merged_df = pd.merge(impute_df, wgs_gt_df, on='merge_cond', how='outer', suffixes=('_impute', '_wgs'))
    # ALT is A0, REF is A1, so ALT_AF is 1 - a1_AF
    merged_df['exp_freq_a0'] = 1.0 - merged_df['exp_freq_a1'].astype(float)
    merged_df['AF_wgs'] = merged_df['AF_wgs'].astype(float)
    print("Merged DataFrame:", merged_df.shape)
    
    # Plot the comparison
    plot_af_comparison(merged_df, x_col='AF_wgs', y_col='exp_freq_a0', title='Comparison of Allele Frequencies', output_file=OUTPUT_PLOT)
    
    # Output the merged DataFrame
    # Reorder the columns
    merged_df.drop(columns=merged_df.columns[merged_df.columns.str.startswith('NGS')], inplace=True)
    non_sample_cols = [col for col in merged_df.columns if not col.startswith('TV2')]
    merged_df = merged_df[non_sample_cols + [col for col in merged_df.columns if col not in non_sample_cols]]
    merged_df.to_csv(OUTPUT_TSV, sep='\t', index=False)

if __name__ == '__main__':
    main()
