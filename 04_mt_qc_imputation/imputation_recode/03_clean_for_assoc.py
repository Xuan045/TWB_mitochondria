import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

INFO_THRESHOLD = 0.7
PROB_THRESHOLD = 0
AF = 0.01
FILE = f'info{INFO_THRESHOLD}_prob{PROB_THRESHOLD}_af_comparison_w_id.tsv'

# Load the data
impute_df = pd.read_csv(FILE, sep='\t')

# Drop unnecessary columns
impute_df.drop(['Chromosome', 'rs_id', 'position', 'a0', 'a1', 'exp_freq_a1', 'info_type0', 'certainty', 'type', 'concord_type0', 'r2_type0', 'AN', 'AC', 'merge_cond', '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'AF_impute'], axis=1, inplace=True)

# Rename columns
impute_df.rename(columns={'Allele_0': 'ALT', 'Allele_1': 'REF', 'exp_freq_a0': 'AF_array'}, inplace=True)

# Recode sample columns
sample_cols = [col for col in impute_df.columns if col.startswith('TV2')]

def convert_values(val):
    if val == '0/0':
        return 0
    elif val == '1/1':
        return 2
    else:
        return None

for col in sample_cols:
    impute_df[col] = impute_df[col].apply(convert_values)

# Filter data
impute_df_filter = impute_df[(impute_df['AF_array'] > AF) & (impute_df['AF_array'] < (1 - AF))]
impute_df_filter_wgs = impute_df_filter[abs(impute_df_filter['AF_array'] - impute_df_filter['AF_wgs']) < 0.1]

# Identify columns that are not sample columns
non_sample_cols = [col for col in impute_df_filter_wgs.columns if not col.startswith('TV2')]

# Merge data
merged_df = pd.merge(impute_df_filter[non_sample_cols], impute_df_filter_wgs[non_sample_cols], how='outer', indicator=True)

# Filter out data
diff_df = merged_df[merged_df['_merge'] == 'left_only']

# Manually add specific rows
manually_add_df = impute_df_filter[impute_df_filter['RS_ID'].isin(['rs527236189', 'Affx-157984095'])]
impute_df_filter_wgs = pd.concat([impute_df_filter_wgs, manually_add_df])

# Format the output summary file for variants
impute_df_out = impute_df_filter_wgs[non_sample_cols]
impute_df_out.rename(columns={'snp_id': 'status', 'RS_ID': 'ID'}, inplace=True)
impute_df_out['status'] = np.where(impute_df_out['status'] == 'chrM', 'genotyped', 'imputed')

impute_df_out = impute_df_out[['ID', 'Position', 'ALT', 'REF', 'status', 'info', 'AF_array', 'AF_wgs']].reset_index(drop=True)
impute_df_out['Position'] = impute_df_out['Position'].astype(int)
impute_df_out['AF_wgs'] = impute_df_out['AF_wgs'].apply(lambda x: f'{x:.3f}')
impute_df_out['AF_array'] = impute_df_out['AF_array'].apply(lambda x: f'{x:.3f}')
impute_df_out.to_csv('twb2_assoc_var_summary.tsv', sep='\t', index=False)

# Reformat for association
impute_df_filter_wgs['Position'] = impute_df_filter_wgs['Position'].astype(int).astype(str)
impute_df_filter_wgs['var'] = impute_df_filter_wgs['Position'] + ':' + impute_df_filter_wgs['REF'] + ':' + impute_df_filter_wgs['ALT']

# Convert DataFrame from wide to long format
impute_melted_df = impute_df_filter_wgs.melt(id_vars=['var'], value_vars=sample_cols, var_name='sample', value_name='genotype')
impute_melted_df['genotype'] = impute_melted_df['genotype'].astype(str)

# Pivot data
impute_pivot_df = impute_melted_df.pivot(index='sample', columns='var', values='genotype')
impute_pivot_df.to_csv(f'impute_for_assoc_info{INFO_THRESHOLD}_prob{PROB_THRESHOLD}_af{AF}_EAS_kingMT.txt', sep='\t')

# Plot allele frequency comparison
def plot_af_comparison(df, x_col, y_col, title, xlim=(0, 1), ylim=(0, 1), figsize=(8, 6), palette='Set2'):
    df_clean = df[[x_col, y_col]].dropna()
    pearson_corr = df_clean[x_col].corr(df_clean[y_col])
    
    df['color_group'] = f'Overlapped (N = {len(df_clean)})'
    df.loc[df[x_col].isna(), 'color_group'] = f'Not in WGS (N = {len(df) - len(df_clean)})'
    df[x_col].fillna(0, inplace=True)

    plt.style.use('default')
    plt.figure(figsize=figsize)
    
    scatter = sns.scatterplot(x=x_col, y=y_col, hue='color_group', data=df, palette=palette, legend='brief')
    plt.legend(title=None)

    plt.plot([0, 1], [0, 1], color='red', linestyle='--')
    plt.xlabel('WGS_AF')
    plt.ylabel('Array_AF')
    plt.xlim(*xlim)
    plt.ylim(*ylim)
    plt.title(title)
    plt.text(0.05, 0.8, f'N = {len(df)}\nPearson r = {pearson_corr:.3f} (Overlapped)')
    plt.savefig(f'af_final_info{INFO_THRESHOLD}_prob{PROB_THRESHOLD}_af{AF}.png', dpi=300, bbox_inches='tight')

# Example usage
plot_af_comparison(impute_df_filter_wgs, 'AF_wgs', 'AF_array', 'AF comparison')
