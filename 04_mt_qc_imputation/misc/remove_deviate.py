#!/usr/bin/python

import pandas as pd
import gzip
from sys import argv

def parse_vcf(vcf_file):
    # Check if the file is gzipped based on its extension
    if vcf_file.endswith('.gz'):
        open_func = gzip.open
    else:
        open_func = open
        
    with open_func(vcf_file, 'rb') as file:
        lines = [line.decode('utf-8').strip() for line in file if not line.startswith(b'##')]
    return pd.DataFrame([line.split('\t') for line in lines[1:]], columns=lines[0].split('\t'))

def parse_gcount(plink_file):
    gcount_df = pd.read_csv(f'{plink_file}.gcount', sep='\t')
    # assign column names while loading bim file
    bim_df = pd.read_csv(f'{plink_file}.bim', sep='\t', header=None, names=['CHROM', 'ID', 'cM', 'POS', 'ALT', 'REF'])
    gcount_df = pd.merge(gcount_df, bim_df, on=['REF', 'ALT', 'ID'])
    # Calculate AF
    gcount_df['AF'] = gcount_df['HAP_ALT_CTS'] / (gcount_df['HAP_REF_CT'] + gcount_df['HAP_ALT_CTS'])
    return gcount_df

def parse_gt(df):
    # Get GT index
    gt_idx = df.iloc[1]['FORMAT'].split(':').index('GT')
    # Get sample columns (all after 'FORMAT')
    sample_cols = df.columns[df.columns.get_loc('FORMAT')+1:]

    # Initialize an empty list to store results
    results = []

    # Process each sample column
    for col in sample_cols:
        # Check each cell in the column
        for idx, cell in df[col].items():
            if pd.notnull(cell):
                gt = cell.split(':')[gt_idx]
                # Keep only GT part
                df.at[idx, col] = gt
                if gt in ['0|1', '1|0']:
                    # Store the sample, variant, and gt
                    results.append({'Sample': col, 'Variant': df.iloc[idx]['#CHROM'] + ':' + str(df.iloc[idx]['POS']), 'GT': gt})
                    # Convert the genotype to '.|.'
                    df.at[idx, col] = '.|.'

    return df, results

def gt_to_dosage(df):
    # Identify columns containing genotype information
    gt_columns = df.columns[df.columns.get_loc('FORMAT')+1:]

    # Define a function to convert GT to dosage
    def to_dosage(gt):
        # Convert GT to dosage
        if gt in ['0/0', '0|0']:
            return 0
        elif gt in ['1/1', '1|1']:
            return 2
        else:
            return float('nan') # 0/1 or 1/0 should be Nan

    # Apply the conversion for each sample column
    for col in gt_columns:
        df[col] = df[col].apply(to_dosage)

    return df

def calculate_AF(df):
    df = parse_gt(df)[0]
    df = gt_to_dosage(df)
    # Count samples with dosage 2 (AC) and non-missing values (AN)
    AC = (df.iloc[:, 9:] == 2).sum(axis=1)
    AN = df.iloc[:, 9:].notna().sum(axis=1)
    AF = AC / AN
    return AF

def find_deviate(wgs_vcf, gcount_df):
    # POS column to string
    wgs_vcf['POS'] = wgs_vcf['POS'].astype(str)
    gcount_df['POS'] = gcount_df['POS'].astype(str)
    # Merge dataframes on ID
    merged_df = pd.merge(wgs_vcf, gcount_df, on=['POS', 'REF', 'ALT'], suffixes=('_wgs', '_gcount'), how='right')
    
    # Calculate the absolute difference in AF
    merged_df['AF_diff'] = abs(merged_df['AF_wgs'] - merged_df['AF'])
    
    # Filter criteria:
    # 1. Keep variants with AF > 0.
    # 2. Keep variants where AF_diff <= 0.2 or AF_wgs is NA (cuz some variants are only in array callset).
    keep_df = merged_df[(merged_df['AF'] > 0) & ((merged_df['AF_diff'] <= 0.2) | merged_df['AF_wgs'].isna())]
    
    # Return cleaned dataframe
    return keep_df

def main():
    WGS_FILE = argv[1]
    plink_file = argv[2]
    
    wgs_vcf = parse_vcf(WGS_FILE)
    gcount_df = parse_gcount(plink_file)
    
    wgs_vcf['AF_wgs'] = calculate_AF(wgs_vcf)
    
    # Find variants with AF difference > 0.2
    gcout_filtered_df = find_deviate(wgs_vcf, gcount_df)

    for id in gcout_filtered_df['ID_gcount']:
        print(id)
    
if __name__ == "__main__":
    main()
    