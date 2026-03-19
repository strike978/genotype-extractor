#!/usr/bin/env python3
import matplotlib.pyplot as plt
import sys
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform

# Function to check if DataFrame is empty


def check_empty(df):
    if df.empty:
        print('No IBS or GBR samples found.')
        sys.exit(1)

# Function to infer major/minor allele and encode as dosage (0,1,2)


def encode_allele_dosage(col):
    alleles = col.dropna().astype(str).str.cat()
    allele_counts = pd.Series(list(alleles)).value_counts()
    if len(allele_counts) < 2:
        return pd.Series([np.nan]*len(col), index=col.index)
    major, minor = allele_counts.index[0], allele_counts.index[1]

    def dosage(gt):
        if pd.isna(gt) or len(gt) != 2:
            return np.nan
        return sum(1 for x in gt if x == minor)
    return col.apply(dosage)


# --- Remove outliers from each group in 1000g and HGDP and save new CSV ---
print("\nRemoving outliers from each group in 1000g and HGDP and saving new CSV...")
df_all = pd.read_csv('genotypes.csv')
cleaned_rows = []
sources_to_filter = ['1000g', 'HGDP']
for source in sources_to_filter:
    df_source = df_all[df_all['source'] == source]
    for group in df_source['group'].unique():
        df_group = df_source[df_source['group'] == group]
        if df_group.empty:
            continue
        genotype_cols = df_group.columns[4:]
        genotypes = df_group[genotype_cols]
        encoded = genotypes.apply(encode_allele_dosage)
        encoded = encoded.dropna(axis=1, thresh=int(0.9*len(encoded)))
        encoded = encoded.fillna(encoded.mean())
        if len(encoded) < 2:
            # Not enough samples to compute distances
            cleaned_rows.append(df_group)
            continue
        dist_matrix = squareform(pdist(encoded.values, metric='euclidean'))
        avg_dist = dist_matrix.mean(axis=1)
        q1 = np.percentile(avg_dist, 25)
        q3 = np.percentile(avg_dist, 75)
        iqr = q3 - q1
        iqr_threshold = q3 + 1.5 * iqr
        outlier_mask = avg_dist > iqr_threshold
        non_outlier_df = df_group.loc[~outlier_mask].copy()
        cleaned_rows.append(non_outlier_df)
# Add all rows from sources not in sources_to_filter (e.g., SGDP) without filtering
other_sources = df_all[~df_all['source'].isin(sources_to_filter)]
if not other_sources.empty:
    cleaned_rows.append(other_sources)
df_cleaned = pd.concat(cleaned_rows, ignore_index=True)
df_cleaned.to_csv('genotypes_no_outliers.csv', index=False)
print("Saved cleaned CSV as 'genotypes_no_outliers.csv'.")
