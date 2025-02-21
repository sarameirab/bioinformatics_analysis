import pandas as pd
from urls import *
from utils import *





# ...rest of your code remains the same...

# Download and read file using gdown
path = download_url(PROTEOMICS_CORRELATION_URL)
depmap_gdscs_merged_ttest_path = download_url(DEPMAP_GDSC1_TTEST_ALL_MERGED)

# Read the downloaded CSV
proteomics_df = pd.read_csv(path)

# Get unique protein IDs from your data
protein_ids = proteomics_df['protein_id'].unique().tolist()

# Map proteins to genes
protein_to_gene = map_uniprot_to_gene(protein_ids)

# Add gene names to your DataFrame
proteomics_df['gene_name'] = proteomics_df['protein_id'].map(protein_to_gene)


pd.set_option('display.max_columns', 50)
pd.set_option('display.width', None)

# Display first 10 rows where cell_line_count > 2
print(proteomics_df.loc[proteomics_df['cell_line_count'] > 2].head(10))

depmap_gdscs_merged_ttest_path = pd.read_csv(depmap_gdscs_merged_ttest_path)
depmap_gdscs_merged_ttest_path['clean_gene'] = depmap_gdscs_merged_ttest_path['mutation'].str.extract(r'(.*?)\s*\(\d+\)')[0]

proteomics_df.columns = ['proteomics_' + col if col != 'gene_name' else col for col in proteomics_df.columns]


# Perform the merge
merged_df = proteomics_df.merge(
    depmap_gdscs_merged_ttest_path,
    left_on=['gene_name', 'proteomics_drug_id'],
    right_on=['clean_gene', 'gdsc1_drug'],
    how='inner'
)


# Display results
print("\nFirst few rows of merged dataframe:")
merged_df = merged_df.drop(['clean_gene'], axis=1)
print(merged_df.head())
print(f"DataFrame dimensions (rows, columns): {merged_df.shape}")



filtered_df = merged_df[
    ((merged_df['proteomics_correlation'] > 0.2) | (merged_df['proteomics_correlation'] < -0.2)) &
    (merged_df['gdsc1_p_value'] < 0.1) &
    (merged_df['depmap_p_value'] < 0.1) & (merged_df['proteomics_p_value'] < 0.1)
]


# Set display options for better visibility
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
print("\nFiltered results:")
print(filtered_df)
print(f"\nNumber of rows meeting criteria: {len(filtered_df)}")

# Write merged DataFrame to CSV
output_path = "proteomics_gdsc1_depmap_lof_merge_all.csv"
merged_df.to_csv(os.path.expanduser(output_path), index=False)
print(f"\nMerged data written to: {output_path}")