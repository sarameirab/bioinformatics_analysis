import pandas as pd
import numpy as np
from scipy import stats
import math

import matplotlib.pyplot as plt
import seaborn as sns

from urls import download_google_drive_url, PROTEOMICS_MATRIX_AVERAGE_ZSCORE, DEPMAP_MODEL_MAPPING

# Read data
proteomics_z_score_path = download_google_drive_url(PROTEOMICS_MATRIX_AVERAGE_ZSCORE)
proteomics_z_score = pd.read_csv(proteomics_z_score_path,
                                sep='\t')

model_mapping_path = download_google_drive_url(DEPMAP_MODEL_MAPPING)
model_mapping = pd.read_csv(model_mapping_path)

# Join proteomics with model mapping
proteomics_z_score_with_depmap_id = proteomics_z_score.merge(
    model_mapping[['SangerModelID', 'ModelID']].rename(columns={'ModelID': 'depmap_id'}),
    left_on='model_id',
    right_on='SangerModelID',
    how='left'
)
proteomics_z_score_with_depmap_id = proteomics_z_score_with_depmap_id.loc[:, ['depmap_id'] + [col for col in proteomics_z_score_with_depmap_id.columns if col != 'depmap_id']]

# Calculate intersection and NA count
intersection = set(proteomics_z_score['model_id']).intersection(set(model_mapping['SangerModelID']))
na_count = proteomics_z_score['model_id'].isna().sum()

# Read GDSC data
gdsc_from_depmap = pd.read_csv("~/workspace/ai_agent_lab_data/GDSC/sanger-dose-response.csv")
gdsc_from_depmap['log2_ic50'] = np.log2(gdsc_from_depmap['IC50_PUBLISHED'])

# Filter GDSC1 data
gdsc1_pre_filtered = gdsc_from_depmap[gdsc_from_depmap['DATASET'] == 'GDSC1']
gdsc1 = gdsc1_pre_filtered[gdsc1_pre_filtered['ARXSPAN_ID'].notna() & (gdsc1_pre_filtered['ARXSPAN_ID'] != '')]

# Pivot drug sensitivity data
gdsc1_drug_sensitivity_for_proteomics = gdsc1.pivot(
    index='DRUG_ID',
    columns='ARXSPAN_ID',
    values='Z_SCORE_PUBLISHED'
).reset_index()


def calculate_protein_drug_correlation_2(proteomics_df, drug_sensitivity_df, protein_id=None, drug_id=None):
    # Make a deep copy to avoid SettingWithCopyWarning
    proteomics_df = proteomics_df.copy()
    drug_sensitivity_df = drug_sensitivity_df.copy()

    # Remove first two rows
    proteomics_df = proteomics_df.iloc[2:].reset_index(drop=True)

    # Convert proteomics values to numeric
    for col in proteomics_df.columns[3:]:
        proteomics_df.loc[:, col] = pd.to_numeric(proteomics_df[col], errors='coerce')

    # Get common depmap IDs
    common_depmap_ids = set(proteomics_df['depmap_id']).intersection(drug_sensitivity_df.columns[1:])
    proteomics_df = proteomics_df[proteomics_df['depmap_id'].isin(common_depmap_ids)].reset_index(drop=True)
    drug_sensitivity_df = drug_sensitivity_df[['DRUG_ID'] + list(common_depmap_ids)]

    results_list = []

    # Filter proteins to process
    proteins_to_process = [protein_id] if protein_id else proteomics_df.columns[3:]

    # Filter drugs to process
    if drug_id:
        drug_sensitivity_df = drug_sensitivity_df[drug_sensitivity_df['DRUG_ID'] == drug_id]
        if drug_sensitivity_df.empty:
            print(f"Drug ID {drug_id} not found in the dataset.")
            return pd.DataFrame()

    # Iterate through proteins
    i = 0
    for protein in proteins_to_process:
        i+=1
        print(f"Calculating protein {protein} number {i}")
        if protein_id:  # Only print progress for all-pairs analysis
            print(f"Processing protein {protein}")

        # Create protein data with depmap_ids
        protein_data = pd.DataFrame({
            'depmap_id': proteomics_df['depmap_id'],
            'protein_value': proteomics_df[protein].astype(float)
        })

        # Iterate through drugs
        for _, drug_row in drug_sensitivity_df.iterrows():
            drug_curr_id = drug_row['DRUG_ID']
            # Create drug data with depmap_ids
            drug_data = pd.DataFrame({
                'depmap_id': drug_sensitivity_df.columns[1:],
                'drug_value': drug_row.iloc[1:].values
            })

            # Merge protein and drug data based on depmap_id
            merged_data = protein_data.merge(drug_data, on='depmap_id', how='inner')
            merged_data = merged_data.dropna()

            if len(merged_data) > 2:
                try:
                    correlation, p_value = stats.pearsonr(
                        merged_data['protein_value'],
                        merged_data['drug_value']
                    )
                    mean_sensitivity = merged_data['drug_value'].mean()
                except (ValueError, TypeError):
                    correlation = p_value = mean_sensitivity = np.nan
            else:
                correlation = p_value = mean_sensitivity = np.nan

            results_list.append({
                'protein_id': protein,
                'drug_id': drug_curr_id,
                'correlation': correlation,
                'p_value': p_value,
                'cell_line_count': len(merged_data),
                'mean_sensitivity': mean_sensitivity
            })

    return pd.DataFrame(results_list)


def calculate_protein_drug_correlation(proteomics_df, drug_sensitivity_df):
    # Make a deep copy to avoid SettingWithCopyWarning
    proteomics_df = proteomics_df.copy()
    drug_sensitivity_df = drug_sensitivity_df.copy()

    # Remove first two rows
    proteomics_df = proteomics_df.iloc[2:].reset_index(drop=True)

    # Convert proteomics values to numeric
    for col in proteomics_df.columns[3:]:
        proteomics_df.loc[:, col] = pd.to_numeric(proteomics_df[col], errors='coerce')

    # Initialize results list instead of DataFrame
    results_list = []

    # Get common depmap IDs
    common_depmap_ids = set(proteomics_df['depmap_id']).intersection(drug_sensitivity_df.columns[1:])
    proteomics_df = proteomics_df[proteomics_df['depmap_id'].isin(common_depmap_ids)].reset_index(drop=True)
    drug_sensitivity_df = drug_sensitivity_df[['DRUG_ID'] + list(common_depmap_ids)]

    # Convert drug sensitivity values to numeric
    for col in drug_sensitivity_df.columns[1:]:
        drug_sensitivity_df[col] = pd.to_numeric(drug_sensitivity_df[col], errors='coerce')

    # Calculate correlations
    i = 0
    for protein in proteomics_df.columns[3:]:
        i += 1
        print(f"Calculating protein {protein} number {i}")
        protein_values = proteomics_df[protein].astype(float)

        for drug_row in range(len(drug_sensitivity_df)):
            drug_id = drug_sensitivity_df.iloc[drug_row]['DRUG_ID']
            drug_values = drug_sensitivity_df.iloc[drug_row, 1:].astype(float)

            # Create merged data with aligned indices
            temp_df = pd.DataFrame({
                'protein_value': protein_values.values,
                'drug_value': drug_values.values
            })
            merged_data = temp_df.dropna()

            # Calculate correlation and statistics
            if len(merged_data) > 2:
                try:
                    correlation, p_value = stats.pearsonr(
                        merged_data['protein_value'].astype(float),
                        merged_data['drug_value'].astype(float)
                    )
                except (ValueError, TypeError):
                    correlation = p_value = np.nan
            else:
                correlation = p_value = np.nan

            # Calculate mean sensitivity
            try:
                mean_sensitivity = merged_data['drug_value'].astype(float).mean()
            except (ValueError, TypeError):
                mean_sensitivity = np.nan

            # Append results to list
            results_list.append({
                'protein_id': protein,
                'drug_id': drug_id,
                'correlation': correlation,
                'p_value': p_value,
                'cell_line_count': len(merged_data),
                'mean_sensitivity': mean_sensitivity
            })

    # Convert results list to DataFrame at the end
    results = pd.DataFrame(results_list)
    return results




def plot_protein_drug_correlation(proteomics_df, drug_sensitivity_df, protein_name, drug_id):
    # Ensure deep copies of the data
    proteomics_df = proteomics_df.copy()
    drug_sensitivity_df = drug_sensitivity_df.copy()

    # Remove first two rows
    proteomics_df = proteomics_df.iloc[2:].reset_index(drop=True)

    # Convert proteomics values to numeric
    for col in proteomics_df.columns[3:]:
        proteomics_df.loc[:, col] = pd.to_numeric(proteomics_df[col], errors='coerce')

    # Get common depmap IDs
    common_depmap_ids = set(proteomics_df['depmap_id']).intersection(drug_sensitivity_df.columns[1:])
    proteomics_df = proteomics_df[proteomics_df['depmap_id'].isin(common_depmap_ids)].reset_index(drop=True)
    drug_sensitivity_df = drug_sensitivity_df[['DRUG_ID'] + list(common_depmap_ids)]

    # Extract drug values for specific drug_id
    drug_row = drug_sensitivity_df[drug_sensitivity_df['DRUG_ID'] == drug_id]
    if drug_row.empty:
        print(f"Drug ID {drug_id} not found in the dataset.")
        return

    # Create a DataFrame with protein values and depmap_ids
    protein_data = pd.DataFrame({
        'depmap_id': proteomics_df['depmap_id'],
        'protein_value': proteomics_df[protein_name].astype(float)
    })

    # Create a DataFrame with drug values and depmap_ids
    drug_data = pd.DataFrame({
        'depmap_id': drug_row.columns[1:],
        'drug_value': drug_row.iloc[0, 1:].values
    })

    # Merge protein and drug data based on depmap_id
    merged_data = protein_data.merge(drug_data, on='depmap_id', how='inner')
    merged_data = merged_data.dropna()

    # Check if there are enough points to calculate correlation
    if len(merged_data) < 3:
        print("Not enough data points to compute correlation.")
        return

    # Calculate correlation and statistics
    correlation, p_value = stats.pearsonr(merged_data['protein_value'], merged_data['drug_value'])
    cell_line_count = len(merged_data)

    # Plot
    plt.figure(figsize=(7, 5))
    sns.scatterplot(data=merged_data, x='protein_value', y='drug_value')
    plt.xlabel(f"{protein_name} Expression (z-score)")
    plt.ylabel(f"Drug Sensitivity (Z-score)")
    plt.title(f"Correlation: {protein_name} vs Drug {drug_id}")

    # Annotate with correlation stats
    plt.text(
        0.05, 0.9,
        f"r = {correlation:.2f}\np = {p_value:.3e}\nCell lines: {cell_line_count}",
        transform=plt.gca().transAxes,
        fontsize=12, bbox=dict(facecolor='white', alpha=0.5)
    )

    plt.show()


# Calculate correlations
#results = calculate_protein_drug_correlation(
#    proteomics_z_score_with_depmap_id,
#    gdsc1_drug_sensitivity_for_proteomics
#)

# Save results
#results.to_csv("~/workspace/proteomics_correlation.csv", index=False)
#
results = calculate_protein_drug_correlation_2(proteomics_z_score_with_depmap_id, gdsc1_drug_sensitivity_for_proteomics, None, None)
results.to_csv("proteomics_drug_sensitivity_correlation.csv", index=False)
