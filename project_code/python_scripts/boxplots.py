import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

from urls import download_google_drive_url, DEPMAP_MUTATION_DATA, PRISM_DRUG_REPURPOSING


def generate_boxplot(gene_name, drug_name):
    """
    Generates a boxplot visualizing drug sensitivity by LOF status of a given gene.

    Args:
        gene_name (str): The name of the gene (e.g., "ABCA1 (19)").
        drug_name (str): The name of the drug (e.g., "VISOMITIN (BRD:BRD-K00003375-004-01-9)").
        data_path (str): The path to the directory containing the Depmap data.
    """

    # Define file paths
    lof_file = download_google_drive_url(DEPMAP_MUTATION_DATA)
    drug_file = download_google_drive_url(PRISM_DRUG_REPURPOSING)

    # Load data
    try:
        lof_data = pd.read_csv(lof_file, index_col=0)
        drug_data = pd.read_csv(drug_file, index_col=0)
    except FileNotFoundError as e:
        print(f"Error loading data: {e}. Please ensure the data files are in the correct path.")
        return

    # Ensure the gene and drug exist in the data
    if gene_name not in lof_data.columns:
        raise ValueError(f"Gene '{gene_name}' not found in LOF data.")
    if drug_name not in drug_data.columns:
        raise ValueError(f"Drug '{drug_name}' not found in drug sensitivity data.")

    # Merge data based on cell line IDs
    common_cell_lines = list(set(lof_data.index).intersection(drug_data.index))
    merged_data = pd.DataFrame({
        'LOF': lof_data.loc[common_cell_lines, gene_name],
        'Sensitivity': drug_data.loc[common_cell_lines, drug_name]
    })

    # Filter out rows with NA in Sensitivity
    merged_data = merged_data.dropna(subset=['Sensitivity'])

    # Classify cell lines based on LOF status
    merged_data['LOF_Status'] = merged_data['LOF'].apply(lambda x: 'LOF' if x == 2 else ('No LOF' if x == 0 else None))

    # Compute statistics
    stats_df = merged_data.groupby('LOF_Status')['Sensitivity'].agg(
        Count='size',
        Mean='mean',
        SD='std'
    ).reset_index()

    # Perform t-test
    lof_sensitivity = merged_data[merged_data['LOF_Status'] == 'LOF']['Sensitivity']
    no_lof_sensitivity = merged_data[merged_data['LOF_Status'] == 'No LOF']['Sensitivity']

    # Check if both groups have data for t-test
    if len(lof_sensitivity) > 1 and len(no_lof_sensitivity) > 1:
        t_test_result = stats.ttest_ind(lof_sensitivity, no_lof_sensitivity)  # Welch's t-test
        p_value = t_test_result.pvalue
    else:
        p_value = float('nan')  # Not enough data for t-test

    # Create the plot
    plt.figure(figsize=(10, 7))
    sns.boxplot(x='LOF_Status', y='Sensitivity', data=merged_data, palette={'LOF': 'skyblue', 'No LOF': 'lightcoral'},
                flierprops={"marker": "None"})
    sns.stripplot(x='LOF_Status', y='Sensitivity', data=merged_data, color='black', size=4, jitter=0.2)

    plt.title(f"Sensitivity to {drug_name} by LOF Status of {gene_name}", fontsize=14)
    plt.xlabel("LOF Status", fontsize=12)
    plt.ylabel("Drug Sensitivity", fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    # Add subtitle with statistics
    lof_count = stats_df[stats_df['LOF_Status'] == 'LOF']['Count'].values[0] if 'LOF' in stats_df[
        'LOF_Status'].values else 0
    no_lof_count = stats_df[stats_df['LOF_Status'] == 'No LOF']['Count'].values[0] if 'No LOF' in stats_df[
        'LOF_Status'].values else 0

    subtitle_text = (f"LOF Cell Lines: {lof_count} | No LOF Cell Lines: {no_lof_count} "
                     f"| p-value: {p_value:.3g}" if not pd.isna(p_value) else f"| p-value: N/A (insufficient data)")
    plt.suptitle(subtitle_text, fontsize=10, y=0.95)

    # Print the statistics
    print("--- Statistics ---")
    print(stats_df)
    print(f"p-value: {p_value:.3g}" if not pd.isna(p_value) else "p-value: N/A (insufficient data)")

    plt.show()


# Example usage
gene_name = "PARD3 (56288)"
drug_name = "BICALUTAMIDE (BRD:BRD-A29485665-001-12-8)"
generate_boxplot(gene_name, drug_name)