import csv
import os

import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path

from project_code.python_scripts.urls import download_google_drive_url

# Constants
GROUPS_MIN_SIZE = 2

# File paths
GDSC_DRUG_IC50_PATH = download_google_drive_url("https://drive.google.com/file/d/1j8aJ4w07coHaMlos3X46V_apYkZXeZUZ")
PRISM_DRUG_SENSITIVITY_PATH = download_google_drive_url("https://drive.google.com/file/d/1DcvKjdJlKE6zy1vHdw-8iHPn3cxF_QvY")

DEPMAP_NUM_LOF_CELL_LINES_PER_GENE = download_google_drive_url("https://drive.google.com/file/d/1Q9NJtvKZbJ4DbfpKuXRWmtdowclgmNfU")
DEPMAP_LOF_PATH = download_google_drive_url("https://drive.google.com/file/d/1iedYFEZoDZxIrBZys8LXAmzD79DMIKsA")
OUTPUT_CSV_DEPMAP_LOF_HOMOZYGOUS = "ttest_homozygous_prism.csv"
OUTPUT_CSV_DEPMAP_LOF_HOMOZYGOUS_GDSC1 = "ttest_homozygous_GDSC1.csv"
DRUG_ROWS_COLUMN_NAME = "string_field_0"



def load_and_prepare_data(drug_path):
    """Load and prepare the input data files."""
    # Load the data
    lof_df = pd.read_csv(DEPMAP_LOF_PATH)
    drug_sensitivity_df = pd.read_csv(
        drug_path)

    # Get filtered genes
    filtered_genes_df = pd.read_csv(DEPMAP_NUM_LOF_CELL_LINES_PER_GENE)
    filtered_genes = filtered_genes_df[filtered_genes_df['count_2'] >= 2]['gene_name'].tolist()

    # Format gene names
    # formatted_genes = [gene.replace(" ", "").replace("(", "..").replace(")", ".")
    #                    for gene in filtered_genes]

    return lof_df, drug_sensitivity_df, filtered_genes


def perform_ttest(lof_df, drug_sensitivity_df, filtered_genes):
    """
    Perform t-test analysis between normal and homozygous cell lines for drug sensitivity.

    Parameters:
    -----------
    lof_df : pandas.DataFrame
        DataFrame containing loss of function data
    drug_sensitivity_df : pandas.DataFrame
        DataFrame containing drug sensitivity data
    filtered_genes : list
        List of genes to analyze

    Returns:
    --------
    pandas.DataFrame
        Results of the t-test analysis
    """
    # Filter genes that exist in lof_df
    filtered_genes_in_lof_df = [gene for gene in filtered_genes if gene in lof_df.columns]

    if not filtered_genes_in_lof_df:
        raise ValueError("No valid genes found in lof_df. Please check the gene names.")

    # Initialize results list
    results = []

    # Set index for faster lookups
    lof_df.set_index('Unnamed: 0', inplace=True)
    drug_sensitivity_df.set_index(DRUG_ROWS_COLUMN_NAME, inplace=True)

    common_cells = set(lof_df.index).intersection(set(drug_sensitivity_df.index))

    print(f"Number of cell lines in lof_df: {len(lof_df.index)}")
    print(f"Number of cell lines in drug_sensitivity_df: {len(drug_sensitivity_df.index)}")
    print(f"Number of common cell lines: {len(common_cells)}")

    # Filter both dataframes to include only common cell lines
    lof_df = lof_df.loc[list(common_cells)]
    drug_sensitivity_df = drug_sensitivity_df.loc[list(common_cells)]

    # Get drug names (all columns except index)
    drug_names = drug_sensitivity_df.columns

    i = 0
    for gene in filtered_genes_in_lof_df:
        i+=1
        print(f"Processing gene {gene} number {i}")
        # Get cell lines for normal and homozygous groups
        normal_cells = lof_df[lof_df[gene] == 0].index
        homozygous_cells = lof_df[lof_df[gene] == 2].index

        for drug in drug_names:
            # Get drug sensitivity values
            normal_sens = drug_sensitivity_df.loc[normal_cells, drug].dropna()
            homozygous_sens = drug_sensitivity_df.loc[homozygous_cells, drug].dropna()

            # Perform t-test if both groups have at least two values
            if len(normal_sens) > 1 and len(homozygous_sens) > 1:
                # Perform t-test
                ttest_result = stats.ttest_ind(normal_sens, homozygous_sens)

                # Calculate statistics
                result = {
                    'Gene': gene,
                    'Drug': drug,
                    'P_Value': ttest_result.pvalue,
                    'Mean_Normal': normal_sens.mean(),
                    'Count_Normal': len(normal_sens),
                    'Std_Normal': normal_sens.std(),
                    'Median_Normal': normal_sens.median(),
                    'Mean_Homozygous': homozygous_sens.mean(),
                    'Count_Homozygous': len(homozygous_sens),
                    'Std_Homozygous': homozygous_sens.std(),
                    'Median_Homozygous': homozygous_sens.median(),
                    'Effect_size_mean': normal_sens.mean() - homozygous_sens.mean(),
                    'Effect_size_median': normal_sens.median() - homozygous_sens.median()
                }
                results.append(result)

    # Convert results to DataFrame
    results_df = pd.DataFrame(results)

    # Reorder columns to match R output
    column_order = ['Gene', 'Drug', 'P_Value', 'Mean_Normal', 'Count_Normal',
                    'Std_Normal', 'Median_Normal', 'Mean_Homozygous',
                    'Count_Homozygous', 'Std_Homozygous', 'Median_Homozygous',
                    'Effect_size_mean', 'Effect_size_median']
    results_df = results_df[column_order]

    return results_df


def write_results_to_csv(results_df, output_path):
    """
    Write the results DataFrame to a CSV file with error handling.

    Parameters:
    -----------
    results_df : pandas.DataFrame
        The DataFrame containing the t-test results
    output_path : str
        The file path where the CSV should be saved
    """
    try:
        # Check if results_df is empty
        if results_df.empty:
            print("Warning: No results to write to CSV")
            return False

        # Create directory if it doesn't exist
        output_dir = os.path.dirname(output_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Write to CSV
        results_df.to_csv(output_path, index=False)
        print(f"Results successfully written to: {output_path}")
        return True

    except PermissionError:
        print(f"Error: Permission denied when writing to {output_path}")
        return False
    except Exception as e:
        print(f"Error writing results to CSV: {str(e)}")
        return False


def main():
    """Main function to run the analysis."""
    # Load and prepare data
    lof_df, drug_sensitivity_df, filtered_genes = load_and_prepare_data(GDSC_DRUG_IC50_PATH)
    # Method 1: Using iloc to get first 5 rows and columns

    # Perform t-test analysis
    ttest_results = perform_ttest(lof_df, drug_sensitivity_df, filtered_genes)

    # Display results
    print(ttest_results)
    write_results_to_csv(ttest_results, OUTPUT_CSV_DEPMAP_LOF_HOMOZYGOUS_GDSC1)

    # Optionally save results
    # ttest_results.to_csv(OUTPUT_CSV_DEPMAP_LOF_HOMOZYGOUS, index=False)


if __name__ == "__main__":
    main()
