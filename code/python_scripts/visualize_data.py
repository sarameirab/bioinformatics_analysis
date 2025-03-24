from urls import *
import pandas as pd
from utils import *

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np
from urls import *
from utils import *


def plot_protein_drug_correlation(protein_id, drug_id):
    """
    Creates a scatter plot with correlation line for a specific protein-drug pair using DepMap IDs
    """
    # Read proteomics data and model mapping
    proteomics_path = download_url(PROTEOMICS_MATRIX_AVERAGE_ZSCORE)
    proteomics_df = pd.read_csv(proteomics_path, sep='\t')
    
    model_mapping = pd.read_csv("~/workspace/ai_agent_lab_data/depmap/Model.csv")
    
    # Map Sanger IDs to DepMap IDs in proteomics data
    proteomics_with_depmap = proteomics_df.merge(
        model_mapping[['SangerModelID', 'ModelID']].rename(columns={'ModelID': 'depmap_id'}),
        left_on='model_id',
        right_on='SangerModelID',
        how='left'
    )
    
    # Read drug sensitivity data
    gdsc_path = download_url(GDSC1_DOSE_RESPONSE)
    gdsc_df = pd.read_csv(gdsc_path)
    
    # Filter for specific drug and GDSC1 dataset
    gdsc1_filtered = gdsc_df[gdsc_df['DATASET'] == 'GDSC1']
    drug_data = gdsc1_filtered[gdsc1_filtered['DRUG_ID'] == drug_id]
    
    if len(drug_data) == 0:
        raise ValueError(f"No data found for drug_id {drug_id}")
    
    # Map ARXSPAN_IDs to DepMap IDs in drug data
    drug_data = drug_data.merge(
        model_mapping[['ARXSPANModelID', 'ModelID']].rename(columns={'ModelID': 'depmap_id'}),
        left_on='ARXSPAN_ID',
        right_on='ARXSPANModelID',
        how='inner'
    )
    
    # Get protein expression for specific protein
    protein_data = proteomics_with_depmap[['depmap_id', protein_id]].dropna()
    
    # Merge protein expression with drug sensitivity data
    merged_data = protein_data.merge(
        drug_data[['depmap_id', 'Z_SCORE_PUBLISHED']],
        on='depmap_id',
        how='inner'
    ).rename(columns={
        protein_id: 'protein_expression',
        'Z_SCORE_PUBLISHED': 'drug_sensitivity'
    })
    
    # Create the plot
    plt.figure(figsize=(10, 8))
    sns.set_style("whitegrid")
    
    # Scatter plot
    sns.scatterplot(
        data=merged_data,
        x='protein_expression',
        y='drug_sensitivity',
        alpha=0.6
    )
    
    # Add regression line and statistics only if enough data points
    if len(merged_data) > 2:
        # Regression line
        sns.regplot(
            data=merged_data,
            x='protein_expression',
            y='drug_sensitivity',
            scatter=False,
            color='red'
        )
        
        # Calculate correlation
        correlation, p_value = stats.pearsonr(
            merged_data['protein_expression'],
            merged_data['drug_sensitivity']
        )
        
        # Add statistics text
        stats_text = (f"Pearson r = {correlation:.3f}\n"
                     f"p-value = {p_value:.2e}\n"
                     f"n = {len(merged_data)}")
        
        plt.text(0.05, 0.95, stats_text,
                transform=plt.gca().transAxes,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Add title and labels
    plt.title(f"Correlation between {protein_id} and Drug {drug_id}", pad=20)
    plt.xlabel(f"{protein_id} Expression (z-score)")
    plt.ylabel("Drug Sensitivity (z-score)")
    
    # Create plots directory if it doesn't exist
    plots_dir = os.path.expanduser("~/workspace/bioinformatics_analysis/plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    # Save and show plot
    output_path = os.path.join(plots_dir, f"{protein_id}_{drug_id}_correlation.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {output_path}")
    
    plt.show()

# Example usage
if __name__ == "__main__":
    plot_protein_drug_correlation("O00750", 1498)