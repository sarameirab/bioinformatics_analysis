from project_code.python_scripts.urls import download_google_drive_url, DEPMAP_MUTATION_DATA, PRISM_DRUG_REPURPOSING, \
    DEPMAP_MODEL_MAPPING, download_url, DEPMAP_EXPRESSION
from random_forest_on_drug import run_random_forest_on_multiple_drugs_and_write_results_to_csv
import pandas as pd
lof_path = download_google_drive_url(DEPMAP_MUTATION_DATA)
prism_path = download_google_drive_url(PRISM_DRUG_REPURPOSING)
model_path = download_google_drive_url(DEPMAP_MODEL_MAPPING)

print("Loading data...")
lof_df = pd.read_csv(lof_path)
prism_df = pd.read_csv(prism_path)
model_df = pd.read_csv(model_path)
lof_df = lof_df.rename(columns={'Unnamed: 0': 'cell_line'})
prism_df = prism_df.rename(columns={'Unnamed: 0': 'cell_line'})
model_df = model_df.rename(columns={'ModelID': 'cell_line'})
cell_lines = prism_df['cell_line'].tolist()
available_cell_lines = list(set(cell_lines) & set(lof_df['cell_line'].tolist()))
len(available_cell_lines)
lof_df = lof_df[lof_df['cell_line'].isin(available_cell_lines)]
prism_df = prism_df[prism_df['cell_line'].isin(available_cell_lines)]
model_df = model_df[model_df['cell_line'].isin(available_cell_lines)]
columns_to_keep = ['cell_line', 'DepmapModelType', 'OncotreePrimaryDisease', 'OncotreeSubtype', 'OncotreeCode', 'Sex', 'PrimaryOrMetastasis', 'OncotreeLineage', 'GrowthPattern', 'SampleCollectionSite']
model_df = model_df[columns_to_keep]

expression_path =download_google_drive_url(DEPMAP_EXPRESSION)
expression_df = pd.read_csv(expression_path)
expression_df = expression_df.rename(columns={'Unnamed: 0': 'ProfileID'})

expression_mapping_path = download_url("https://42basepairs.com/download/s3/d3b-openaccess-us-east-1-prd-pbta/pbta-splicing/v8/OmicsDefaultModelProfiles.csv")
expression_mapping_df = pd.read_csv(expression_mapping_path)
expression_mapping_df = expression_mapping_df.rename(columns={'Unnamed: 0': 'ProfileID'})
expression_merge = pd.merge(expression_mapping_df, expression_df,
                          left_on="ProfileID",
                          right_on="ProfileID",
                          how='inner')
expression_merge = expression_merge.rename(columns={'ModelID': 'cell_line'})
expression_merge = expression_merge.drop(['ProfileType', 'ProfileID'], axis=1)


feature_df = pd.merge(model_df, lof_df, on='cell_line', how='inner')
feature_df = pd.merge(feature_df, expression_merge, on='cell_line', how='inner')

merged_df = pd.merge(feature_df, prism_df, on='cell_line', how='inner')

import pandas as pd
random_forest_results_path = download_google_drive_url("https://drive.google.com/file/d/1n8hV1mfLbsqR-zABud4JX1LccxhVhRH9")
random_forest_results = pd.read_csv(random_forest_results_path)

random_forest_results_sorted = random_forest_results.sort_values(by="pearson_correlation", ascending=False)
top_5_drugs = random_forest_results_sorted.head(5)

print("Running random forest on " +str(top_5_drugs))

run_random_forest_on_multiple_drugs_and_write_results_to_csv(feature_df, prism_df, merged_df, top_5_drugs["drug"], "random_forest_results_with_more_data_2.csv")