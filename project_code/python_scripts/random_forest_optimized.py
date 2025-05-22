
from urls import download_google_drive_url, PRISM_DRUG_REPURPOSING, DEPMAP_MUTATION_DATA
import pandas as pd
import numpy as np
from random_forest_on_drug import run_random_forest_on_drug
from sklearn.metrics import root_mean_squared_error
from sklearn.metrics import root_mean_squared_error

from machine_learning_utils import create_plots_for_model


DRUG_1 = "MILADEMETAN (BRD:BRD-K00003406-001-01-9)"
DRUG_2 = "IDASANUTLIN (BRD:BRD-K62627508-001-01-5)"
DRUG_3 = "SAR405838 (BRD:BRD-A16035238-001-01-7)"
DRUG_4 = "CGM097 (BRD:BRD-K79584249-001-01-3)"
DRUG_5 = "AMG-232 (BRD:BRD-K64925568-001-01-8)"
DRUG_6 = "RG7112 (BRD:BRD-A78210457-001-01-5)"
TARGET_DRUG = DRUG_1

lof_path = download_google_drive_url(DEPMAP_MUTATION_DATA)
prism_path = download_google_drive_url(PRISM_DRUG_REPURPOSING)

lof_df = pd.read_csv(lof_path)
prism_df = pd.read_csv(prism_path)

lof_df = lof_df.rename(columns={'Unnamed: 0': 'cell_line'})
prism_df = prism_df.rename(columns={'Unnamed: 0': 'cell_line'})


cell_lines = prism_df['cell_line'].tolist()
available_cell_lines = list(set(cell_lines) & set(lof_df['cell_line'].tolist()))

#%%
lof_df = lof_df[lof_df['cell_line'].isin(available_cell_lines)]
prism_df = prism_df[prism_df['cell_line'].isin(available_cell_lines)]



# Merge the dataframes
merged_df = pd.merge(lof_df, prism_df, on='cell_line', how='inner')

random_forest_model, X_train, X_test, X_val, y_train, y_test, y_val, x, y = run_random_forest_on_drug(lof_df, prism_df, merged_df, TARGET_DRUG)
random_forest_predictions_test = random_forest_model.predict(X_val)
forest_rmse_test = root_mean_squared_error(y_val, random_forest_predictions_test)

random_forest_predictions_train = random_forest_model.predict(X_train)
forest_rmse_train = root_mean_squared_error(y_train, random_forest_predictions_train)


baseline_prediction = np.full_like(y_val, y_train.mean())
baseline_rmse = root_mean_squared_error(y_val, baseline_prediction)

print(f"Random Forest RMSE test: {forest_rmse_test}")
print(f"Random Forest RMSE training: {forest_rmse_train}")
print(f"Baseline RMSE: {baseline_rmse}")



random_forest_predictions_training = random_forest_model.predict(X_train)
forest_rmse_training = root_mean_squared_error(y_train, random_forest_predictions_training)
random_forest_predictions_val = random_forest_model.predict(X_val)
forest_rmse_val = root_mean_squared_error(y_val, random_forest_predictions_val)

create_plots_for_model(y_val, y_train, random_forest_predictions_val, random_forest_predictions_training, y)


# feature_importances = zip(X_train.columns, random_forest_model.feature_importances_)
#
# # Sort the features by their importance scores in descending order
# sorted_features = sorted(feature_importances, key=lambda x: x[1], reverse=True)
#
# # Print the sorted features
# for name, score in sorted_features:
#     print(f"{name}: {score:.4f}")
