from sqlite3.dbapi2 import paramstyle

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import make_scorer, mean_squared_error
from scipy.stats import randint, uniform


def run_random_forest_classifier_on_drug(lof_df, prism_df, drug_name):
    print(f"Running Random Forest on {drug_name}...")

    # Separate the 'cell_line' column
    cell_line_col = prism_df[['cell_line']]

    # Convert the rest to strings
    df_strings = prism_df.drop(columns='cell_line').astype(str)

    # Combine back together
    prism_df_str = pd.concat([cell_line_col, df_strings], axis=1)

    cell_lines = prism_df_str['cell_line'].tolist()
    available_cell_lines = list(set(cell_lines) & set(lof_df['cell_line'].tolist()))

    lof_df = lof_df[lof_df['cell_line'].isin(available_cell_lines)]
    prism_df_str = prism_df_str[prism_df['cell_line'].isin(available_cell_lines)]

    prism_df_target_drug = prism_df_str[['cell_line', drug_name]]

    # Merge the dataframes
    merged_drug_with_lof = pd.merge(lof_df, prism_df_target_drug, on='cell_line', how='inner')
    merged_drug_with_lof = merged_drug_with_lof.dropna(subset=[drug_name])

    # Separate features and target
    #X = merged_df.drop(columns=['DepMap_ID', TARGET_DRUG])
    y = merged_drug_with_lof[drug_name]
    x = merged_drug_with_lof.drop(columns=['cell_line', drug_name])


    # Create base model
    #rf_base = RandomForestClassifier(random_state=42, n_jobs=-1, bootstrap=True, max_depth=10, max_features=None, min_samples_leaf=4, min_samples_split=8, n_estimators=100, verbose=1)
    rf_base = RandomForestClassifier(random_state=42, n_jobs=-1)


    # Best params
    # Best parameters: {'bootstrap': True, 'max_depth': 10, 'max_features': None, 'min_samples_leaf': 4, 'min_samples_split': 8, 'n_estimators': 100}







    # Split the data into training and testing sets
    # First split: separate test set (20% of data)
    # X_temp, X_test, y_temp, y_test = train_test_split(x, y, test_size=0.1, random_state=42)

    # Second split: divide remaining data into train and validation (80% train, 20% validation of remaining data)
    # This results in a 64-16-20 split of the original data
    # X_train, X_val, y_train, y_val = train_test_split(X_temp, y_temp, test_size=0.2, random_state=42)
    X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=42)

    print(f"Training set size: {X_train.shape[0]} samples")
    print(f"Validation set size: {X_test.shape[0]} samples")
    print(f"Test set size: {X_test.shape[0]} samples")
    #random_search.fit(X_train, y_train)
    #best_rf_model = random_search.best_estimator_

    rf_base.fit(X_train, y_train)

    return rf_base, X_train, X_test, X_test, y_train, y_test, y_test, x, y