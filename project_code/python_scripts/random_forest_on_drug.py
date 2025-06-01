from sqlite3.dbapi2 import paramstyle

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix, root_mean_squared_error
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import make_scorer, mean_squared_error
from scipy.stats import randint, uniform
import scipy.stats as stats
import traceback

from sklearn.pipeline import make_pipeline, Pipeline
from sklearn.preprocessing import OneHotEncoder

ONE_HOT_ENCODING_COLUMNS = ['DepmapModelType', 'OncotreePrimaryDisease', 'OncotreeSubtype', 'OncotreeCode', 'Sex', 'PrimaryOrMetastasis', 'OncotreeLineage', 'GrowthPattern', 'SampleCollectionSite']
def run_random_forest_on_drug(feature_df, prism_df, merged_df, drug_name):
    print(f"Running Random Forest on {drug_name}...")

    cell_lines = prism_df['cell_line'].tolist()
    available_cell_lines = list(set(cell_lines) & set(feature_df['cell_line'].tolist()))

    feature_df = feature_df[feature_df['cell_line'].isin(available_cell_lines)]
    prism_df = prism_df[prism_df['cell_line'].isin(available_cell_lines)]

    prism_df_target_drug = merged_df[['cell_line', drug_name]]

    # Merge the dataframes
    merged_drug_with_lof = pd.merge(feature_df, prism_df_target_drug, on='cell_line', how='inner')
    merged_drug_with_lof = merged_drug_with_lof.dropna(subset=[drug_name])

    # Separate features and target
    #X = merged_df.drop(columns=['DepMap_ID', TARGET_DRUG])
    y = merged_drug_with_lof[drug_name]
    x = merged_drug_with_lof.drop(columns=['cell_line', drug_name])

    # Define parameter grid
    param_grid = {
        'n_estimators': [100, 200, 300],
        'max_depth': [None, 10, 20, 30],
        'min_samples_split': [2, 5, 10],
        'min_samples_leaf': [1, 2, 4],
        'max_features': ['auto', 'sqrt', None],
        'bootstrap': [True, False]
    }
    param_distributions = {
        'n_estimators': [100, 200, 300],  # you can keep as list
        'max_depth': [None, 10, 20, 30],  # you can keep as list
        'min_samples_split': randint(2, 11),  # will sample integers from 2 to 10
        'min_samples_leaf': randint(1, 5),  # will sample integers from 1 to 4
        'max_features': ['auto', 'sqrt', None],  # keep as list for categorical
        'bootstrap': [True, False]  # keep as list for boolean
    }

    # Create base model
    #RMSE: 1.0122604602030016. Pearson correaltion 0.64
    #rf_base = RandomForestRegressor(random_state=42, n_jobs=-1, bootstrap=True, max_depth=10, min_samples_leaf=4, min_samples_split=8, n_estimators=100, verbose=1)
    #RMSE: 1.189419455947339. Pearson correlation 0.57

    numeric_cols = [c for c in x.columns if c not in ONE_HOT_ENCODING_COLUMNS]

    preprocess = ColumnTransformer(
        transformers=[
            ('cat', OneHotEncoder(handle_unknown='ignore'), ONE_HOT_ENCODING_COLUMNS),
            ('num', 'passthrough', numeric_cols)
        ]
    )
    rf_base = RandomForestRegressor(random_state=42, n_jobs=-1, n_estimators=100)

    pipeline = Pipeline(steps=[('prep', preprocess),
                               ('rf', rf_base)])

    # Create RMSE scorer (negative because GridSearchCV tries to maximize the score)
    rmse_scorer = make_scorer(lambda y_true, y_pred: -np.sqrt(mean_squared_error(y_true, y_pred)))

    # Initialize RandomizedSearchCV
    random_search = RandomizedSearchCV(
        estimator=rf_base,
        param_distributions=param_distributions,
        n_iter=100,  # number of parameter settings that are sampled
        cv=5,
        scoring=rmse_scorer,
        n_jobs=-1,
        verbose=2,
        random_state=42,
        return_train_score=True
    )

    # Fit GridSearchCV
    # Best params
    # Best parameters: {'bootstrap': True, 'max_depth': 10, 'max_features': None, 'min_samples_leaf': 4, 'min_samples_split': 8, 'n_estimators': 100}







    #Split the data into training and testing sets
    #First split: separate test set (20% of data)
    X_temp, X_test, y_temp, y_test = train_test_split(x, y, test_size=0.1, random_state=42)

    #Second split: divide remaining data into train and validation (80% train, 20% validation of remaining data)
    #This results in a 64-16-20 split of the original data
    X_train, X_val, y_train, y_val = train_test_split(X_temp, y_temp, test_size=0.2, random_state=42)
    #X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=42)

    print(f"Training set size: {X_train.shape[0]} samples")
    print(f"Validation set size: {X_test.shape[0]} samples")
    print(f"Test set size: {X_test.shape[0]} samples")
    #random_search.fit(X_train, y_train)
    #best_rf_model = random_search.best_estimator_

    pipeline.fit(X_train, y_train)

    # print("Best parameters:", random_search.best_params_)
    # print("Best cross-validation score:", random_search.best_score_)

    #
    # random_forest_model = RandomForestRegressor(n_estimators=100, verbose=1, n_jobs=-1, random_state=42)
    # random_forest_model.fit(X_train, y_train)
    return pipeline, X_train, X_test, X_val, y_train, y_test, y_val, x, y

def run_random_forest_on_multiple_drugs_and_write_results_to_csv(lof_df, prism_df, merged_df, drug_names, csv_name="random_forest_results.csv"):

    # Create empty lists to store results
    results = []

    # Iterate through each drug
    for drug in drug_names:
        try:
            # Run random forest for the current drug
            pipeline, X_train, X_test, X_val, y_train, y_test, y_val, x, y= run_random_forest_on_drug(lof_df, prism_df,
                                                                                                 merged_df, drug)

            # Get predictions
            test_predictions = pipeline.predict(X_val)

            # Calculate metrics
            rmse = root_mean_squared_error(y_val, test_predictions)
            baseline_pred = np.full_like(y_val, y_train.mean())
            baseline_rmse = root_mean_squared_error(y_val, baseline_pred)
            pearson_corr, _ = stats.pearsonr(y_val, test_predictions)

            rf = pipeline.named_steps['rf']  # RandomForestRegressor
            preprocessor = pipeline.named_steps['prep']  # ColumnTransformer

            # --- 2.  Recover the feature names *after* one-hot expansion -----------------
            try:
                # scikit-learn ≥ 1.0
                feature_names = preprocessor.get_feature_names_out()
            except AttributeError:
                # scikit-learn < 1.0 – fall back to manual assembly
                feature_names = []
                for name, trans, cols in preprocessor.transformers_:
                    if trans == 'passthrough':
                        feature_names.extend(cols)
                    else:  # e.g. OneHotEncoder
                        fn = trans.get_feature_names_out(cols)
                        feature_names.extend(fn)

            # --- 3.  Pair names with importances and tidy them ---------------------------
            # 2) Get the post-encoding feature names
            feature_names = preprocessor.get_feature_names_out()

            # 3) Pair every name with its importance
            importances = rf.feature_importances_

            # 4) Build a *sorted* dict  →  {feature: importance, ...}
            importance_dict = dict(
                sorted(
                    zip(feature_names, importances),  # (name, value) pairs
                    key=lambda item: item[1],  # sort by importance
                    reverse=True  # biggest first
                )
            )

            # 5) Use it
            print(list(importance_dict.items())[:10])  # top-10 as a sanity check

            # Create a result dictionary with all metrics
            result = {
                'drug': drug,
                'rmse': rmse,
                'baseline_rmse': baseline_rmse,
                'pearson_correlation': pearson_corr,
                **importance_dict  # Unpack feature importances
            }

            results.append(result)
            print(f"Processed drug: {drug}")

        except Exception as e:
            print(f"Error processing drug {drug}: {str(e)} ")
            traceback.print_exc()
            continue

    # Create DataFrame from results
    results_df = pd.DataFrame(results)

    # Save results to CSV
    results_df.to_csv(csv_name, index=False)

    # Display first few rows of the results
    return results_df
