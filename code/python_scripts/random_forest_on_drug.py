import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.ensemble import RandomForestRegressor

def run_random_forest_on_drug(lof_df, prism_df, merged_df, drug_name):
    print(f"Running Random Forest on {drug_name}...")
    cell_lines = prism_df['cell_line'].tolist()
    available_cell_lines = list(set(cell_lines) & set(lof_df['cell_line'].tolist()))

    lof_df = lof_df[lof_df['cell_line'].isin(available_cell_lines)]
    prism_df = prism_df[prism_df['cell_line'].isin(available_cell_lines)]

    prism_df_target_drug = merged_df[['cell_line', drug_name]]

    # Merge the dataframes
    merged_drug_with_lof = pd.merge(lof_df, prism_df_target_drug, on='cell_line', how='inner')
    merged_drug_with_lof = merged_drug_with_lof.dropna(subset=[drug_name])

    # Separate features and target
    #X = merged_df.drop(columns=['DepMap_ID', TARGET_DRUG])
    y = merged_drug_with_lof[drug_name]
    x = merged_drug_with_lof.drop(columns=['cell_line', drug_name])




    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=42)

    random_forest_model = RandomForestRegressor(n_estimators=100, verbose=1, n_jobs=-1, random_state=42)
    random_forest_model.fit(X_train, y_train)
    return random_forest_model, X_train, X_test, y_train, y_test