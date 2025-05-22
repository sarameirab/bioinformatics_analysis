
import os
import pandas as pd
import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import numpy as np
from sklearn.model_selection import KFold

# from project_code.python_scripts.urls import download_url, DEPMAP_MUTATION_DATA, PRISM_DRUG_REPURPOSING

DRUG_1 = "MILADEMETAN _BRD:BRD-K00003406-001-01-9_"
DRUG_2 = "IDASANUTLIN _BRD:BRD-K62627508-001-01-5_"
DRUG_3 = "SAR405838 _BRD:BRD-A16035238-001-01-7_"
DRUG_4 = "CGM097 _BRD:BRD-K79584249-001-01-3_"
DRUG_5 = "AMG-232 _BRD:BRD-K64925568-001-01-8_"
DRUG_6 = "RG7112 _BRD:BRD-A78210457-001-01-5_"

# Set your Google Cloud project ID
PROJECT_ID = "cobalt-logic-447008-p7"

# Initialize BigQuery client

# Define BigQuery table names
LOF_TABLE = "depmap.OmicsSomaticMutationsMatrixDamaging"
PRISM_TABLE = "depmap.PRISM_Repurposing_Public_24Q2_subsetted"

# lof_path = download_url(DEPMAP_MUTATION_DATA)
# prism_path =download_url(PRISM_DRUG_REPURPOSING)
#
# lof_df = pd.read_csv(lof_path)
# prism_df = pd.read_csv(prism_path)

# # Define the target drug
# TARGET_DRUG = DRUG_1
#
# query_cell_lines = f"""
#     SELECT string_field_0 FROM `{PRISM_TABLE}`
# """
# cell_lines_df = client.query(query_cell_lines).to_dataframe()
# cell_lines = cell_lines_df['string_field_0'].tolist()
#
# # Fetch LOF data
# query_lof = f"""
#     SELECT * FROM `{LOF_TABLE}`
# """
# lof_df = client.query(query_lof).to_dataframe()
# lof_df = lof_df.set_index('sample_id')
#
# # Get the intersection of cell lines from both tables
# available_cell_lines = list(set(cell_lines) & set(lof_df.columns))
#
# # Filter LOF data to only include cell lines present in both tables
# lof_df = lof_df[available_cell_lines]
# lof_df = lof_df.transpose()
# lof_df.index.name = 'DepMap_ID'
# lof_df = lof_df.reset_index()
#
# # Fetch drug sensitivity data
# query_drug = f"""
#     SELECT string_field_0, `{TARGET_DRUG}` FROM `{PRISM_TABLE}`
# """
# drug_df = client.query(query_drug).to_dataframe()
# drug_df = drug_df.rename(columns={'string_field_0': 'DepMap_ID'})
#
# # Merge the dataframes
# merged_df = pd.merge(lof_df, drug_df, on='DepMap_ID', how='inner')
# merged_df = merged_df.dropna(subset=[TARGET_DRUG])
#
# # Separate features and target
# X = merged_df.drop(columns=['DepMap_ID', TARGET_DRUG])
# y = merged_df[TARGET_DRUG]
# scaler = StandardScaler()
# X_scaled = scaler.fit_transform(X)
#
# # Split data into training and testing sets
# X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)
# input_shape = X_train.shape[1]
# model = keras.Sequential([
#         layers.Input(shape=(input_shape,)),
#         layers.BatchNormalization(),
#         layers.Dense(256, activation='relu',
#                      kernel_regularizer=tf.keras.regularizers.l1_l2(l1=0.001, l2=0.001)),
#         layers.Dropout(0.4),
#         layers.Dense(128, activation='relu',
#                      kernel_regularizer=tf.keras.regularizers.l1_l2(l1=0.001, l2=0.001)),
#         layers.Dropout(0.3),
#         layers.Dense(64, activation='relu'),
#         layers.Dense(1, activation='linear')
#     ])
#
# model.compile(
#     optimizer=tf.keras.optimizers.Adam(learning_rate=0.0001),
#     loss='mse',
#     metrics=['mae', tf.keras.metrics.R2Score()]
# )
# early_stopping = tf.keras.callbacks.EarlyStopping(
#     monitor='val_loss',
#     patience=10,
#     restore_best_weights=True
# )
#
# model_checkpoint = tf.keras.callbacks.ModelCheckpoint(
#     'my_model.keras',
#     monitor='val_loss',
#     save_best_only=True
# )
#
#
#
# # K-Fold Cross-Validation
# def cross_validate_model(X, y, input_shape, n_splits=5):
#     # Convert to numpy arrays to ensure indexing works
#     X = np.array(X)
#     y = np.array(y)
#
#     kfold = KFold(n_splits=n_splits, shuffle=True, random_state=42)
#     cv_scores = []
#
#     for train_index, val_index in kfold.split(X):
#         X_train, X_val = X[train_index], X[val_index]
#         y_train, y_val = y[train_index], y[val_index]
#
#         history = model.fit(
#             X_train, y_train,
#             validation_data=(X_val, y_val),
#             epochs=100,
#             batch_size=32,
#             callbacks=[early_stopping, model_checkpoint],
#             verbose=1
#         )
#
#         val_r2 = model.evaluate(X_val, y_val)[2]  # R2 score index
#         cv_scores.append(val_r2)
#
#     return np.mean(cv_scores), np.std(cv_scores)
#
#
# mean_cv_score, std_cv_score = cross_validate_model(X_train, y_train, input_shape)
# print(f"Cross-Validation R2: {mean_cv_score:.4f} ± {std_cv_score:.4f}")
#
#
# X_test_np = np.array(X_test)
# y_test_np = np.array(y_test)
# model.evaluate(X_test_np, y_test_np)
#
# predictions = []
# for current in X_test:
#     reshape_current = current.reshape(1, -1)
#     prediction = model.predict(reshape_current)
#     predictions.append(prediction)
