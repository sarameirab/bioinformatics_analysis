import matplotlib.pyplot as plt
from scipy.stats import pearsonr

def create_plots_for_model(y_test, y_train, predictions, predictions_on_training, y):
    """
    Create plots for model evaluation.
    
    Parameters:
    - y_test: True values for the test set.
    - y_train: True values for the training set.
    - predictions: Model predictions for the test set.
    - predictions_on_training: Model predictions for the training set.
    - y: All sensitivities.
    
    Returns:
    None
    """
    
    # Convert predictions to lists
    true_values_list = y_test.tolist()
    true_values_list = y_test.tolist()
    true_values_on_training_list = y_train.tolist()
    all_sensitivities = y.tolist()

    r, p_value = pearsonr(true_values_list, predictions)
    plt.title(f'Scatter Plot with R={r:.2f}')
    plt.scatter(true_values_list, predictions)
    plt.show()


    plt.hist(true_values_list, bins=5, color='skyblue', edgecolor='black')
    plt.title("True values distribution")
    plt.show()

    plt.hist(predictions, bins=5, color='skyblue', edgecolor='black')
    plt.title("Prediction distribution")
    plt.show()

    plt.hist(true_values_on_training_list, bins=5, color='skyblue', edgecolor='black')
    plt.title("True values on training distribution")
    plt.show()


    plt.hist(predictions_on_training, bins=5, color='skyblue', edgecolor='black')
    plt.title("Prediction on training distribution")
    plt.show()


    plt.hist(all_sensitivities, bins=5, color='skyblue', edgecolor='black')
    plt.title("True All sensitivities distribution")
    plt.show()
