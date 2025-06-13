
# final_models_22features

This repository contains a Python script that trains and evaluates machine learning models using a set of 22 selected genomic features. The workflow includes data loading, feature selection, model training, evaluation, and saving the trained models.

## Script Overview

### Script: `final_models_22features_english.py`

This script is organized into the following sections:

1. **Instructions**
   - Describes the 22 selected features used in the analysis.

2. **Data Loading**
   - Loads training and testing datasets.
   - Features are extracted and labels are defined.

3. **Model Training**
   - Trains several classifiers (e.g., Random Forest, XGBoost) using the selected features.
   - Includes hyperparameter tuning and cross-validation.

4. **Model Evaluation**
   - Calculates performance metrics like accuracy, precision, recall, F1-score, and ROC-AUC.

5. **Model Saving**
   - Saves the trained models and performance results for future use.

6. **Prediction & Result Output**
   - Applies trained models to a separate prediction dataset.
   - Outputs model predictions and feature importances.

## Requirements

- Python 3.7+
- `pandas`
- `numpy`
- `scikit-learn`
- `xgboost`
- `joblib`
- `matplotlib` (for optional plotting)

Install dependencies using pip:

```bash
pip install pandas numpy scikit-learn xgboost joblib matplotlib
```

## Usage

```bash
python final_models_22features_english.py
```

Make sure the data files referenced in the script are available in the specified paths.

## Authors

- [Your Name or Lab]
- Originally written in Chinese, translated and restructured for GitHub sharing.

## License

MIT License (or specify your own)
