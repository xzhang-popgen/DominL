#!/usr/bin/env python3
"""
DominL Training Pipeline - XGBoost model training.
Trains XGBoost classifier with 22 features to distinguish recessive vs additive 1MB windows.
Uses normalization: (data - train_mean) / train_std
Reference: Zhang et al. (2025) bioRxiv doi:10.1101/2025.05.07.652751
"""

import argparse
import os
import pandas as pd
import xgboost as xgb
import joblib
import numpy as np

# Add parent for config
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from config.features import SELECTED_FEATURES_22, SELECTED_FEATURES_20, COLUMNS_TO_DROP


def load_train_test(train_path, test_path):
    """Load training and test data, drop irrelevant columns."""
    train_data = pd.read_csv(train_path) if train_path.endswith('.csv') else pd.read_excel(train_path)
    test_data = pd.read_csv(test_path) if test_path.endswith('.csv') else pd.read_excel(test_path)

    for col in COLUMNS_TO_DROP:
        if col in train_data.columns:
            train_data = train_data.drop(columns=[col])
        if col in test_data.columns:
            test_data = test_data.drop(columns=[col])

    return train_data, test_data


def train_model(train_data, feature_list=SELECTED_FEATURES_22, output_dir=None):
    """
    Train XGBoost model with selected features.
    Returns model, train_mean, train_std for use in prediction.
    """
    # Ensure all features exist
    missing = [f for f in feature_list if f not in train_data.columns]
    if missing:
        raise ValueError(f"Missing features in training data: {missing}")

    train_subset = train_data[feature_list + ['dominance']].dropna()
    feature_cols = [f for f in feature_list if f in train_subset.columns]
    X_train = train_subset[feature_cols]
    y_train = train_subset['dominance'].astype(int)

    # Normalization (population std ddof=0)
    train_mean = X_train.mean()
    train_std = X_train.std(ddof=0).replace(0, 1)

    X_train_norm = (X_train - train_mean[feature_cols]) / train_std[feature_cols]

    model = xgb.XGBClassifier(eval_metric='logloss')
    model.fit(X_train_norm, y_train)

    # Save model and normalization params
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        joblib.dump(model, os.path.join(output_dir, 'xgboost_model_22features.joblib'))
        norm_params = pd.DataFrame({
            'feature': feature_cols,
            'mean': train_mean[feature_cols].values,
            'std': train_std[feature_cols].values
        })
        norm_params.to_csv(os.path.join(output_dir, 'normalization_params.csv'), index=False)

    return model, train_mean, train_std, feature_cols


def main():
    parser = argparse.ArgumentParser(description='Train DominL XGBoost model')
    parser.add_argument('--train', required=True, help='Path to training data (csv or xlsx)')
    parser.add_argument('--test', help='Path to test data (optional, for evaluation)')
    parser.add_argument('--output', default='./models', help='Output directory for model and params')
    parser.add_argument('--features', choices=['20', '22'], default='22',
                        help='Number of features (20 or 22)')
    args = parser.parse_args()

    feature_list = SELECTED_FEATURES_22 if args.features == '22' else SELECTED_FEATURES_20

    train_data = pd.read_csv(args.train) if args.train.endswith('.csv') else pd.read_excel(args.train)
    for col in COLUMNS_TO_DROP:
        if col in train_data.columns:
            train_data = train_data.drop(columns=[col])

    test_data = None
    if args.test and os.path.exists(args.test):
        test_data = pd.read_csv(args.test) if args.test.endswith('.csv') else pd.read_excel(args.test)
        for col in COLUMNS_TO_DROP:
            if col in test_data.columns:
                test_data = test_data.drop(columns=[col])

    print(f"Training samples: {len(train_data)}")
    model, train_mean, train_std, feature_cols = train_model(
        train_data, feature_list=feature_list, output_dir=args.output
    )
    print(f"Model saved to {args.output}")

    if test_data is not None:
        from sklearn.metrics import roc_auc_score, accuracy_score
        test_subset = test_data[feature_cols + ['dominance']].dropna()
        X_test = (test_subset[feature_cols] - train_mean[feature_cols]) / train_std[feature_cols]
        y_test = test_subset['dominance'].astype(int)
        y_prob = model.predict_proba(X_test)[:, 1]
        y_pred = (y_prob >= 0.5).astype(int)
        print(f"Test ROC-AUC: {roc_auc_score(y_test, y_prob):.3f}")
        print(f"Test Accuracy: {accuracy_score(y_test, y_pred):.3f}")


if __name__ == '__main__':
    main()
