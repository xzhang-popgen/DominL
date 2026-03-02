#!/usr/bin/env python3
"""
DominL Empirical Prediction Pipeline.
Apply trained XGBoost model to 1KG non-African population empirical summary statistics.
Reference: Zhang et al. (2025) bioRxiv doi:10.1101/2025.05.07.652751
"""

import argparse
import os
import pandas as pd
import numpy as np
import joblib

# Add parent for config
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from config.features import SELECTED_FEATURES_22, EMPIRICAL_COLUMN_MAPPING


def load_model_and_params(model_dir):
    """Load trained model and normalization parameters."""
    model_path = os.path.join(model_dir, 'xgboost_model_22features.joblib')
    params_path = os.path.join(model_dir, 'normalization_params.csv')
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Model not found: {model_path}")
    model = joblib.load(model_path)
    norm_params = pd.read_csv(params_path)
    train_mean = norm_params.set_index('feature')['mean']
    train_std = norm_params.set_index('feature')['std'].replace(0, 1)
    feature_cols = norm_params['feature'].tolist()
    return model, train_mean, train_std, feature_cols


def prepare_empirical_data(empirical_path, sep='\t'):
    """Load and prepare empirical data with column renaming."""
    df = pd.read_csv(empirical_path, sep=sep)
    df = df.rename(columns=EMPIRICAL_COLUMN_MAPPING)
    return df


def predict_empirical(empirical_path, model_dir, output_path=None, sep='\t',
                     exon_thresholds=[200, 400, 600, 800]):
    """
    Run prediction on empirical data.
    Returns DataFrame with prediction probabilities and optional exon-density subsets.
    """
    model, train_mean, train_std, feature_cols = load_model_and_params(model_dir)
    empirical = prepare_empirical_data(empirical_path, sep=sep)

    # Check for required features
    missing = [f for f in feature_cols if f not in empirical.columns]
    if missing:
        raise ValueError(f"Missing features in empirical data: {missing}. "
                         f"Available: {list(empirical.columns)}")

    X = empirical[feature_cols].copy()
    X = X.fillna(X.median())  # Fill NaN with median
    X_norm = (X - train_mean[feature_cols]) / train_std[feature_cols]

    prob_recessive = model.predict_proba(X_norm)[:, 1]
    pred_recessive = (prob_recessive >= 0.5).astype(int)

    result = empirical.copy()
    result['prediction_probability'] = prob_recessive
    result['predicted_recessive'] = pred_recessive

    if output_path:
        out_dir = os.path.dirname(output_path)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        result.to_csv(output_path, index=False)
        print(f"Full predictions saved to {output_path}")

        # Subsets by exon density
        if 'exon_density' in result.columns:
            for t in exon_thresholds:
                subset = result[result['exon_density'] >= t]
                sub_path = output_path.replace('.csv', f'_exondensity_{t}.csv')
                subset[['prediction_probability', 'predicted_recessive']].to_csv(sub_path, index=False)
                print(f"  exon_density >= {t}: {len(subset)} windows -> {sub_path}")

    return result


def main():
    parser = argparse.ArgumentParser(description='DominL empirical prediction')
    parser.add_argument('--empirical', required=True, help='Path to empirical stats file (tab or csv)')
    parser.add_argument('--model-dir', required=True, help='Directory containing trained model')
    parser.add_argument('--output', help='Output path for predictions')
    parser.add_argument('--sep', default='\t', help='Separator for empirical file')
    args = parser.parse_args()

    output = args.output or args.empirical.replace('.txt', '_predictions.csv').replace('.csv', '_predictions.csv')
    result = predict_empirical(args.empirical, args.model_dir, output_path=output, sep=args.sep)

    # Summary
    pct_recessive = result['predicted_recessive'].mean() * 100
    print(f"\nSummary: {len(result)} windows, {pct_recessive:.1f}% predicted as recessive (prob >= 0.5)")
    print(result['prediction_probability'].describe())


if __name__ == '__main__':
    main()
