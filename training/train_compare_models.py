#!/usr/bin/env python3
"""
DominL Model Comparison - Train and compare 13 ML models for recessive vs additive classification.
Used for model selection; XGBoost was selected for best performance.
Reference: Zhang et al. (2025) bioRxiv doi:10.1101/2025.05.07.652751
"""

import argparse
import os
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, accuracy_score, confusion_matrix
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC
from sklearn.ensemble import GradientBoostingClassifier

# Add parent for config
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from config.features import SELECTED_FEATURES_22, COLUMNS_TO_DROP

# Optional: CatBoost, LightGBM, ElasticNet
try:
    import lightgbm as lgb
    HAS_LIGHTGBM = True
except ImportError:
    HAS_LIGHTGBM = False
try:
    import catboost as cb
    HAS_CATBOOST = True
except ImportError:
    HAS_CATBOOST = False

import xgboost as xgb


def get_all_models():
    """Return dict of (name, model_instance) for the 13 models."""
    models = {
        "ETC": ExtraTreesClassifier(n_estimators=100, random_state=0, max_features="sqrt",
                                    min_samples_split=10, min_samples_leaf=10),
        "RF": RandomForestClassifier(n_estimators=100, random_state=0, max_depth=2),
        "L0LR": LogisticRegression(solver='saga', multi_class='multinomial', penalty='l1',
                                   tol=0.01, C=1e10, random_state=0),
        "L1LR": LogisticRegression(solver='saga', multi_class='multinomial', penalty='l1', random_state=0),
        "L2LR": LogisticRegression(solver='lbfgs', multi_class='multinomial', penalty='l2', random_state=0),
        "MLP": MLPClassifier(solver='lbfgs', alpha=1e-5, hidden_layer_sizes=(512, 256, 128), random_state=1),
        "RBF_SVM": SVC(kernel='rbf', class_weight='balanced', probability=True, random_state=0),
        "XGBoost": xgb.XGBClassifier(eval_metric='logloss', random_state=0),
        "GradientBoosting": GradientBoostingClassifier(n_estimators=100, random_state=0),
        "Logistic": LogisticRegression(solver='lbfgs', max_iter=1000, random_state=0),
        "Lasso": LogisticRegression(solver='saga', penalty='l1', random_state=0),
        "ElasticNet": LogisticRegression(solver='saga', penalty='elasticnet', l1_ratio=0.5, random_state=0),
    }
    if HAS_LIGHTGBM:
        models["LightGBM"] = lgb.LGBMClassifier(n_estimators=100, random_state=0, verbose=-1)
    if HAS_CATBOOST:
        models["CatBoost"] = cb.CatBoostClassifier(iterations=100, random_state=0, verbose=0)
    return models


def load_data(train_path, test_path, feature_list):
    """Load and prepare train/test data."""
    train = pd.read_csv(train_path) if train_path.endswith('.csv') else pd.read_excel(train_path)
    test = pd.read_csv(test_path) if test_path.endswith('.csv') else pd.read_excel(test_path)
    for col in COLUMNS_TO_DROP:
        if col in train.columns:
            train = train.drop(columns=[col])
        if col in test.columns:
            test = test.drop(columns=[col])
    train = train[feature_list + ['dominance']].dropna()
    test = test[feature_list + ['dominance']].dropna()
    return train, test


def train_and_evaluate(train_data, test_data, feature_list, output_dir=None):
    """
    Train all models, evaluate on test set, return comparison DataFrame.
    """
    X_train = train_data[feature_list]
    y_train = train_data['dominance'].astype(int)
    X_test = test_data[feature_list]
    y_test = test_data['dominance'].astype(int)

    # Normalization
    train_mean = X_train.mean()
    train_std = X_train.std(ddof=0).replace(0, 1)
    X_train_norm = (X_train - train_mean) / train_std
    X_test_norm = (X_test - train_mean) / train_std

    results = []
    trained_models = {}

    for name, model in get_all_models().items():
        try:
            model.fit(X_train_norm, y_train)
            y_prob = model.predict_proba(X_test_norm)[:, 1]
            y_pred = (y_prob >= 0.5).astype(int)
            roc_auc = roc_auc_score(y_test, y_prob)
            acc = accuracy_score(y_test, y_pred)
            cm = confusion_matrix(y_test, y_pred)
            if cm.size != 4:
                fpr, tpr = np.nan, np.nan
            else:
                tn, fp, fn, tp = cm.ravel()
                fpr = fp / (fp + tn) if (fp + tn) > 0 else 0
                tpr = tp / (tp + fn) if (tp + fn) > 0 else 0
            results.append({
                'model': name,
                'roc_auc': roc_auc,
                'accuracy': acc,
                'TPR': tpr,
                'FPR': fpr
            })
            trained_models[name] = {
                'model': model,
                'train_mean': train_mean,
                'train_std': train_std,
                'features': feature_list
            }
        except Exception as e:
            print(f"  {name}: failed ({e})")
            results.append({'model': name, 'roc_auc': np.nan, 'accuracy': np.nan, 'TPR': np.nan, 'FPR': np.nan})

    df = pd.DataFrame(results)
    df = df.sort_values('roc_auc', ascending=False).reset_index(drop=True)

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        df.to_csv(os.path.join(output_dir, 'model_comparison_results.csv'), index=False)
        # Save best model (XGBoost) for downstream use
        if 'XGBoost' in trained_models:
            import joblib
            b = trained_models['XGBoost']
            joblib.dump(b['model'], os.path.join(output_dir, 'xgboost_model_22features.joblib'))
            norm_params = pd.DataFrame({
                'feature': b['features'],
                'mean': b['train_mean'][b['features']].values,
                'std': b['train_std'][b['features']].values
            })
            norm_params.to_csv(os.path.join(output_dir, 'normalization_params.csv'), index=False)

    return df, trained_models


def main():
    parser = argparse.ArgumentParser(description='Train and compare 13 ML models')
    parser.add_argument('--train', required=True, help='Path to training data')
    parser.add_argument('--test', required=True, help='Path to test data')
    parser.add_argument('--output', default='./models', help='Output directory')
    parser.add_argument('--features', choices=['20', '22'], default='22')
    args = parser.parse_args()

    from config.features import SELECTED_FEATURES_20
    feature_list = SELECTED_FEATURES_22 if args.features == '22' else SELECTED_FEATURES_20

    train_data, test_data = load_data(args.train, args.test, feature_list)
    print(f"Training: {len(train_data)}, Test: {len(test_data)}")

    df, _ = train_and_evaluate(train_data, test_data, feature_list, args.output)
    print("\n=== Model comparison (sorted by ROC-AUC) ===")
    print(df.to_string(index=False))
    print(f"\nBest model: {df.iloc[0]['model']} (ROC-AUC={df.iloc[0]['roc_auc']:.3f})")
    print(f"Results saved to {args.output}/model_comparison_results.csv")


if __name__ == '__main__':
    main()
