"""
DominL Intermediate Analysis Functions
For publication reporting - implements all analysis functions described in METHODS_SUPPLEMENT.md.
Reference: Zhang et al. (2025) bioRxiv doi:10.1101/2025.05.07.652751
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, confusion_matrix, precision_recall_curve, average_precision_score
from scipy.stats import entropy


# =============================================================================
# 1. KL DIVERGENCE FOR FEATURE SELECTION
# =============================================================================

def compute_kl_divergence(train_data, real_data, feature_cols, n_bins=50):
    """
    Compute KL divergence between training and empirical (real) distributions per feature.
    Used for feature selection: retain features where KL < threshold (e.g., 2.0).

    Interpretation:
      KL < 0.5:  High similarity
      0.5-2.0:   Moderate differences
      KL >= 2.0: Large differences, consider excluding
    """
    kl_divergences = {}
    for col in feature_cols:
        if col not in train_data.columns or col not in real_data.columns:
            continue
        t = train_data[col].dropna()
        r = real_data[col].dropna()
        if len(t) < 2 or len(r) < 2:
            continue
        # Common bins across both
        all_vals = np.concatenate([t.values, r.values])
        bins = np.histogram_bin_edges(all_vals, bins=n_bins)
        train_hist, _ = np.histogram(t, bins=bins, density=True)
        real_hist, _ = np.histogram(r, bins=bins, density=True)
        # Smooth to avoid zeros
        train_hist = train_hist + 1e-10
        real_hist = real_hist + 1e-10
        train_hist = train_hist / train_hist.sum()
        real_hist = real_hist / real_hist.sum()
        kl_divergences[col] = entropy(train_hist, real_hist)
    return kl_divergences


# =============================================================================
# 2. ROC AND MODEL COMPARISON
# =============================================================================

def _unpack_model(model_obj):
    """Support dict {model, features} or fitted model with feature_names_in_."""
    if isinstance(model_obj, dict):
        return model_obj["model"], model_obj["features"]
    if not hasattr(model_obj, "feature_names_in_"):
        raise ValueError("Model needs feature_names_in_. Refit or save with feature list.")
    return model_obj, list(model_obj.feature_names_in_)


def compare_models_roc(model_new, model_old, new_name, old_name, test_data, mu, sigma,
                       fig_name="ROC comparison", save_path=None, threshold=0.9):
    """
    Compare two models via ROC curves. Report TPR/FPR at fixed probability threshold.
    """
    model_new, feat_new = _unpack_model(model_new)
    model_old, feat_old = _unpack_model(model_old)
    y_test = test_data["dominance"]

    plt.figure(figsize=(8, 7))

    def eval_and_plot(model, features, label):
        X = (test_data[features] - mu[features]) / sigma[features]
        y_prob = model.predict_proba(X)[:, 1]
        fpr, tpr, _ = roc_curve(y_test, y_prob)
        roc_auc = auc(fpr, tpr)
        y_pred_thr = (y_prob >= threshold).astype(int)
        tn, fp, fn, tp = confusion_matrix(y_test, y_pred_thr).ravel()
        tpr_thr = tp / (tp + fn) if (tp + fn) > 0 else 0
        fpr_thr = fp / (fp + tn) if (fp + tn) > 0 else 0
        plt.plot(fpr, tpr, linewidth=2,
                 label=f"{label} (AUC={roc_auc:.3f})\nthr={threshold}: TPR={tpr_thr:.3f}, FPR={fpr_thr:.3f}")

    eval_and_plot(model_new, feat_new, new_name)
    eval_and_plot(model_old, feat_old, old_name)
    plt.plot([0, 1], [0, 1], "k--", alpha=0.6)
    plt.xlim([0, 1])
    plt.ylim([0, 1.05])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(fig_name)
    plt.legend(loc="lower right", fontsize=10)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
    plt.close()


def compare_multiple_models_roc(model_bundles, test_data, mu, sigma, fig_name="ROC comparison",
                                save_path=None, threshold=0.9):
    """
    model_bundles: list of (model_bundle, model_name)
    """
    plt.figure(figsize=(8, 7))
    for model_bundle, model_name in model_bundles:
        model, features = _unpack_model(model_bundle)
        X = (test_data[features] - mu[features]) / sigma[features]
        y_prob = model.predict_proba(X)[:, 1]
        fpr, tpr, _ = roc_curve(test_data["dominance"], y_prob)
        roc_auc = auc(fpr, tpr)
        y_pred_thr = (y_prob >= threshold).astype(int)
        tn, fp, fn, tp = confusion_matrix(test_data["dominance"], y_pred_thr).ravel()
        tpr_thr = tp / (tp + fn) if (tp + fn) > 0 else 0
        fpr_thr = fp / (fp + tn) if (fp + tn) > 0 else 0
        plt.plot(fpr, tpr, linewidth=2,
                 label=f"{model_name} (AUC={roc_auc:.3f})\nthr={threshold}: TPR={tpr_thr:.3f}, FPR={fpr_thr:.3f}")
    plt.plot([0, 1], [0, 1], "k--", alpha=0.6)
    plt.xlim([0, 1])
    plt.ylim([0, 1.05])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(fig_name)
    plt.legend(loc="lower right", fontsize=9)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
    plt.close()


def compute_tpr_fpr_at_threshold(y_true, y_prob, threshold):
    """Return TPR and FPR at a fixed probability threshold."""
    y_pred = (y_prob >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    tpr = tp / (tp + fn) if (tp + fn) > 0 else 0
    fpr = fp / (fp + tn) if (fp + tn) > 0 else 0
    return tpr, fpr


# =============================================================================
# 3. RECESSIVE PERCENTAGE BY DOMINANCE
# =============================================================================

def calculate_recessive_percentage(data, model, features, mu, sigma, thresholds=[0.5, 0.9]):
    """
    For each dominance value, compute % of windows predicted as recessive.
    Returns dict: {threshold: {dominance_values, recessive_percentages}}
    """
    results = {}
    X = (data[features] - mu[features]) / sigma[features]
    dominance_values = sorted(data["dominance"].unique())
    for thresh in thresholds:
        pred = (model.predict_proba(X)[:, 1] >= thresh).astype(int)
        pcts = []
        for v in dominance_values:
            subset = pred[data["dominance"] == v]
            pcts.append(subset.mean() * 100)
        results[thresh] = {"dominance_values": dominance_values, "recessive_percentages": pcts}
    return results


def plot_recessive_percentages(results, save_path="recessive_percentage_plot.png",
                               figure_title="Percentage of Recessive Predictions at Different Thresholds"):
    plt.figure(figsize=(10, 8))
    for threshold, data in results.items():
        plt.scatter(data["dominance_values"], data["recessive_percentages"], label=f"Prob >= {threshold}")
        plt.plot(data["dominance_values"], data["recessive_percentages"])
    plt.xlabel("Dominance Values")
    plt.ylabel("Percentage Predicted as Recessive")
    plt.title(figure_title)
    plt.legend()
    plt.savefig(save_path)
    plt.close()


# =============================================================================
# 4. POWER, FPR, SPECIFICITY
# =============================================================================

def calculate_metrics(y_true, y_pred):
    """Compute power (TPR), FPR, specificity (TNR)."""
    cm = confusion_matrix(y_true, y_pred)
    if cm.size != 4:
        return np.nan, np.nan, np.nan  # Need both classes
    tn, fp, fn, tp = cm.ravel()
    power = tp / (tp + fn) if (tp + fn) > 0 else 0
    fpr = fp / (fp + tn) if (fp + tn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    return power, fpr, specificity


# =============================================================================
# 5. PERCENTILE-BASED ANCESTRY MODELS
# =============================================================================

def make_percentile_predictor(train_data, test_data, feature_name):
    """
    Create a predictor using percentile rank of a single feature.
    Used for Neanderthal ancestry-only comparison models.
    """
    train_data = train_data.copy()
    test_data = test_data.copy()
    train_data[f"{feature_name}_pct"] = train_data[feature_name].rank(pct=True)
    test_data[f"{feature_name}_pct"] = test_data[feature_name].rank(pct=True)
    return train_data, test_data, f"{feature_name}_pct"


# =============================================================================
# 6. CROSS-DATASET ROC/PR
# =============================================================================

def compare_model_across_datasets_roc_pr(model_bundle, datasets_dict, original_name, mu, sigma,
                                         fig_name="XGBoost across datasets",
                                         roc_save_path=None, pr_save_path=None, threshold=0.9):
    """
    Apply one model to multiple datasets, plot ROC and Precision-Recall.
    datasets_dict: {name: DataFrame with dominance and feature cols}
    """
    model, features = _unpack_model(model_bundle)
    plt.figure(figsize=(8, 7))
    for name, df in datasets_dict.items():
        X = (df[features] - mu[features]) / sigma[features]
        y = df["dominance"]
        y_prob = model.predict_proba(X)[:, 1]
        fpr, tpr, _ = roc_curve(y, y_prob)
        roc_auc = auc(fpr, tpr)
        tpr_thr, fpr_thr = compute_tpr_fpr_at_threshold(y, y_prob, threshold)
        plt.plot(fpr, tpr, linewidth=2,
                 label=f"{name} (AUC={roc_auc:.3f})\nthr={threshold}: TPR={tpr_thr:.3f}, FPR={fpr_thr:.3f}")
    plt.plot([0, 1], [0, 1], "k--", alpha=0.6)
    plt.xlim([0, 1])
    plt.ylim([0, 1.05])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"{fig_name} — ROC")
    plt.legend(loc="lower right", fontsize=9)
    plt.tight_layout()
    if roc_save_path:
        plt.savefig(roc_save_path, dpi=300)
    plt.close()

    # Precision-Recall
    plt.figure(figsize=(8, 7))
    for name, df in datasets_dict.items():
        X = (df[features] - mu[features]) / sigma[features]
        y = df["dominance"]
        y_prob = model.predict_proba(X)[:, 1]
        precision, recall, _ = precision_recall_curve(y, y_prob)
        ap = average_precision_score(y, y_prob)
        plt.plot(recall, precision, linewidth=2, label=f"{name} (AP={ap:.3f})")
    plt.xlim([0, 1])
    plt.ylim([0, 1.05])
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title(f"{fig_name} — Precision-Recall")
    plt.legend(loc="lower left", fontsize=9)
    plt.tight_layout()
    if pr_save_path:
        plt.savefig(pr_save_path, dpi=300)
    plt.close()
