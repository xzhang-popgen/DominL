# DominL Methods Supplement: Intermediate Analysis Functions

This document catalogues all analysis functions performed for the DominL method, including those not part of the core reproducible pipeline. These are intended for **Methods** or **Supplementary Materials** documentation.

**Reference:** Zhang et al. (2025) *Neanderthal introgressed ancestry reveals human genomic regions enriched with recessive deleterious mutations*  
bioRxiv doi: 10.1101/2025.05.07.652751

---

## 1. Data Preprocessing

### 1.1 Missing value handling

- **Columns with missing values:** `introg_anc_window`, `D`, `Q95`, `patterson_f3`
- **Procedure:** For `introg_anc_window`, `D`, `Q95`, missing values were filled within 5-row windows by the mean of that window
- **patterson_f3:** Excluded from feature set due to high missingness (750,972 / 750,972 rows)
- **Implemented in:** `training/train_xgboost.py` (fillna handling)

### 1.2 Train–test split

- **Stratification:** By `exon_density` (e.g., ≥800 normalized value)
- **Split:** 80% train, 20% test (random state 42)
- **Additional processing:** Correlation-based feature removal (drop features with correlation > 0.8)
- **Implemented in:** `training/train_xgboost.py`, `training/train_compare_models.py`

### 1.3 Normalization

- **Training/test:** Z-score normalization using training data: `(x - μ_train) / σ_train` (population std, ddof=0)
- **Replace σ=0 with 1** to avoid division by zero
- **Implemented in:** `training/train_xgboost.py`, `training/train_compare_models.py`

---

## 2. Feature Selection

### 2.1 KL divergence for distribution overlap

**Purpose:** Retain features where empirical and simulation distributions overlap; exclude those with large distributional mismatch.

**Method:**
1. Discretize each continuous feature into histograms (bins) for simulation (training) and empirical data
2. Compute Kullback–Leibler divergence: `KL(P||Q) = Σ P(x) log(P(x)/Q(x))` using `scipy.stats.entropy`
3. Interpret thresholds:
   - **KL < 0.5:** High similarity; features likely contribute consistently
   - **0.5 ≤ KL < 2.0:** Moderate differences; may introduce some variability
   - **KL ≥ 2.0:** Large differences; consider excluding or further investigation

**Features with large KL (typically excluded):** `df_p3_p2`, `hap_diversity_p3`, `garud_h1`, `garud_h12`, `garud_h2_h1`, `fD`, `RD`, `mean_recrate`, `num_variant_window` (when comparing raw simulation vs empirical)

**Implemented in:** `analysis_functions/analysis_functions.py` (`compute_kl_divergence`)

### 2.2 Feature importance (XGBoost)

- Rank features by XGBoost `feature_importances_`
- Used to prioritize features and validate KL-based selection
- **Implemented in:** `training/train_xgboost.py`, `training/train_compare_models.py`

### 2.3 Other feature selection methods (stat2ML.py)

- **SelectKBest (chi2):** Univariate feature selection
- **RFE (Recursive Feature Elimination):** Wrapper method
- **Ridge/Lasso regularization:** Coefficient-based selection
- **Implemented in:** `slim/stat2ML.py` (legacy; includes SelectKBest, RFE, Ridge/Lasso)

---

## 3. Model Comparison

### 3.1 22- vs 20-feature XGBoost

- **20 features:** Exclude `mean_recrate`, `recrate_window`
- **22 features:** Include recombination rate features
- **Comparison:** ROC curves, TPR/FPR at probability threshold 0.9
- **Stratification:** By exon density (≥200, ≥400, ≥600, ≥800 per 5MB)
- **Implemented in:** `training/train_compare_models.py`, `analysis_functions/analysis_functions.py`

### 3.2 Neanderthal ancestry–only models

**Percentile-based predictors:**
- `mean_introg_anc_pct` = percentile rank of `mean_introg_anc`
- `introg_anc_window_pct` = percentile rank of `introg_anc_window`
- Used as single-feature “models” for comparison with full XGBoost
- **Implemented in:** `analysis_functions/analysis_functions.py` (`make_percentile_predictor`, `compare_multiple_models_roc`)

### 3.3 Four-model ROC comparison

- XGBoost (22 features)
- Percentile model: `mean_introg_anc`
- Percentile model: `introg_anc_window`
- Compare via `compare_multiple_models_roc`
- **Implemented in:** `analysis_functions/analysis_functions.py` (`make_percentile_predictor`, `compare_multiple_models_roc`)

### 3.4 Multi-dataset ROC and precision–recall

**Datasets:** Original test, high recombination (`highr`), low recombination (`lowr`), positive selection (`possel`), neutral, Varyh (varying dominance)

**Functions:**
- `compare_model_across_datasets_roc`: ROC curves across datasets
- `compare_model_across_datasets_roc_pr`: ROC and Precision–Recall curves
- **Implemented in:** `analysis_functions/analysis_functions.py` (`compare_model_across_datasets_roc_pr`)

### 3.5 FPR on neutral data

- `compute_fpr_neutral_only`: Compute FPR when applying model to neutral (no selection) simulations
- Used to mark FPR of neutral data on ROC plots
- **Implemented in:** `analysis_functions/analysis_functions.py` (`compute_tpr_fpr_at_threshold`, `calculate_metrics`)

### 3.6 TPR/FPR at fixed threshold

- `compute_tpr_fpr_at_threshold(y_true, y_prob, threshold)`: Returns TPR and FPR at a given probability cutoff
- **Implemented in:** `analysis_functions/analysis_functions.py` (`compute_tpr_fpr_at_threshold`, `calculate_metrics`)

---

## 4. Stratification Analyses

### 4.1 Exon density stratification

- Subset test/empirical data by `exon_density` thresholds: 200, 400, 600, 800 (exons per 5MB)
- High exon density (≥600) corresponds to best model accuracy
- **Implemented in:** `analysis_functions/analysis_functions.py`, `empirical_prediction/predict_empirical.py`

### 4.2 Power, FPR, and specificity

- **Power (TPR):** tp / (tp + fn)
- **FPR:** fp / (fp + tn)
- **Specificity (TNR):** tn / (tn + fp)
- Computed at probability threshold 0.9
- **Implemented in:** `analysis_functions/analysis_functions.py` (`calculate_metrics`)

---

## 5. ML Model Comparison (13 models)

**Models evaluated in `training/train_compare_models.py`:**

| Abbreviation | Model |
|--------------|-------|
| ETC | Extra Trees Classifier |
| RF | Random Forest |
| L0LR | L0-penalized (multinomial) Logistic Regression |
| L1LR | L1-penalized Logistic Regression |
| L2LR | L2-penalized Logistic Regression |
| MLP | Multi-layer Perceptron |
| RBF | RBF-kernel SVM |
| XGBoost | XGBoost (selected for final model) |
| CatBoost | (from empirical prediction outputs) |
| LightGBM | (from empirical prediction outputs) |
| Elastic Net, Lasso, Gradient Boosting | (from final_models) |

**Procedure:** Train each on same features, compare ROC, precision–recall, confusion matrices.

---

## 6. Recessive percentage by dominance value

- **Function:** `calculate_recessive_percentage(empirical_data, model, features, mu, sigma, thresholds=[0.5, 0.9])`
- For each dominance value (e.g., in Varyh test data), compute % of windows predicted as recessive
- **Plot:** `plot_recessive_percentages(results, save_path, figure_title)`
- **Implemented in:** `analysis_functions/analysis_functions.py` (`calculate_recessive_percentage`, `plot_recessive_percentages`)

---

## 7. Memory and timing monitoring

- **Function:** `monitor_memory(pid, interval, stats)` – records peak RSS during training
- Used for computational resource reporting
- **Implemented in:** `analysis_functions/analysis_functions.py` (`calculate_recessive_percentage`, `plot_recessive_percentages`)

---

## 8. Confusion matrix and precision–recall

- **Function:** `prec_recall_mod(cm, num_classes)` – precision and recall per class
- **Plots:** `plot_prcurve_average_mod`, `plot_roc_mod` – multi-class ROC and precision–recall
- **Implemented in:** `slim/stat2ML.py` (legacy; includes SelectKBest, RFE, Ridge/Lasso)

---

## 9. Feature importance plotting

- **Function:** `plot_featurescore(scores_file, savefile)` – horizontal bar plot of feature importances
- **Function:** `write_importance(outfile, feature_list, score_list)` – write importance to file
- **Implemented in:** `slim/stat2ML.py` (legacy; includes SelectKBest, RFE, Ridge/Lasso)

---

## 10. Pipeline files and analyses

| File | Main analyses |
|------|----------------|
| `training/train_xgboost.py` | XGBoost training, normalization, 22-feature model |
| `training/train_compare_models.py` | 13-model comparison (ETC, RF, LR, MLP, SVM, XGBoost, etc.) |
| `empirical_prediction/predict_empirical.py` | Apply model to empirical data, column mapping |
| `analysis_functions/analysis_functions.py` | KL divergence, ROC/PR comparison, stratification metrics |
| `slim/sim2vcf2stats.py` | SLiM simulation, VCF extraction, stat computation |
| `slim/stat2ML.py` | Legacy: SelectKBest, RFE, Ridge/Lasso, ETC, RF, MLP, SVM |

---

## 12. Implemented analysis functions (runnable)

The module `analysis_functions/analysis_functions.py` provides executable implementations of the main analysis functions:

| Function | Description |
|----------|-------------|
| `compute_kl_divergence` | KL divergence per feature between train and empirical |
| `compare_models_roc` | ROC comparison of two models with TPR/FPR at threshold |
| `compare_multiple_models_roc` | ROC comparison of multiple models |
| `compute_tpr_fpr_at_threshold` | TPR and FPR at fixed probability cutoff |
| `calculate_recessive_percentage` | % recessive by dominance value |
| `plot_recessive_percentages` | Plot recessive % vs dominance |
| `calculate_metrics` | Power, FPR, specificity from confusion matrix |
| `make_percentile_predictor` | Percentile-rank single-feature predictor |
| `compare_model_across_datasets_roc_pr` | ROC and PR across multiple datasets |

---

## 11. Output files for publication

| Analysis | Typical output |
|----------|----------------|
| 22 vs 20 feature ROC | `1.1_roc_22_vs_20.png`, `1.1_roc_22_vs_20_exondensity_{200,400,600,800}.png` |
| Ancestry-only ROC | `1.2_roc_introg_anc_window_vs_mean_introg_anc*.png`, `1.2_roc_4_models*.png` |
| Cross-dataset ROC/PR | `2.1_roc_across_datasets*.png`, `2.2_roc_across_datasets_w_FPR_neutral.png` |
| Empirical predictions | `2.5_model_22_vs_20_pred_{POP}.csv`, `*_exondensity_{200,400,600,800}.csv` |
| KL divergence | Printed/saved per feature |
| Model metrics | `model_comparison_metrics.xlsx` (power, FPR, specificity) |
