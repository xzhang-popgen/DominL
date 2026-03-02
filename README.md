# DominL Integrated Pipeline

Full pipeline for the **DominL** method: a machine learning approach to identify genomic regions enriched with recessive deleterious mutations using Neanderthal introgression patterns.

**Reference:** Zhang et al. (2025) *Neanderthal introgressed ancestry reveals human genomic regions enriched with recessive deleterious mutations*  
bioRxiv doi: [10.1101/2025.05.07.652751](https://www.biorxiv.org/content/10.1101/2025.05.07.652751v1)

## Overview

1. **Simulations** (SLiM): Simulate 5MB genomic segments under two dominance regimes (h=0.5 additive, h=0 recessive) with Eurasian demography + Neanderthal introgression
2. **Training**: Extract 1MB-window summary statistics, train XGBoost to classify recessive vs additive
3. **Empirical prediction**: Apply trained model to 1000 Genomes non-African population statistics

### Key specifications

- **22 features** (includes `mean_recrate`, `recrate_window`) – recommended model
- **Normalization**: `(x - train_mean) / train_std` with train statistics
- **XGBoost** as primary classifier (best of 13 models tested; see `training/train_compare_models.py`)

---

## Directory structure

```
├── config/
│   ├── features.py         # Feature lists, column mappings
│   └── paths.py            # Path configuration
├── slim/
│   ├── Gravel_Eurasian_varyh.slim
│   ├── sim2vcf2stats.py    # Run SLiM, extract stats
│   ├── stat2ML.py          # Legacy ML script (use training/ instead)
│   ├── 1kseg_byexon/       # Segment definition files
│   └── README.md
├── training/
│   ├── train_xgboost.py    # XGBoost training (primary)
│   └── train_compare_models.py  # 13-model comparison
├── empirical_prediction/
│   └── predict_empirical.py
├── analysis_functions/
│   └── analysis_functions.py   # KL divergence, ROC comparison, etc.
├── run_pipeline.sh
├── requirements.txt
├── METHODS_SUPPLEMENT.md
└── README.md
```

---

## Data availability

**Google Drive supplementary data:** Training data used in DominL training, and empirical predictions on 7 1000 Genomes populations can be found:

https://drive.google.com/drive/u/2/folders/1sOg0dvZ4KUqurTOFA8gsKldPEJiKBEf8

---

## Requirements

```
pandas numpy scikit-learn xgboost joblib openpyxl
```

For simulations: `msprime`, `pyslim`, `scikit-allel`, `bcftools`, SLiM 3.x

---

## Usage

### 1. Training

```bash
python training/train_xgboost.py \
  --train /path/to/train_data_all.xlsx \
  --test /path/to/test_data_all.xlsx \
  --output models \
  --features 22
```

Training data must contain columns: `dominance` (0=additive, 1=recessive) plus the 22 feature columns.

### 2. Empirical prediction

```bash
python empirical_prediction/predict_empirical.py \
  --empirical /path/to/hg19_CHB_empirical-stats_all.txt \
  --model-dir models \
  --output predictions/CHB_predictions.csv
```

Empirical file columns should match (or use `EMPIRICAL_COLUMN_MAPPING` in `config/features.py`):
- `nea_anc` → `introg_anc_window`
- `exon_5mb` → `exon_density`
- `exon` → `exon_window`
- `nea_anc_5mb` → `mean_introg_anc`
- `recomb_5mb` → `mean_recrate`
- `recomb` → `recrate_window`

### 3. Run pipeline (training + prediction)

```bash
export TRAIN_DATA=/path/to/train_data_all.xlsx
export TEST_DATA=/path/to/test_data_all.xlsx
export EMPIRICAL_DATA=/path/to/hg19_CHB_empirical-stats_all.txt
./run_pipeline.sh
```

---

## Reproducibility: file structure, data format, and workflow

### File structure and data dependencies

| Component | Purpose | Data dependencies |
|-----------|---------|-------------------|
| `slim/` | Generate simulation stats | Segment files in `1kseg_byexon/`; outputs stats CSVs |
| `training/` | Train XGBoost model | `train_data.xlsx`, `test_data.xlsx` (or from simulations) |
| `empirical_prediction/` | Apply model to empirical data | Trained model + empirical stats (1KG, HGDP, etc.) |
| `config/features.py` | Feature names, column mapping | Edit `EMPIRICAL_COLUMN_MAPPING` for new datasets |

Pre-computed training data and empirical predictions are available on [Google Drive](#data-availability).

### Workflow overview

```
Simulations (slim/)  →  Training (training/)  →  Empirical prediction (empirical_prediction/)
         │                      │                              │
         │                      │                              │
   stat CSVs per window    train_data.xlsx              empirical_stats.txt
   (or use Google Drive)   test_data.xlsx               (1KG, HGDP, etc.)
```

1. **Simulations:** Run SLiM to generate synthetic data; `sim2vcf2stats.py` extracts per-window stats.
2. **Training:** Train XGBoost on simulation stats with labels (0=additive, 1=recessive).
3. **Empirical prediction:** Apply the trained model to empirical summary statistics from any population.

### Input format for training

Training and test files (CSV or Excel) must include:

| Column | Description |
|--------|-------------|
| `dominance` | 0 = additive, 1 = recessive (2 = partial, if using 3-class) |
| Plus all 22 feature columns (see below) | Per 1MB window |

**Required 22 features:** `exon_density`, `mean_introg_anc`, `exon_window`, `introg_anc_window`, `divergence_p3_p1`, `divergence_p3_p2`, `watterson_theta_p3`, `windowed_tajima_d_p3`, `D`, `Het`, `Q95`, `U0`, `U20`, `U50`, `U80`, `num_seg_p1`, `num_seg_p3`, `num_private_seg_p1`, `num_private_seg_p2`, `num_private_seg_p3`, `mean_recrate`, `recrate_window`

Optional metadata columns (`segment`, `rep`, `start`, `end`, `patterson_f3`) are dropped automatically.

### Input format for empirical prediction

Empirical data (e.g., from 1000 Genomes or HGDP) must have summary statistics per 1MB window. Column names can differ; the pipeline maps them via `config/features.py`:

| If your file has | It maps to |
|------------------|------------|
| `nea_anc` | `introg_anc_window` |
| `exon_5mb` | `exon_density` |
| `exon` | `exon_window` |
| `nea_anc_5mb` | `mean_introg_anc` |
| `recomb_5mb` | `mean_recrate` |
| `recomb` | `recrate_window` |

Other features must match the 22 names exactly. Add mappings in `EMPIRICAL_COLUMN_MAPPING` in `config/features.py` for new column names.

### Training on a new dataset (e.g., HGDP)

To train DominL on a different dataset:

1. **Prepare summary statistics** per 1MB window for your populations:
   - Neanderthal introgression ancestry (`introg_anc_window`, `mean_introg_anc`)
   - Exon counts (`exon_window`, `exon_density`)
   - Recombination (`mean_recrate`, `recrate_window`)
   - D, fD, diversity, divergence, etc. (see `SELECTED_FEATURES_22` in `config/features.py`)

2. **Add labels:** For simulation-based training, labels come from the simulation (additive vs recessive). For empirical-only applications, use the pre-trained model from our training data.

3. **Organize inputs:** One CSV/Excel file per split (e.g., `train_HGDP.csv`, `test_HGDP.csv`) with columns: `dominance` plus the 22 features.

4. **Run training:**
   ```bash
   python training/train_xgboost.py --train train_HGDP.csv --test test_HGDP.csv --output models_hgdp --features 22
   ```

5. **Run empirical prediction** on new populations:
   ```bash
   python empirical_prediction/predict_empirical.py --empirical hgdp_pop_stats.txt --model-dir models_hgdp --output predictions/HGDP_pop.csv
   ```

---

## Simulations

Simulations are run via `slim/sim2vcf2stats.py`. Paths are configured in `config/paths.py` (or overridden with `--dir-master`, `--path-slim`).

```bash
# Example: one batch, one segment, additive dominance
python slim/sim2vcf2stats.py -s 1 -g 1 -d 1

# dominance: 1=additive, 2=partial, 3=recessive
# Output: slim/stat_1kseg_byexon/{segment}_{dominance}_rep{N}_1MB_stat.csv
```

See `slim/README.md` for details.

---

## Output interpretation

- **prediction_probability**: P(recessive) for each 1MB window
- **predicted_recessive**: 1 if prob ≥ 0.5, else 0
- ~9% of genome predicted recessive overall; ~3% in high exon-density (≥600/5MB) where model performs best

---

## Publication reporting: intermediate analyses

**METHODS_SUPPLEMENT.md** catalogues all intermediate analysis functions for Methods/Supplementary Materials, including:

- Data preprocessing (missing values, train-test split, normalization)
- Feature selection (KL divergence, feature importance, SelectKBest/RFE)
- Model comparison (22 vs 20 features, ancestry-only, 13 ML models)
- Stratification (exon density, power/FPR/specificity)
- ROC, precision-recall, confusion matrix analyses

**analysis_functions/analysis_functions.py** provides runnable implementations of these functions for reproducibility.

---
## Google Drive Supplementary Data

- Training data used in DominL training, and empirical predictions on 7 1KG populations can be found: https://drive.google.com/drive/u/2/folders/1sOg0dvZ4KUqurTOFA8gsKldPEJiKBEf8
