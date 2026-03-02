#!/bin/bash
# DominL Integrated Pipeline - Master run script
# Reference: Zhang et al. (2025) bioRxiv doi:10.1101/2025.05.07.652751

set -e
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PIPELINE_DIR"

echo "=== DominL Integrated Pipeline ==="

# Step 1: Training (requires train_data_all.xlsx or train_data_all.csv)
if [ -n "$TRAIN_DATA" ] && [ -f "$TRAIN_DATA" ]; then
    echo "[1/3] Training XGBoost model..."
    python training/train_xgboost.py --train "$TRAIN_DATA" \
        ${TEST_DATA:+--test "$TEST_DATA"} \
        --output models
else
    echo "[1/3] Skipping training (set TRAIN_DATA, TEST_DATA for training)"
    if [ ! -d models ] || [ ! -f models/xgboost_model_22features.joblib ]; then
        echo "  WARNING: No trained model found. Run training first or copy pre-trained model to models/"
    fi
fi

# Step 2: Empirical prediction (requires empirical stats and trained model)
if [ -n "$EMPIRICAL_DATA" ] && [ -f "$EMPIRICAL_DATA" ] && [ -f models/xgboost_model_22features.joblib ]; then
    echo "[2/3] Running empirical prediction..."
    POP=$(basename "$EMPIRICAL_DATA" | sed 's/.*_\([A-Z]*\)_.*/\1/' || echo "results")
    python empirical_prediction/predict_empirical.py \
        --empirical "$EMPIRICAL_DATA" \
        --model-dir models \
        --output "predictions/${POP}_predictions.csv"
else
    echo "[2/3] Skipping empirical prediction (set EMPIRICAL_DATA and ensure model exists)"
fi

echo "[3/3] Done. See README.md for full pipeline (simulations, etc.)."
