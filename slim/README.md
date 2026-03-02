# DominL Simulation Pipeline (SLiM → Stats)

This folder contains the simulation and feature-extraction steps of the DominL pipeline. It generates training data (summary statistics tables) that are then fed into the training pipeline.

## Workflow

```
Baseline SLiM script  →  Python modifies script (dominance, region, rep)
                              ↓
                    Run SLiM simulation
                              ↓
                    Tree sequence + VCF outputs
                              ↓
                    Extract features (ancestry, D/fD, diversity, etc.)
                              ↓
                    Write stats tables (CSV)
                              ↓
                    After enough replicates → training/
```

## Files

| File | Description |
|------|-------------|
| `sim2vcf2stats.py` | Main script: modifies SLiM template, runs simulations, extracts stats |
| `Gravel_Eurasian_varyh.slim` | Baseline SLiM script (Eurasian demography, dominance varies) |
| `stat2ML.py` | **Legacy** – original ML training; use `training/` instead |
| `1kseg_byexon/` | Segment info files (`sim_seq_info_seg*.txt`) with exon and recombination annotations |

## Setup

1. **Edit paths** in `../config/paths.py`:
   - `DIR_MASTER`: Simulation root (default: this `slim/` folder)
   - `PATH_SLIM`: Path to SLiM executable (e.g. `slim` or `/path/to/slim`)

2. **Dependencies**: `msprime`, `pyslim`, `scikit-allel`, `bcftools` (for VCF merging)

3. **SLiM**: Install from [Messer Lab](https://messerlab.org/slim/) or conda:
   ```bash
   conda install -c conda-forge slim
   ```

## Running

```bash
# From pipeline root
python slim/sim2vcf2stats.py -g 1 -d 1 -s 1

# With custom paths
python slim/sim2vcf2stats.py -g 1 -d 1 --dir-master /path/to/sim_root --path-slim /path/to/slim
```

**Arguments:**
- `-g, --gene`: Genomic region index (1–N, from `1kseg_byexon/sim_seq_info_seg*.txt`)
- `-d, --dominance`: 1=additive, 2=partial, 3=recessive
- `-s, --sim`: Batch ID (for parallel runs)
- `--dir-master`: Override simulation root
- `--path-slim`: Override SLiM executable path

Each run performs 100 replicates (configurable in script), writes stats to `stat_1kseg_byexon/{segment}_{dominance}_rep{N}_1MB_stat.csv`.

## Output Columns

Stats tables include: `segment`, `rep`, `start`, `end`, `dominance`, `exon_density`, `mean_recrate`, `mean_introg_anc`, `exon_window`, `recrate_window`, plus population genetics stats (D, fD, diversity, divergence, etc.). These feed into `training/train_xgboost.py` and `training/train_compare_models.py`.

## Population Labels

- **p1**: Archaic (Neanderthal)
- **p2**: African (YRI)
- **p3**: European (CEU)
