"""
DominL Pipeline - Path configuration.
Edit these paths for your environment.
"""

import os

# Base directory - parent of this pipeline folder
PIPELINE_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Simulation root: where segment files (1kseg_byexon), SLiM template, and outputs live.
# Default: slim/ folder in pipeline (contains 1kseg_byexon, Gravel_Eurasian_varyh.slim)
# For production, set to a separate dir (e.g. simulation_data) and copy segments there.
DIR_MASTER = os.path.join(PIPELINE_ROOT, 'slim')

# SLiM executable path (use full path if not in PATH)
PATH_SLIM = 'slim'

# Paths relative to DIR_MASTER
SEGMENTS_DIR = '1kseg_byexon'
STATS_OUTPUT_DIR = 'stat_1kseg_byexon'
