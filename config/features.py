"""
DominL Pipeline - Feature configuration.
22 features selected based on KL divergence overlap between empirical and simulation data.
Reference: Zhang et al. (2025) bioRxiv doi:10.1101/2025.05.07.652751
"""

# 20 features (without recombination rate)
SELECTED_FEATURES_20 = [
    'exon_density', 'mean_introg_anc', 'exon_window', 'introg_anc_window',
    'divergence_p3_p1', 'divergence_p3_p2', 'watterson_theta_p3',
    'windowed_tajima_d_p3', 'D', 'Het', 'Q95', 'U0', 'U20', 'U50', 'U80',
    'num_seg_p1', 'num_seg_p3', 'num_private_seg_p1', 'num_private_seg_p2',
    'num_private_seg_p3'
]

# 22 features (adds mean_recrate, recrate_window)
SELECTED_FEATURES_22 = [
    'exon_density', 'mean_introg_anc', 'exon_window', 'introg_anc_window',
    'divergence_p3_p1', 'divergence_p3_p2', 'watterson_theta_p3',
    'windowed_tajima_d_p3', 'D', 'Het', 'Q95', 'U0', 'U20', 'U50', 'U80',
    'num_seg_p1', 'num_seg_p3', 'num_private_seg_p1', 'num_private_seg_p2',
    'num_private_seg_p3', 'mean_recrate', 'recrate_window'
]

# Columns to drop from data (metadata, not features)
COLUMNS_TO_DROP = ['segment', 'rep', 'start', 'end', 'patterson_f3']

# Mapping from empirical data column names to pipeline names
EMPIRICAL_COLUMN_MAPPING = {
    "nea_anc": "introg_anc_window",
    "exon_5mb": "exon_density",
    "exon": "exon_window",
    "nea_anc_5mb": "mean_introg_anc",
    "recomb_5mb": "mean_recrate",
    "recomb": "recrate_window"
}
