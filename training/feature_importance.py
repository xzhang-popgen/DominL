import pandas as pd
import xgboost as xgb
import matplotlib.pyplot as plt

df = pd.read_csv("/home/zhulx/lab/summary_statistics/data/trainingdata_features.csv")
print(df.head())
df.shape
data = df.values

features = [
    "num_seg_p2", "num_private_seg_p2", "exon_window",
    "num_variant_window", "mean_introg_anc", "watterson_theta_p3",
    "divergence_p3_p1", "exon_density", "num_private_seg_p3",
    "divergence_p3_p2", "num_seg_p1", "recrate_window", "mean_recrate",
    "df_p3_p1", "U50", "U0", "Het", "num_seg_p3", "garud_h1", "RD", "D",
    "hap_diversity_p3", "df_p3_p2", "Q95", "introg_anc_window",
    "U20", "fD", "num_private_seg_p1", "windowed_tajima_d_p3",
    "U80", "garud_h2_h1", "garud_h12"
]
features_cnn = [
    "num_seg_p2",
    "num_private_seg_p2",
    "exon_window",
    "num_variant_window",
    "mean_introg_anc",
    "watterson_theta_p3",
    "divergence_p3_p1",
    "exon_density",
    "num_private_seg_p3",
    "divergence_p3_p2",
    "num_seg_p1",
    "U50",
    "U0",
    "Het",
    "num_seg_p3",
    "D",
    "Q95",
    "introg_anc_window",
    "U20",
    "num_private_seg_p1",
    "windowed_tajima_d_p3",
    "U80"
]
features_other = [
    "num_private_seg_p2",
    "exon_window",
    "mean_introg_anc",
    "watterson_theta_p3",
    "divergence_p3_p1",
    "exon_density",
    "num_private_seg_p3",
    "divergence_p3_p2",
    "num_seg_p1",
    "U50",
    "U0",
    "Het",
    "num_seg_p3",
    "D",
    "Q95",
    "introg_anc_window",
    "U20",
    "num_private_seg_p1",
    "windowed_tajima_d_p3",
    "U80"
]


X = df[features]
y = df["dominance"]

# --- Train XGB model ---
model = xgb.XGBClassifier(
    # n_estimators=200,
    # max_depth=4,
    # learning_rate=0.1,
    random_state=42 #,
    # verbosity=1,
    # use_label_encoder=False,    
    # eval_metric="auc",      # or "error", "auc" metrics
)
model.fit(X, y)
model_cnn = xgb.XGBClassifier(random_state=42)
model_cnn.fit(df[features_cnn], y)
model_other = xgb.XGBClassifier(random_state=42)
model_other.fit(df[features_other], y)



# --- 3. extract and order feature importances ---
# model.feature_importances_ 
importances = pd.Series(model.feature_importances_, index=features)
importances = importances.sort_values(ascending=False)
importances_cnn = pd.Series(model_cnn.feature_importances_, index=features_cnn).sort_values(ascending=False)
importances_other = pd.Series(model_other.feature_importances_, index=features_other).sort_values(ascending=False)


print("=== Feature importances (descending) ===")
print(importances)

print("=== Feature importances for cnn featrues (descending) ===")
print(importances_cnn)

print("=== Feature importances for other models' features (descending) ===")
print(importances_other)

# --- 4. visualization ---
plt.figure(figsize=(12, 6))
importances.plot.bar()
plt.title("XGBoost Feature Importance")
plt.ylabel("Importance (gain)")
plt.tight_layout()
plt.show()


# # Log loss
# model = xgb.XGBClassifier(
#     n_estimators=200,
#     max_depth=4,
#     learning_rate=0.1,
#     random_state=42,
#     verbosity=1,
#     use_label_encoder=False,    
#     eval_metric="logloss",      # or "error", "auc" metrics
# )
# num_private_seg_p2      0.251872
# num_seg_p2              0.169422
# mean_introg_anc         0.086522
# exon_window             0.081105
# num_private_seg_p3      0.069366
# watterson_theta_p3      0.067818
# num_variant_window      0.060386
# divergence_p3_p1        0.032404
# num_seg_p1              0.028284
# exon_density            0.026295
# df_p3_p1                0.018336
# divergence_p3_p2        0.017644
# recrate_window          0.012117
# mean_recrate            0.011597
# U50                     0.010235
# U0                      0.010150
# Het                     0.008139
# D                       0.007118
# num_seg_p3              0.006056
# RD                      0.005819
# introg_anc_window       0.005312
# hap_diversity_p3        0.003214
# fD                      0.003022
# num_private_seg_p1      0.002895
# windowed_tajima_d_p3    0.002060
# Q95                     0.001694
# U20                     0.001119
# garud_h2_h1             0.000000
# garud_h1                0.000000
# df_p3_p2                0.000000
# U80                     0.000000
# garud_h12               0.000000

model = xgb.XGBClassifier(random_state=42)
num_private_seg_p2      0.183392
num_seg_p2              0.167184
exon_window             0.109873
num_variant_window      0.084310
watterson_theta_p3      0.080463
mean_introg_anc         0.061320
divergence_p3_p1        0.052295
num_private_seg_p3      0.037139
exon_density            0.029557
divergence_p3_p2        0.023804
num_seg_p1              0.022224
recrate_window          0.016700
mean_recrate            0.014689
df_p3_p1                0.012344
U0                      0.010175
garud_h1                0.009383
Het                     0.008242
U50                     0.007580
RD                      0.007134
D                       0.007120
hap_diversity_p3        0.007008
num_seg_p3              0.006580
fD                      0.006206
introg_anc_window       0.006195
U20                     0.005794
Q95                     0.005272
num_private_seg_p1      0.005138
windowed_tajima_d_p3    0.004820
df_p3_p2                0.004106
U80                     0.003950
garud_h2_h1             0.000000
garud_h12               0.000000


num_private_seg_p2      0.206274
num_seg_p2              0.179453
exon_window             0.109913
num_variant_window      0.094737
watterson_theta_p3      0.088442
mean_introg_anc         0.063915
divergence_p3_p1        0.048051
num_private_seg_p3      0.035945
exon_density            0.030501
divergence_p3_p2        0.025848
num_seg_p1              0.024043
U50                     0.011718
U0                      0.011571
num_private_seg_p1      0.010037
Het                     0.009900
num_seg_p3              0.008754
D                       0.007920
introg_anc_window       0.007226
U80                     0.006851
Q95                     0.006785
U20                     0.006687
windowed_tajima_d_p3    0.005429



=== Feature importances for other models' features (descending) ===
num_private_seg_p2      0.268178
watterson_theta_p3      0.170303
exon_window             0.144755
num_private_seg_p3      0.073946
mean_introg_anc         0.070597
U0                      0.045001
exon_density            0.039558
divergence_p3_p2        0.033227
num_seg_p3              0.024013
divergence_p3_p1        0.021283
U50                     0.015356
num_private_seg_p1      0.014529
num_seg_p1              0.011854
U80                     0.011336
introg_anc_window       0.010291
Het                     0.010245
Q95                     0.010027
U20                     0.009503
D                       0.008169
windowed_tajima_d_p3    0.007828
