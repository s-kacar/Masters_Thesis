


#####################################################

#!/usr/bin/env python3
import re, numpy as np, pandas as pd
import matplotlib.pyplot as plt

# 1) load
df = pd.read_csv(r'E:\Guido\sibel\Masters_Thesis\Germination_Data\Germinationdata_Inventory_found.tsv',
                 sep='\t').fillna("")

# Normalize column names a bit (helps with stray spaces)
df.columns = (
    df.columns
      .str.replace('\u00A0', ' ', regex=False)  # non-breaking space
      .str.strip()
)

# Easiest: prefix check (fast, no regex)
intensity_cols = [c for c in df.columns if c.startswith("iBAQ_")]
if not intensity_cols:
    raise ValueError("No columns matched regex ^iBAQ_norm_ch_\\d+_[^_]+$. Example columns:\n"
                     + "\n".join(df.columns[:25]))

M = df[intensity_cols].apply(pd.to_numeric, errors="coerce").to_numpy()

# 3) optional log2
APPLY_LOG2 = True
if APPLY_LOG2:
    with np.errstate(divide="ignore", invalid="ignore"):
        M = np.log2(M)

# 4) mask NaNs and zeros
M_masked = np.ma.masked_where((M == 0) | np.isnan(M), M)

# 5) plot
plt.figure(figsize=(max(8, len(intensity_cols) * 0.35), 8))
cmap = plt.cm.RdBu_r
cmap.set_bad(color="lightgray")
im = plt.imshow(M_masked, aspect="auto", cmap=cmap)
plt.colorbar(im, label=("log2 intensity" if APPLY_LOG2 else "intensity"))

plt.xticks(range(len(intensity_cols)), intensity_cols, rotation=90, ha="center", fontsize=8)
if len(df) <= 60 and "orthogroup" in df.columns:
    plt.yticks(range(len(df)), df["orthogroup"], fontsize=7)
else:
    plt.yticks([])

plt.xlabel("Time Points")
plt.ylabel("orthogroup")
plt.title("Ortho Group × Tissue (6 Poaceae species)")
plt.tight_layout()

outpath = r"E:\Guido\sibel\Masters_Thesis\data\Plot_germination.png"
plt.savefig(outpath, dpi=300)
print(f"✓ saved heatmap to {outpath}")

# Show a window if you expect one; otherwise omit.
# plt.show()

#######################################################
#CLUSTER MAP
# --- 3) optional log2 (keep your code) ---
APPLY_LOG2 = True
if APPLY_LOG2:
    with np.errstate(divide="ignore", invalid="ignore"):
        M = np.log2(M)

# --- 3.5) optional row z-score before clustering (recommended) ---
ROW_ZSCORE = True
if ROW_ZSCORE:
    row_mu = np.nanmean(M, axis=1, keepdims=True)
    row_sd = np.nanstd(M, axis=1, keepdims=True)
    Z = (M - row_mu) / (row_sd + 1e-8)
else:
    Z = M.copy()

# --- 3.6) compute row order (cluster rows only) ---
def cluster_row_order(Z):
    # fill NaNs with 0 after z-scoring (0 ~ row mean)
    Z_filled = np.nan_to_num(Z, nan=0.0, posinf=0.0, neginf=0.0)
    try:
        # Prefer hierarchical clustering if SciPy is available
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import pdist
        # Euclidean on z-scored rows is standard
        L = linkage(pdist(Z_filled, metric="euclidean"), method="average")
        order = leaves_list(L)
        return order
    except Exception:
        # Fallback: spectral/correlation seriation (no SciPy)
        # Build correlation; turn NaNs to 0; use first eigenvector to sort
        C = np.corrcoef(Z_filled)
        C = np.nan_to_num(C, nan=0.0)
        # ensure symmetric
        C = (C + C.T) / 2
        # leading eigenvector of correlation (largest eigenvalue)
        w, v = np.linalg.eigh(C)
        lead = v[:, -1]
        order = np.argsort(lead)
        return order

row_order = cluster_row_order(Z)

# --- 4) mask NaNs and zeros AFTER computing the order ---
M_masked = np.ma.masked_where((M == 0) | np.isnan(M), M)
M_masked = M_masked[row_order, :]

# also reorder labels
df_ordered = df.iloc[row_order].reset_index(drop=True)

plt.figure(figsize=(max(8, len(intensity_cols) * 0.35), 8))
cmap = plt.cm.RdBu_r
cmap.set_bad(color="lightgray")

im = plt.imshow(M_masked, aspect="auto", cmap=cmap)
plt.colorbar(im, label=("log2 intensity" if APPLY_LOG2 else "intensity"))

plt.xticks(range(len(intensity_cols)), intensity_cols, rotation=90, ha="center", fontsize=8)

if len(df_ordered) <= 60 and "orthogroup" in df_ordered.columns:
    plt.yticks(range(len(df_ordered)), df_ordered["orthogroup"], fontsize=7)
else:
    plt.yticks([])

plt.xlabel("Time Points")
plt.ylabel("orthogroup")
plt.title("Clsutering Orthogroups Germination data")
plt.tight_layout()
plt.savefig(r"E:\Guido\sibel\Masters_Thesis\data\Clustered_Germination.png", dpi=300)

#####################################################################################
import numpy as np
from matplotlib.colors import TwoSlopeNorm
import re, numpy as np, pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(r'E:\Guido\sibel\Masters_Thesis\Germination_Data\AmylaseiGErminationa_Inventory_found.tsv',
                 sep='\t').fillna("")

# Normalize column names a bit (helps with stray spaces)
df.columns = (
    df.columns
      .str.replace('\u00A0', ' ', regex=False)  # non-breaking space
      .str.strip()
)

# Easiest: prefix check (fast, no regex)
intensity_cols = [c for c in df.columns if c.startswith("iBAQ_")]
if not intensity_cols:
    raise ValueError("No columns matched regex ^iBAQ_norm_ch_\\d+_[^_]+$. Example columns:\n"
                     + "\n".join(df.columns[:25]))

# Build raw numeric matrix
# --- build numeric, log2, column-zscore as you already do ---
# --- build numeric, log2 ---
X = df[intensity_cols].apply(pd.to_numeric, errors="coerce").to_numpy()
with np.errstate(divide="ignore", invalid="ignore"):
    M = np.log2(X)

# Make a finite-only copy for stats (turn -inf/+inf into NaN)
M_clean = M.copy()
M_clean[~np.isfinite(M_clean)] = np.nan

# --- column-wise z-score on the cleaned matrix ---
col_mu = np.nanmean(M_clean, axis=1, keepdims=True)
col_sd = np.nanstd(M_clean,  axis=1, keepdims=True)
col_sd = np.where((~np.isfinite(col_sd)) | (col_sd == 0), 1.0, col_sd)
Z = (M_clean - col_mu) / col_sd

# --- masking: mimic earlier behavior (mask only NaNs) ---
# (i.e., DON'T mask zeros; zeros -> -inf got turned to NaN in M_clean, so they’ll be masked here only if truly missing)
mask = np.isnan(M_clean)

# --- cluster rows on z-scored data (fill NaNs as 0 just for distances) ---
Z_filled = np.nan_to_num(Z, nan=0.0, posinf=0.0, neginf=0.0)
try:
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import pdist
    L = linkage(pdist(Z_filled, metric="euclidean"), method="average")
    row_order = leaves_list(L)
except Exception:
    C = np.corrcoef(Z_filled)
    C = np.nan_to_num((C + C.T) / 2, nan=0.0)
    w, v = np.linalg.eigh(C)
    row_order = np.argsort(v[:, -1])

# reorder
Z_masked = np.ma.array(Z, mask=mask)[row_order, :]
df_ordered = df.iloc[row_order].reset_index(drop=True)

# --- plot (unchanged) ---
from matplotlib.colors import TwoSlopeNorm
vmin, vmax = -2.0, 2.0
norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)

fig, ax = plt.subplots(figsize=(max(8, len(intensity_cols)*0.35), 8))
cmap = plt.cm.RdBu_r
cmap.set_bad(color="lightgray")
im = ax.imshow(Z_masked, aspect="auto", cmap=cmap, norm=norm)
fig.colorbar(im, ax=ax, label="column z-score of log2(iBAQ)")

ax.set_xticks(range(len(intensity_cols)))
ax.set_xticklabels(intensity_cols, rotation=90, ha="center", fontsize=8)

labels = df_ordered["orthogroup"].values if "orthogroup" in df_ordered.columns else np.arange(len(df_ordered))
max_labels = 60
n = len(labels)
if n <= max_labels:
    yt_pos, yt_lab = np.arange(n), labels
else:
    step = int(np.ceil(n / max_labels))
    yt_pos, yt_lab = np.arange(0, n, step), labels[::step]

ax.set_yticks(yt_pos)
ax.set_yticklabels(yt_lab, fontsize=7)

ax.set_xlabel("Timepoint for each species")
ax.set_ylabel("Orthogroup")
ax.set_title("Clustermap based on Orthogroups after Row-wise Z-score Normalization")

# give margins so labels don't get cut
fig.subplots_adjust(left=0.25, bottom=0.35)
plt.tight_layout()

plt.savefig(r"E:\Guido\sibel\Masters_Thesis\data\row_Zscore_Amylase.png", dpi=300)
