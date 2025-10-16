
#!/usr/bin/env python3
# mean_per_cell_ignore_multival.py
# Usage: python mean_per_cell_ignore_multival.py  input.tsv  output.tsv

import sys, re, math
import pandas as pd
from pathlib import Path
##### MEAN CALCULATION OF COLLAPSED FINAL TSV
if len(sys.argv) != 3:
    print("Plotting.py  input.tsv  output.tsv")
    sys.exit(0)

IN  = Path(sys.argv[1])
OUT = Path(sys.argv[2])

df = pd.read_csv(IN, sep="\t", dtype=str).fillna("")

# 1) pick intensity columns only
intensity_cols = [c for c in df.columns if re.search(r"_P\d", c) or "MaxLFQ" in c]

def process_cell(cell: str) -> str:
    if not isinstance(cell, str) or cell.strip() == "":
        return ""
    parts = [p.strip() for p in cell.split(",") if p.strip()]

    # single value → keep as-is
    if len(parts) == 1:
        return parts[0]

    # multiple values
    vals = []
    nonzero_vals = []
    for p in parts:
        try:
            x = float(p)
        except ValueError:
            continue
        vals.append(x)
        if x != 0.0 and not math.isnan(x):
            nonzero_vals.append(x)

    if nonzero_vals:
        # average only nonzero/non-NaN values
        return str(sum(nonzero_vals) / len(nonzero_vals))
    elif all(v == 0.0 for v in vals):
        # all zeros → return 0
        return "0"
    else:
        # only NaNs or unparsable values → blank
        return ""

# 3) apply to intensity columns only
for col in intensity_cols:
    df[col] = df[col].apply(process_cell)

# 4) save
df.to_csv(OUT, sep="\t", index=False)
print(f"✓ Wrote {OUT}  ({df.shape[0]} rows × {df.shape[1]} cols)")


#####################################################

#!/usr/bin/env python3
# plot_heatmap.py
# Usage: python plot_heatmap.py input.tsv heatmap.png

import sys, re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ── 1) load ───────────────────────────────────────────────────────────
df = pd.read_csv(r"E:\Guido\sibel\Masters_Thesis\data\Poacea_Wheat_Based_Orthogroups_Intensity_MEAN.tsv", sep="\t", dtype=str).fillna("")
# change intensity cols if named differently
intensity_cols = [
    c for c in df.columns 
    if "_" not in c and c.lower() != "hog"
]
M = df[intensity_cols].apply(pd.to_numeric, errors="coerce").to_numpy()

# OPTIONAL: log2
APPLY_LOG2 = True
if APPLY_LOG2:
    with np.errstate(divide="ignore", invalid="ignore"):
        M = np.log2(M)

# mask NaNs and zeros
M_masked = np.ma.masked_where((M == 0) | np.isnan(M), M)

plt.figure(figsize=(max(8, len(intensity_cols)*0.35), 8))

# set colormap so masked cells appear gray
cmap = plt.cm.RdBu_r
cmap.set_bad(color="lightgray")

im = plt.imshow(M_masked, aspect="auto", cmap=cmap)
plt.colorbar(im, label=("log2 intensity" if APPLY_LOG2 else "intensity"))

# X labels
plt.xticks(range(len(intensity_cols)), intensity_cols, rotation=90, ha="center", fontsize=8)

# Y labels (only if not too many HOGs)
if len(df) <= 60:
    plt.yticks(range(len(df)), df["HOG"], fontsize=7)
else:
    plt.yticks([])

plt.xlabel("Tissues (MaxLFQ)")
plt.ylabel("HOG")
plt.title("Ortho Group × Tissue (6 Poaceae species)")
plt.tight_layout()
plt.savefig(r"E:\Guido\sibel\Masters_Thesis\data\Plot.png", dpi=300)
print(f"✓ saved heatmap to {r"E:\Guido\sibel\Masters_Thesis\data\Plot.png"}")
#######################################################
#!/usr/bin/env python3
# grouped_heatmap.py
import re, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# ── paths ────────────────────────────────────────────────────────────
IN  = Path(r"E:\Guido\sibel\Masters_Thesis\data\Poacea_Wheat_Based_Orthogroups_Intensity_MEAN.tsv")
OUT = Path(r"E:\Guido\sibel\Masters_Thesis\data\Plot.png")
# ── load ─────────────────────────────────────────────────────────────
df = pd.read_csv(IN, sep="\t", dtype=str).fillna("")
if "HOG" not in df.columns:
    sys.exit("HOG column not found")
cols = list(df.columns)
# --- build species blocks from Entry_<species> columns ---
cols = list(df.columns)
entry_idx = [(i, c) for i, c in enumerate(cols) if c.startswith("Entry_")]
if not entry_idx:
    sys.exit("No 'Entry_<species>' columns found; cannot infer species blocks.")

is_tissue = lambda c: ("_" not in c) and (c.lower() != "hog")
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.patches import Rectangle
from matplotlib.colors import ListedColormap

# --- assume blocks, X_labels, X_species, block_spans, and M_raw (with triple NaN gap cols) are already built ---

# 1) MASKS FROM RAW (pre-log2)
gap_cols   = np.array([sp == "_gap_" for sp in X_species])     # shape (n_cols,)
gap_mask   = np.tile(gap_cols, (len(df), 1))                    # shape (n_rows, n_cols)

data_nan_mask = np.isnan(M_raw) & ~gap_mask                     # real missing values
zero_mask     = np.isfinite(M_raw) & (M_raw == 0) & ~gap_mask   # true zeros
valid_mask    = np.isfinite(M_raw) & (M_raw > 0) & ~gap_mask    # values to plot (pos)

# 2) LOG2 ONLY ON VALID CELLS
M = np.full(M_raw.shape, np.nan, dtype=float)
with np.errstate(divide="ignore", invalid="ignore"):
    M[valid_mask] = np.log2(M_raw[valid_mask])

# 3) SCALING WITH A REAL MIDPOINT (avoid “all red”)
if np.any(valid_mask):
    vmin    = np.nanquantile(M[valid_mask], 0.02)
    vmax    = np.nanquantile(M[valid_mask], 0.98)
    vcenter = np.nanmedian(M[valid_mask])   # use 0.0 if you mean-centered earlier
else:
    vmin, vmax, vcenter = 0.0, 1.0, 0.5
norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)

# 4) PLOT LAYERS: NaN(gray) < zeros(gray+hatched) < valid(RdBu)
fig = plt.figure(figsize=(max(10, len(X_labels)*0.32), 8))
ax = plt.gca()
ax.set_facecolor("white")  # gutters remain white

n_rows, n_cols = M.shape

# A) data-NaN as solid gray (but NOT gaps)
nan_layer = np.ma.masked_where(~data_nan_mask, np.ones_like(M))
ax.imshow(nan_layer, aspect="auto", interpolation="none",
          cmap=ListedColormap(["#c9c9c9"]), vmin=0, vmax=1, zorder=1)

# B) zeros as gray with stripes
zi, zj = np.where(zero_mask)
for r, c in zip(zi, zj):
    ax.add_patch(Rectangle((c-0.5, r-0.5), 1, 1,
                           facecolor="#c9c9c9", edgecolor="#888888",
                           hatch="///", linewidth=0.0, zorder=2))

# C) valid heatmap
M_valid = np.ma.masked_where(~valid_mask, M)
cmap = plt.cm.RdBu_r.copy()
im = ax.imshow(M_valid, aspect="auto", interpolation="none",
               cmap=cmap, norm=norm, zorder=3)
plt.colorbar(im, label="log2 intensity")

# ticks
plt.xticks(range(len(X_labels)), X_labels, rotation=90, ha="center", fontsize=8)
if len(df) <= 60:
    plt.yticks(range(len(df)), df["HOG"], fontsize=7)
else:
    plt.yticks([])

# species-colored tick labels (greys)
uniq = list(dict.fromkeys([sp for sp, _ in blocks]))
greys = ["#111","#2c2c2c","#474747","#636363","#7e7e7e","#999","#b3b3b3"]
sp_to_color = {sp: greys[i % len(greys)] for i, sp in enumerate(uniq)}
for tick, sp in zip(ax.get_xticklabels(), X_species):
    tick.set_color("lightgray" if sp == "_gap_" else sp_to_color.get(sp, "#555"))

# optional: light lines at gap edges + black frames around species
for i, sp in enumerate(X_species):
    if sp == "_gap_":
        ax.axvline(i - 0.5, color="#e6e6e6", lw=0.8, zorder=4)

for x0, x1 in block_spans:
    ax.add_patch(Rectangle((x0 - 0.5, -0.5), x1 - x0 + 1, n_rows,
                           fill=False, ec="black", lw=1.2, zorder=5))

plt.xlabel("Tissues (grouped by species)")
plt.ylabel("HOG")
plt.title("Orthogroup × Tissue — RdBu (median-centered), NaNs gray, zeros hatched, triple gutters")
plt.tight_layout()
plt.savefig(r"E:\Guido\sibel\Masters_Thesis\data\Plot.png", dpi=300)
print("✓ saved heatmap")
