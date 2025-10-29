#!/usr/bin/env python3
# family_plot_HOG_clean.py

import re
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as mpatches

# ========= CONFIG =========
IN_PATH  = r"E:\Guido\sibel\Masters_Thesis\data\All23102025_Orthogroups_intensity_mean.MEDIAN_NORM.xlsx"
SHEET    = "Sheet2"

OUT_SIMPLE         = r"E:\Guido\sibel\Masters_Thesis\data\all_cluster23102025.png"
OUT_CLUSTER        = r"E:\Guido\sibel\Masters_Thesis\data\all_cluster23102025_1_.png"
OUT_CLUSTER_FAMILY = r"E:\Guido\sibel\Masters_Thesis\data\all_zscore_specieswise23102025.png"

# ========= LOAD =========
df = pd.read_excel(IN_PATH, sheet_name=SHEET).fillna("")
# normalize headers
df.columns = (
    df.columns.astype(str)
      .str.replace("\u00A0", " ", regex=False)
      .str.replace(r"\s+", " ", regex=True)
      .str.strip()
)

# ========= INTENSITY COLUMNS (any 'P' + digits) =========
# --- find intensity columns: any column that has 'P' followed by digits ---
META_EXACT  = {"hog", "orthogroup"}   # exclude exact metadata cols
META_PREFIX = ("entry_",)             # exclude prefixes like Entry_*

# don't use \b; underscores break it. This matches P + digits,
# and stops before a non-digit (underscore, letter) or end-of-string.
P_CH_PATTERN = re.compile(r"(?i)P\d+(?=[^0-9]|$)")

def is_intensity_col(col: str) -> bool:
    c = col.strip()
    cl = c.lower()
    if cl in META_EXACT:
        return False
    if any(cl.startswith(p) for p in META_PREFIX):
        return False
    return bool(P_CH_PATTERN.search(c))

intensity_cols = [c for c in df.columns if is_intensity_col(c)]
if not intensity_cols:
    print("[DEBUG] No intensity columns found. First 40 cols:")
    for c in df.columns[:40]:
        print(" -", c)
    raise ValueError("No intensity columns matched (must contain 'P' followed by digits).")
else:
    print(f"[INFO] Matched {len(intensity_cols)} intensity columns (showing first 10):")
    for c in intensity_cols[:10]:
        print(" -", c)
from collections import defaultdict

SPECIES_TOKEN = re.compile(r"([A-Z][a-z]+_[a-z]+)")

def species_key(col: str) -> str:
    hits = SPECIES_TOKEN.findall(col)
    return hits[-1] if hits else ""

species_to_idx = defaultdict(list)
for j, c in enumerate(intensity_cols):
    sp = species_key(c)
    if sp:
        species_to_idx[sp].append(j)

# ========= HELPERS =========
def ensure_dir(path_str):
    Path(path_str).parent.mkdir(parents=True, exist_ok=True)

def log2_matrix_from_df(df_in, cols):
    X = df_in[cols].apply(pd.to_numeric, errors="coerce").to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        M = np.log2(X)
    M_clean = M.copy()
    M_clean[~np.isfinite(M_clean)] = np.nan
    return X, M, M_clean

def cluster_row_order(Z):
    Z_filled = np.nan_to_num(Z, nan=0.0, posinf=0.0, neginf=0.0)
    try:
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import pdist
        L = linkage(pdist(Z_filled, metric="euclidean"), method="average")
        return leaves_list(L)
    except Exception:
        C = np.corrcoef(Z_filled)
        C = np.nan_to_num((C + C.T) / 2, nan=0.0)
        w, v = np.linalg.eigh(C)
        return np.argsort(v[:, -1])

# --- family extraction helpers ---
STOPWORDS = {
    "protein","putative","like","fragment","chain","isoform","precursor",
    "uncharacterized","probable","possible","subunit","domain","containing"
}
FAMILY_RULES = [
    (r"dimeric\s+alpha[- ]?amylase\s+inhibitor", "alpha-amylase inhibitor (dimeric)"),
    (r"alpha[- ]?amylase/trypsin\s+inhibitor\s+cm\d+", "ATI CM"),
    (r"alpha[- ]?amylase/trypsin\s+inhibitor", "ATI"),
    (r"\bcm\d+\b|\bcm\b", "ATI CM"),
    (r"\bkunitz\b", "Kunitz inhibitor"),
    (r"\bserpin\b", "Serpin"),
    (r"\bgamma[- ]?gliadin\b", "γ-gliadin"),
    (r"\bdelta[- ]?gliadin\b", "δ-gliadin"),
    (r"\bgliadin\b", "gliadin"),
    (r"glutenin.*high\s+molecular\s+weight|h?mw[-\s]?gs", "HMW glutenin"),
    (r"glutenin.*low\s+molecular\s+weight|lmw[-\s]?gs", "LMW glutenin"),
    (r"\bglutenin\b", "glutenin"),
    (r"\bglobulin[- ]?3a\b", "Globulin-3A"),
    (r"\bglobulin[- ]?1\b", "Globulin-1"),
    (r"\bglobulin\b", "Globulin"),
    (r"\boleosin\b", "Oleosin"),
    (r"\bseipin\b", "Seipin"),
    (r"hydrophobic seed protein", "Hydrophobic seed protein"),
    (r"\bferritin\b", "Ferritin"),
    (r"\bgrain softness protein|\bgsp(-| )?1\b|puroindoline", "Grain Softness (Puroindoline)"),
    (r"\bavenin[- ]?like\b", "Avenin-like"),
    (r"cupin type-?1", "Cupin type-1"),
    (r"\bprolamin\b", "Prolamin"),
    (r"protease inhibitor", "Protease inhibitor-like"),
]
def _norm(s: str) -> str:
    s = (s or "").lower()
    s = s.replace("α","alpha").replace("γ","gamma").replace("δ","delta")
    s = re.sub(r"[()/_-]+", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s
def name_to_family(name: str) -> str:
    s = _norm(name)
    for rx, fam in FAMILY_RULES:
        if re.search(rx, s):
            return fam
    toks = [t for t in re.split(r"[^a-z0-9]+", s) if t and t not in STOPWORDS]
    return "—" if not toks else " ".join(toks[:2])

# ========= SPECIES GROUPING FOR PER-SPECIES Z =========
SPECIES_TOKEN = re.compile(r"([A-Z][a-z]+_[a-z]+)")  # e.g., Cenchrus_americanus
def species_key(col: str) -> str:
    hits = SPECIES_TOKEN.findall(col)
    return hits[-1] if hits else ""

species_to_idx = defaultdict(list)
for j, c in enumerate(intensity_cols):
    sp = species_key(c)
    if sp:
        species_to_idx[sp].append(j)

# ========= PLOT 1: simple log2 heatmap =========
_, M, _ = log2_matrix_from_df(df, intensity_cols)
M_masked = np.ma.masked_where((M == 0) | np.isnan(M), M)

fig1, ax1 = plt.subplots(figsize=(max(8, len(intensity_cols)*0.35), 8))
cmap = plt.cm.RdBu_r
cmap.set_bad(color="lightgray")
im1 = ax1.imshow(M_masked, aspect="auto", cmap=cmap)
fig1.colorbar(im1, ax=ax1, label="log2 intensity")

ax1.set_xticks(range(len(intensity_cols)))
ax1.set_xticklabels(intensity_cols, rotation=90, ha="center", fontsize=8)
if len(df) <= 60 and "HOG" in df.columns:
    ax1.set_yticks(range(len(df)))
    ax1.set_yticklabels(df["HOG"], fontsize=7)
else:
    ax1.set_yticks([])

ax1.set_xlabel("Time Points / Channels")
ax1.set_ylabel("Orthogroup")
ax1.set_title("Ortho Group × Tissue (log2)")
fig1.tight_layout()
ensure_dir(OUT_SIMPLE)
fig1.savefig(OUT_SIMPLE, dpi=300, bbox_inches="tight")
plt.close(fig1)
print(f"[OK] saved {OUT_SIMPLE}")

# ========= PLOT 2: clustered heatmap using per-species row-wise z-scores =========
X, M, M_clean = log2_matrix_from_df(df, intensity_cols)

# per-species row-wise z-score
Z = np.full_like(M_clean, np.nan)
for sp, idxs in species_to_idx.items():
    block = M_clean[:, idxs]
    mu = np.nanmean(block, axis=1, keepdims=True)
    sd = np.nanstd(block,  axis=1, keepdims=True)
    sd = np.where((~np.isfinite(sd)) | (sd == 0), 1.0, sd)
    Z[:, idxs] = (block - mu) / sd

row_order = cluster_row_order(np.nan_to_num(Z, nan=0.0))
Z_masked = np.ma.array(Z, mask=np.isnan(M_clean))[row_order, :]
df_ordered = df.iloc[row_order].reset_index(drop=True)

fig2, ax2 = plt.subplots(figsize=(max(8, len(intensity_cols)*0.35), 8))
norm = TwoSlopeNorm(vmin=-2, vcenter=0.0, vmax=2)
cmap = plt.cm.RdBu_r; cmap.set_bad(color="lightgray")
im2 = ax2.imshow(Z_masked, aspect="auto", cmap=cmap, norm=norm)
fig2.colorbar(im2, ax=ax2, label="per-species row z-score of log2(intensity)")

ax2.set_xticks(range(len(intensity_cols)))
ax2.set_xticklabels(intensity_cols, rotation=90, ha="center", fontsize=8)
if len(df_ordered) <= 60 and "HOG" in df_ordered.columns:
    ax2.set_yticks(range(len(df_ordered)))
    ax2.set_yticklabels(df_ordered["HOG"], fontsize=7)
else:
    ax2.set_yticks([])

ax2.set_xlabel("Time Points / Channels")
ax2.set_ylabel("Orthogroup")
ax2.set_title("Clustering Orthogroups (per-species row-wise z-score)")
fig2.tight_layout()
ensure_dir(OUT_CLUSTER)
fig2.savefig(OUT_CLUSTER, dpi=300, bbox_inches="tight")
plt.close(fig2)
print(f"[OK] saved {OUT_CLUSTER}")

# ========= PLOT 3: per-species z-score + family strip on far right =========
# reuse Z_masked and df_ordered computed above

# derive family from FIRST name before ';' in ProteinName_Triticum_aestivum
name_col = None
for c in df_ordered.columns:
    cc = c.replace("\u00A0"," ").strip()
    if re.search(r"^ProteinName[_\s]*Triticum[_\s]*aestivum$", cc, flags=re.I):
        name_col = c; break
if name_col is None:
    for c in df_ordered.columns:
        if re.search(r"protein.*name", c, re.I) and re.search(r"triticum.*aestivum", c, re.I):
            name_col = c; break
if name_col is None:
    raise KeyError("Could not find a 'ProteinName_Triticum_aestivum' column in the input file.")

main_names = (
    df_ordered[name_col].fillna("").astype(str)
      .str.split(";", n=1).str[0].str.strip()
)
family_col = "ProteinFamily_Triticum_aestivum"
df_ordered[family_col] = main_names.map(name_to_family)


# ---------------- LOCK FAMILY→COLOR MAP + DEBUG PRINTS ----------------
# Use the SAME labels that were used to build the row strip (df_ordered-aligned)
labels_series = (
    df_ordered[family_col]
      .fillna("—").replace("", "—")
      .astype(str)
)

# 1) LOCK category order (deterministic)
#    Option A: alphabetical, but put "—" at the end if present
cats = sorted(set(labels_series))
if "—" in cats:
    cats = [c for c in cats if c != "—"] + ["—"]  # move placeholder last

#    (If you prefer “first appearance” order, use the next line instead of the two above)
# cats = list(dict.fromkeys(labels_series))  # preserves first observed order

# 2) Build code map ONCE and reuse everywhere
cat_to_code = {c: i for i, c in enumerate(cats)}
codes = labels_series.map(cat_to_code).to_numpy()

# 3) Color LUT (stable, repeats tab20 if needed)
base = mpl.colormaps["tab20"].colors  # 20 distinct colors
reps = max(1, int(np.ceil(len(cats) / 20)))
colors_lut = np.vstack([base] * reps)[:len(cats)]
row_colors = colors_lut[codes]  # <- used for the strip

# 4) DEBUG: print locked map + counts to terminal
vc = labels_series.value_counts()
print("\n=== LOCKED family → code map (alphabetical; '—' last) ===")
for fam in cats:
    code = cat_to_code[fam]
    rgb  = tuple(np.round(colors_lut[code], 3))
    cnt  = int(vc.get(fam, 0))
    print(f"[{code:02d}] {fam:40s}  n={cnt:4d}  colorRGB={rgb}")

# (Optional) write to disk for provenance
_map_out = r"E:\Guido\sibel\Masters_Thesis\data\family_color_map.tsv"
pd.DataFrame({
    "family": cats,
    "code": [cat_to_code[c] for c in cats],
    "color_rgb": [",".join(map(str, np.round(colors_lut[cat_to_code[c]], 6))) for c in cats],
    "count": [int(vc.get(c, 0)) for c in cats],
}).to_csv(_map_out, sep="\t", index=False)
print(f"[info] wrote family→code map to {_map_out}")

# 5) DRAW heatmap (Z_masked already computed above)
fig3, ax3 = plt.subplots(figsize=(max(8, len(intensity_cols)*0.35), 8))
norm = TwoSlopeNorm(vmin=-2, vcenter=0.0, vmax=2)
cmap = plt.cm.RdBu_r; cmap.set_bad(color="lightgray")
im3 = ax3.imshow(Z_masked, aspect="auto", cmap=cmap, norm=norm)

ax3.set_xticks(range(len(intensity_cols)))
ax3.set_xticklabels(intensity_cols, rotation=90, ha="center", fontsize=8)
if len(df_ordered) <= 60 and "HOG" in df_ordered.columns:
    ax3.set_yticks(range(len(df_ordered)))
    ax3.set_yticklabels(df_ordered["HOG"], fontsize=7)
else:
    ax3.set_yticks([])

ax3.set_xlabel("Tissues")
ax3.set_ylabel("Orthogroup")
ax3.set_title("Clustering on intensities, per-species row z-score")

# --- RIGHT: thin heatmap colorbar + ultra-thin family strip ---
from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax3)

cax_cb = divider.append_axes("right", size="2.0%", pad=0.10)
fig3.colorbar(im3, cax=cax_cb, label="per-species row z-score of log2(intensity)")

cax_family = divider.append_axes("right", size="1.0%", pad=0.08)
cax_family.imshow(row_colors.reshape(-1, 1, 3), aspect="auto")
cax_family.set_xticks([]); cax_family.set_yticks([])
cax_family.set_ylim(ax3.get_ylim())

# 6) Legend built from the SAME LUT (so colors match the strip)
legend_labels = list(vc.index[:12])  # top 12 to avoid clutter
patches = [
    mpatches.Patch(facecolor=colors_lut[cat_to_code[label]], edgecolor='none', label=label)
    for label in legend_labels
]

# put legend a bit further right so it never overlaps the strip
ax3.legend(
    handles=patches,
    title=family_col,
    loc="upper left",
    bbox_to_anchor=(1.18, 1.0),
    frameon=False,
    fontsize=8,
    borderaxespad=0.0,
)

# Save
fig3.subplots_adjust(right=0.80)
fig3.tight_layout()
ensure_dir(OUT_CLUSTER_FAMILY)
fig3.savefig(OUT_CLUSTER_FAMILY, dpi=300, bbox_inches="tight")
plt.close(fig3)
print(f"[OK] saved {OUT_CLUSTER_FAMILY}")
