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
IN_PATH  = r"E:\Guido\sibel\Masters_Thesis\data\Wheat_Inventory30102025_Orthogroups_intensity_MEDIAN_NORM.xlsx"
SHEET    = "Wheat_Inventory30102025_Orthogr"

OUT_SIMPLE         = r"E:\Guido\sibel\Masters_Thesis\data\testall_cluster31102025.png"
OUT_CLUSTER        = r"E:\Guido\sibel\Masters_Thesis\data\testall_cluster31102025_1_.png"
OUT_CLUSTER_FAMILY = r"E:\Guido\sibel\Masters_Thesis\data\testall_zscore_31102025.png"

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

# ---------- whitelist ----------
SPECIES_ALLOW = {
    'GCF_904849725.1_Hordeum_vulgare',
    'GCF_902167145.1_Zea_mays',
    'GCA_963924085.1_Cenchrus_americanus.helixer',
    'GCF_034140825.1_Oryza_sativa',
    'GCF_000003195.3_Sorghum_bicolor',
    'GCF_018294505.1_Triticum_aestivum',
}

# ---------- main ----------

# load
df = pd.read_excel(IN_PATH, sheet_name=SHEET)
print(f"[chk] loaded: shape={df.shape}, ncols={len(df.columns)}")

# species cleanup + whitelist
df['Species'] = df['Species'].astype(str).str.strip()
species_values = [sp for sp in df['Species'].unique() if sp in SPECIES_ALLOW]
unknown = sorted(set(df['Species'].unique()) - SPECIES_ALLOW)
if unknown:
    print(f"[warn] {len(unknown)} species in file but not whitelisted (skipped): {unknown}")

print(f"[chk] will normalize across {len(species_values)} species (whitelisted)")

# intensity columns
META_EXACT  = {"HOG", "orthogroup", "Species", "ProteinID"}
META_PREFIX = ("Entry_",)
intensity_cols = [
    c for c in df.columns
    if ("MaxLFQ" in c) and (c not in META_EXACT) and (not any(c.startswith(p) for p in META_PREFIX))
]
if not intensity_cols:
    raise ValueError("No intensity columns matched '*MaxLFQ*'.")

print(f"[chk] intensity_cols={len(intensity_cols)}")
print("[chk] first 10 intensity cols:", intensity_cols[:10])
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
    # If empty or mostly non-letters (e.g., numeric IDs), return placeholder
    if not s or not re.search(r"[a-z]", s):
        return "—"
    for rx, fam in FAMILY_RULES:
        if re.search(rx, s):
            return fam
    toks = [t for t in re.split(r"[^a-z0-9]+", s) if t and t not in STOPWORDS]
    return "—" if not toks else " ".join(toks[:2])

    for rx, fam in FAMILY_RULES:
        if re.search(rx, s):
            return fam
    toks = [t for t in re.split(r"[^a-z0-9]+", s) if t and t not in STOPWORDS]
    return "—" if not toks else " ".join(toks[:2])
# ========= ADD-ON: split by species → per-species DF + log2/Z + PNG =========
from matplotlib.colors import TwoSlopeNorm

# 0) Where to write outputs
PER_SPECIES_DIR = Path(OUT_CLUSTER).parent / "per_species_pages"
PER_SPECIES_DIR.mkdir(parents=True, exist_ok=True)
SAVE_SPECIES_TABLES = True  # set False if you don't want TSVs for each species

def _safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s))

# 1) Pick the correct name column ONCE, then reuse for every species DF
def _norm_header(s: str) -> str:
    return re.sub(r"\s+", " ", s.replace("\u00A0"," ")).strip().lower()

def pick_name_col(df_any: pd.DataFrame) -> str | None:
    targets = {
        "proteinname triticum aestivum",
        "proteinname_triticum_aestivum",
        "protein name triticum aestivum",
    }
    for c in df_any.columns:
        nh = _norm_header(c).replace("_", " ")
        if nh in targets:
            return c
    # no strict match → don’t guess; return None
    return None

NAME_COL = pick_name_col(df)
if NAME_COL:
    print(f"[info] Family mapping will use column: {NAME_COL}")
else:
    print("[warn] No 'ProteinName_Triticum_aestivum' column found; families will be '—'.")

# 2) Build species list (keep your whitelist if you’re using one)
species_values = df["Species"].astype(str).str.strip().unique().tolist()
try:
    SPECIES_ALLOW
except NameError:
    SPECIES_ALLOW = None
if SPECIES_ALLOW:
    species_values = [sp for sp in species_values if sp in SPECIES_ALLOW]

print(f"[info] Creating per-species pages for {len(species_values)} species.")

# 3) Helpers for matrices
def make_log2_and_z(df_in: pd.DataFrame, cols: list[str]):
    X = df_in[cols].apply(pd.to_numeric, errors="coerce").to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        M = np.log2(X)
    M[~np.isfinite(M)] = np.nan
    mu = np.nanmean(M, axis=1, keepdims=True)
    sd = np.nanstd(M,  axis=1, keepdims=True)
    sd = np.where((~np.isfinite(sd)) | (sd == 0), 1.0, sd)
    Z = (M - mu) / sd
    return X, M, Z

# 4) Add family mapping to a DF slice
def add_family_column(df_slice: pd.DataFrame) -> pd.DataFrame:
    out = df_slice.copy()
    family_col = "ProteinFamily_Triticum_aestivum"
    if NAME_COL is None:
        out[family_col] = "—"
        return out
    main_names = (
        out[NAME_COL].fillna("").astype(str)
          .str.split(";", n=1).str[0].str.strip()
    )
    out[family_col] = main_names.map(name_to_family)
    return out

# 5) Plot helpers (PNG only)
def plot_heatmap_png(matrix_masked, cols, ylabels, title, cbar_label, out_path):
    w = max(8, len(cols) * 0.35)
    fig, ax = plt.subplots(figsize=(w, 8))
    cmap = plt.cm.RdBu_r; cmap.set_bad(color="lightgray")
    im = ax.imshow(matrix_masked, aspect="auto", cmap=cmap)
    fig.colorbar(im, ax=ax, label=cbar_label)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=90, ha="center", fontsize=8)
    ax.set_yticks(range(len(ylabels))); ax.set_yticklabels(ylabels, fontsize=7)
    ax.set_xlabel("Tissues")
    ax.set_ylabel("Hierarchical Orthogroup (HOG)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] saved PNG → {out_path}")

def plot_z_png(Z, M, cols, ylabels, title, out_path):
    w = max(8, len(cols) * 0.35)
    fig, ax = plt.subplots(figsize=(w, 8))
    cmap = plt.cm.RdBu_r; cmap.set_bad(color="lightgray")
    im = ax.imshow(np.ma.masked_where(np.isnan(M), Z), aspect="auto", cmap=cmap,
                   norm=TwoSlopeNorm(vmin=-2, vcenter=0.0, vmax=2))
    fig.colorbar(im, ax=ax, label="row-wise z-score within species (log2)")
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=90, ha="center", fontsize=8)
    ax.set_yticks(range(len(ylabels))); ax.set_yticklabels(ylabels, fontsize=7)
    ax.set_xlabel("Tissues")
    ax.set_ylabel("Hierarchical Orthogroup (HOG)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] saved PNG → {out_path}")

# 6) Drive: one DF + two PNGs per species (log2 & Z), plus optional TSV
for sp in species_values:
    df_sp = df.loc[df["Species"].astype(str).str.strip() == sp]
    if df_sp.empty:
        print(f"[warn] {sp}: no rows; skipping")
        continue

    # add family mapping to the species DF
    df_sp = add_family_column(df_sp)

    # HOG labels (or row indices if missing)
    if "HOG" in df_sp.columns:
        ylabels = df_sp["HOG"].astype(str).tolist()
    else:
        ylabels = [str(i) for i in range(len(df_sp))]

    # Using the SAME intensity_cols you defined earlier
    X, M, Z = make_log2_and_z(df_sp, intensity_cols)

    # Mask NaNs only (keeps zeros visible unless they became NaN after log2)
    M_masked = np.ma.masked_where(np.isnan(M), M)
    Z_masked = np.ma.masked_where(np.isnan(M), Z)

    # Save PNGs
    base = PER_SPECIES_DIR / _safe_name(sp)
    plot_heatmap_png(
        matrix_masked=M_masked,
        cols=intensity_cols,
        ylabels=ylabels,
        title=f"{sp} — Storage Protein Abundance across Tissues",
        cbar_label="log2 Intensity",
        out_path=str(base) + "_log2.png"
    )
    plot_z_png(
        Z=Z_masked,
        M=M,
        cols=intensity_cols,
        ylabels=ylabels,
        title=f"{sp} — Per-species row z-score",
        out_path=str(base) + "_z.png"
    )

    # Optional: save the per-species table (with family column)
    if SAVE_SPECIES_TABLES:
        out_tsv = str(base) + ".tsv"
        df_sp.to_csv(out_tsv, sep="\t", index=False)
        print(f"[OK] wrote per-species table → {out_tsv}")


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

ax1.set_xlabel("Tissues")
ax1.set_ylabel("Hierarchical Orthogroup (HOG)")
ax1.set_title("Storage Protein Abundance (log2 intensity)")
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

ax2.set_xlabel("Tissues")
ax2.set_ylabel("Hierarchical Orthogroup (HOG)")
ax2.set_title("Clustering Orthogroups (per-species row-wise z-score)")
fig2.tight_layout()
ensure_dir(OUT_CLUSTER)
fig2.savefig(OUT_CLUSTER, dpi=300, bbox_inches="tight")
plt.close(fig2)
print(f"[OK] saved {OUT_CLUSTER}")

# ========= PLOT 3: per-species z-score + family strip on far right =========
# reuse Z_masked and df_ordered computed above

# derive family from FIRST name before ';' in ProteinName_Triticum_aestivum
def _norm_header(s: str) -> str:
    return re.sub(r"\s+", " ", s.replace("\u00A0"," ")).strip().lower()

target_headers = {
    "proteinname triticum aestivum",
    "proteinname_triticum_aestivum",
    "protein name triticum aestivum",
}
name_col = None
for c in df_ordered.columns:
    nh = _norm_header(c).replace("_", " ")
    if nh in target_headers:
        name_col = c
        break

if name_col is None:
    # strict fallback OFF: we won't guess — avoid wrong mapping
    print("[warn] Could not find 'ProteinName_Triticum_aestivum' column; setting families to '—'.")
    df_ordered["ProteinFamily_Triticum_aestivum"] = "—"
else:
    print(f"[info] Using name column for family mapping: {name_col}")
    # show a few samples to verify it's really names
    print("[info] Sample name values:", df_ordered[name_col].astype(str).fillna("").head(5).tolist())

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
# ========= ADD-ON: per-species PNG pages (HOG on y-axis) =========
# ========= ADD-ON: split by species → per-species DF + log2/Z + PNG + TSV =========
from matplotlib.colors import TwoSlopeNorm

# output folder for per-species assets
PER_SPECIES_DIR = Path(OUT_CLUSTER).parent / "per_species_pages"
PER_SPECIES_DIR.mkdir(parents=True, exist_ok=True)
SAVE_SPECIES_TABLES = True  # also write the raw per-species table with family column

def _safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s))

# 1) pick exact name column once (for family mapping)
def _norm_header(s: str) -> str:
    return re.sub(r"\s+", " ", s.replace("\u00A0"," ")).strip().lower()

def pick_name_col(df_any: pd.DataFrame) -> str | None:
    targets = {
        "proteinname triticum aestivum",
        "proteinname_triticum_aestivum",
        "protein name triticum aestivum",
    }
    for c in df_any.columns:
        nh = _norm_header(c).replace("_", " ")
        if nh in targets:
            return c
    return None  # no strict match → do not guess

NAME_COL = pick_name_col(df)
print(f"[info] Family mapping column: {NAME_COL if NAME_COL else '(none → families set to —)'}")

# 2) species list (honor whitelist if present)
species_values = df["Species"].astype(str).str.strip().unique().tolist()
try:
    SPECIES_ALLOW
except NameError:
    SPECIES_ALLOW = None
if SPECIES_ALLOW:
    species_values = [sp for sp in species_values if sp in SPECIES_ALLOW]
print(f"[info] Creating per-species assets for {len(species_values)} species.")

# 3) helpers
def make_log2_and_z(df_in: pd.DataFrame, cols: list[str]):
    X = df_in[cols].apply(pd.to_numeric, errors="coerce").to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        M = np.log2(X)
    M[~np.isfinite(M)] = np.nan
    mu = np.nanmean(M, axis=1, keepdims=True)
    sd = np.nanstd(M,  axis=1, keepdims=True)
    sd = np.where((~np.isfinite(sd)) | (sd == 0), 1.0, sd)
    Z = (M - mu) / sd
    return X, M, Z

def add_family_column(df_slice: pd.DataFrame) -> pd.DataFrame:
    out = df_slice.copy()
    fam_col = "ProteinFamily_Triticum_aestivum"
    if NAME_COL is None:
        out[fam_col] = "—"
        return out
    main_names = (
        out[NAME_COL].fillna("").astype(str)
          .str.split(";", n=1).str[0].str.strip()
    )
    out[fam_col] = main_names.map(name_to_family)
    return out

def _style_axes_white_grid(ax, nrows: int, ncols: int, grid_w=0.5):
    # white background and white “cell lines”
    ax.set_facecolor("white")
    for sp in ax.spines.values():
        sp.set_color("white")
        sp.set_linewidth(0.0)
    ax.set_xticks(np.arange(-.5, ncols, 1), minor=True)
    ax.set_yticks(np.arange(-.5, nrows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=grid_w)
    ax.tick_params(which="minor", bottom=False, left=False)

def plot_heatmap_png(matrix_masked, cols, ylabels, title, cbar_label, out_path):
    nrows, ncols = matrix_masked.shape
    w = max(8, ncols * 0.35)
    fig, ax = plt.subplots(figsize=(w, 8))
    cmap = plt.cm.RdBu_r; cmap.set_bad(color="lightgray")
    im = ax.imshow(matrix_masked, aspect="auto", cmap=cmap)
    fig.colorbar(im, ax=ax, label=cbar_label)
    _style_axes_white_grid(ax, nrows, ncols)

    ax.set_xticks(range(ncols)); ax.set_xticklabels(cols, rotation=90, ha="center", fontsize=8)
    ax.set_yticks(range(nrows)); ax.set_yticklabels(ylabels, fontsize=7)
    ax.set_xlabel("Channels / Time Points")
    ax.set_ylabel("Orthogroup (HOG)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] saved PNG → {out_path}")

def plot_z_png(Z, M, cols, ylabels, title, out_path):
    nrows, ncols = Z.shape
    w = max(8, ncols * 0.35)
    fig, ax = plt.subplots(figsize=(w, 8))
    cmap = plt.cm.RdBu_r; cmap.set_bad(color="lightgray")
    im = ax.imshow(np.ma.masked_where(np.isnan(M), Z), aspect="auto", cmap=cmap,
                   norm=TwoSlopeNorm(vmin=-2, vcenter=0.0, vmax=2))
    fig.colorbar(im, ax=ax, label="row-wise z-score within species (log2)")
    _style_axes_white_grid(ax, nrows, ncols)

    ax.set_xticks(range(ncols)); ax.set_xticklabels(cols, rotation=90, ha="center", fontsize=8)
    ax.set_yticks(range(nrows)); ax.set_yticklabels(ylabels, fontsize=7)
    ax.set_xlabel("Channels / Time Points")
    ax.set_ylabel("Orthogroup (HOG)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] saved PNG → {out_path}")

# helper to derive a species token ("Genus_species") from df["Species"] value
SPECIES_TOKEN_RX = re.compile(r"([A-Z][a-z]+_[a-z]+)")
def species_token_from_value(sp_val: str) -> str | None:
    m = SPECIES_TOKEN_RX.search(str(sp_val))
    return m.group(1) if m else None

# 4) drive per species
for sp in species_values:
    df_sp = df.loc[df["Species"].astype(str).str.strip() == sp]
    if df_sp.empty:
        print(f"[warn] {sp}: no rows; skipping")
        continue

    # add family mapping
    df_sp = add_family_column(df_sp)

    # labels
    ylabels = df_sp["HOG"].astype(str).tolist() if "HOG" in df_sp.columns else [str(i) for i in range(len(df_sp))]

    # --- A) Use ALL intensity_cols (your original request)
    X_all, M_all, Z_all = make_log2_and_z(df_sp, intensity_cols)
    M_all_masked = np.ma.masked_where(np.isnan(M_all), M_all)
    Z_all_masked = np.ma.masked_where(np.isnan(M_all), Z_all)

    # TSVs: after log2 and after z (ALL cols)
    base = PER_SPECIES_DIR / _safe_name(sp)
    df_log2_all = df_sp.copy()
    df_z_all    = df_sp.copy()
    df_log2_all.loc[:, intensity_cols] = pd.DataFrame(M_all, index=df_sp.index, columns=intensity_cols)
    df_z_all.loc[:,    intensity_cols] = pd.DataFrame(Z_all, index=df_sp.index, columns=intensity_cols)
    df_log2_all.to_csv(str(base) + "_log2_all.tsv", sep="\t", index=False)
    df_z_all.to_csv(   str(base) + "_z_all.tsv",    sep="\t", index=False)
    print(f"[OK] wrote TSVs → {base.name}_log2_all.tsv, {base.name}_z_all.tsv")

    # PNGs (ALL cols)
    plot_heatmap_png(
        matrix_masked=M_all_masked,
        cols=intensity_cols,
        ylabels=ylabels,
        title=f"{sp} — Ortho Group × Tissue (log2, all MaxLFQ cols)",
        cbar_label="log2 intensity",
        out_path=str(base) + "_log2_all.png"
    )
    plot_z_png(
        Z=Z_all_masked,
        M=M_all,
        cols=intensity_cols,
        ylabels=ylabels,
        title=f"{sp} — Per-species row z-score (all MaxLFQ cols)",
        out_path=str(base) + "_z_all.png"
    )

    # --- B) Use SPECIES-SPECIFIC intensity cols only
    token = species_token_from_value(sp)
    sp_cols = [c for c in intensity_cols if token and token in c]
    if len(sp_cols) == 0:
        print(f"[warn] {sp}: no species-specific columns matched by token={token!r}; skipping species-only plots/TSVs.")
    else:
        X_spc, M_spc, Z_spc = make_log2_and_z(df_sp, sp_cols)
        M_spc_masked = np.ma.masked_where(np.isnan(M_spc), M_spc)
        Z_spc_masked = np.ma.masked_where(np.isnan(M_spc), Z_spc)

        # TSVs: after log2 and after z (SPECIES cols)
        df_log2_spc = df_sp.copy()
        df_z_spc    = df_sp.copy()
        df_log2_spc.loc[:, sp_cols] = pd.DataFrame(M_spc, index=df_sp.index, columns=sp_cols)
        df_z_spc.loc[:,    sp_cols] = pd.DataFrame(Z_spc, index=df_sp.index, columns=sp_cols)
        df_log2_spc.to_csv(str(base) + "_log2_species.tsv", sep="\t", index=False)
        df_z_spc.to_csv(   str(base) + "_z_species.tsv",    sep="\t", index=False)
        print(f"[OK] wrote TSVs → {base.name}_log2_species.tsv, {base.name}_z_species.tsv")

        # PNGs (SPECIES cols)
        plot_heatmap_png(
            matrix_masked=M_spc_masked,
            cols=sp_cols,
            ylabels=ylabels,
            title=f"{sp} — Ortho Group × Tissue (log2, species-only cols)",
            cbar_label="log2 intensity",
            out_path=str(base) + "_log2_species.png"
        )
        plot_z_png(
            Z=Z_spc_masked,
            M=M_spc,
            cols=sp_cols,
            ylabels=ylabels,
            title=f"{sp} — Per-species row z-score (species-only cols)",
            out_path=str(base) + "_z_species.png"
        )

    # optional: also save the per-species raw table (with family column) for provenance
    if SAVE_SPECIES_TABLES:
        out_tsv = str(base) + ".tsv"
        df_sp.to_csv(out_tsv, sep="\t", index=False)
        print(f"[OK] wrote per-species table → {out_tsv}")
# ========= ADD-ON: per-species Z plots using only non-empty columns =========
from matplotlib.colors import TwoSlopeNorm

# output folder
PER_SPECIES_DIR = Path(OUT_CLUSTER).parent / "per_species_z_subset"
PER_SPECIES_DIR.mkdir(parents=True, exist_ok=True)

def _safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s))

def _style_axes_white_grid(ax, nrows: int, ncols: int, grid_w=0.5):
    ax.set_facecolor("white")
    for spn in ax.spines.values():
        spn.set_color("white")
        spn.set_linewidth(0.0)
    ax.set_xticks(np.arange(-.5, ncols, 1), minor=True)
    ax.set_yticks(np.arange(-.5, nrows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=grid_w)
    ax.tick_params(which="minor", bottom=False, left=False)

def _rowwise_z(M: np.ndarray) -> np.ndarray:
    M2 = M.copy()
    M2[~np.isfinite(M2)] = np.nan
    mu = np.nanmean(M2, axis=1, keepdims=True)
    sd = np.nanstd(M2,  axis=1, keepdims=True)
    sd = np.where((~np.isfinite(sd)) | (sd == 0), 1.0, sd)
    return (M2 - mu) / sd

# species list (respect SPECIES_ALLOW if you already use one)
_species_vals = df["Species"].astype(str).str.strip().unique().tolist()
try:
    SPECIES_ALLOW
except NameError:
    SPECIES_ALLOW = None
if SPECIES_ALLOW:
    _species_vals = [sp for sp in _species_vals if sp in SPECIES_ALLOW]

for sp in _species_vals:
    mask_sp = (df["Species"].astype(str).str.strip() == sp)
    df_sp = df.loc[mask_sp]
    if df_sp.empty:
        print(f"[warn] {sp}: no rows; skipping")
        continue

    # pick only columns that have at least one numeric, non-NaN value for this species
    keep_cols = []
    for c in intensity_cols:
        vals = pd.to_numeric(df_sp[c], errors="coerce")
        if vals.notna().any():
            keep_cols.append(c)

    if not keep_cols:
        print(f"[warn] {sp}: no non-empty MaxLFQ columns; skipping")
        continue

    # log2 → z within kept columns only
    X = df_sp[keep_cols].apply(pd.to_numeric, errors="coerce").to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        M = np.log2(X)
    Z = _rowwise_z(M)

    # mask NaNs (keep zeros visible)
    Z_masked = np.ma.masked_where(~np.isfinite(M), Z)

    # y labels (HOG if present)
    ylabels = df_sp["HOG"].astype(str).tolist() if "HOG" in df_sp.columns else [str(i) for i in range(len(df_sp))]

    # plot
    nrows, ncols = Z_masked.shape
    w = max(8, ncols * 0.35)
    fig, ax = plt.subplots(figsize=(w, 8))
    cmap = plt.cm.RdBu_r; cmap.set_bad(color="lightgray")
    im = ax.imshow(Z_masked, aspect="auto", cmap=cmap, norm=TwoSlopeNorm(vmin=-2, vcenter=0.0, vmax=2))
    fig.colorbar(im, ax=ax, label="row-wise z-score within species (log2)")

    _style_axes_white_grid(ax, nrows, ncols)

    ax.set_xticks(range(ncols))
    ax.set_xticklabels(keep_cols, rotation=90, ha="center", fontsize=8)
    ax.set_yticks(range(nrows))
    ax.set_yticklabels(ylabels, fontsize=7)

    ax.set_xlabel("Tissues")
    ax.set_ylabel("Hierarchical Orthogroup (HOG)")
    ax.set_title(f"{sp} — Storage Protein Abundance across Tissues")

    out_png = PER_SPECIES_DIR / f"{_safe_name(sp)}_z_species_nonempty.png"
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] saved → {out_png}")


