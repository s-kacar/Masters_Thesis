#!/usr/bin/env python3
import re
import numpy as np
import pandas as pd
from pathlib import Path
from decimal import Decimal, ROUND_HALF_UP, getcontext

getcontext().prec = 28  # plenty for reporting

def round10(x: float) -> float:
    # for DISPLAY/EXPORT only; not used in math
    return float(Decimal(str(x)).quantize(Decimal("0.0000000000"), rounding=ROUND_HALF_UP))

# --- paths ---
IN_PATH  = r"E:\Guido\sibel\Masters_Thesis\data\melted_OG_matrix_non_comma_seperated_based_oninventory28.10.tsv"
OUT_PATH = r"E:\Guido\sibel\Masters_Thesis\data\All23102025_Orthogroups_intensity_28.10mean.MEDIAN_NORM.tsv"

# --- load ---
df = pd.read_csv(IN_PATH, sep="\t")
df.columns = (
    df.columns.astype(str)
      .str.replace("\u00A0", " ", regex=False)
      .str.replace(r"\s+", " ", regex=True)
      .str.strip()
)

# --- pick intensity columns: *_P<digits> only ---
META_EXACT  = {"hog", "orthogroup"}
META_PREFIX = ("entry_",)
P_CH_PATTERN = re.compile(r"(?i)[_-]p\d+(?=[_-]|$)")

def is_intensity_col(col: str) -> bool:
    cl = col.strip().lower()
    if cl in META_EXACT:
        return False
    if any(cl.startswith(p) for p in META_PREFIX):
        return False
    return bool(P_CH_PATTERN.search(col))

intensity_cols = [c for c in df.columns if is_intensity_col(c)]
if not intensity_cols:
    raise ValueError("No intensity columns matched '*_P<digits>'.")

# --- species key = Genus_species token in the column name ---
SPECIES_TOKEN = re.compile(r"([A-Z][a-z]+_[a-z]+)")  # e.g., Cenchrus_americanus

def species_key(col: str) -> str:
    hits = SPECIES_TOKEN.findall(col)
    return hits[-1] if hits else ""

# Set to None to auto-discover; or keep a set to enforce exactly these
SPECIES_ALLOW = {
    "Cenchrus_americanus",
    "Sorghum_bicolor",
    "Triticum_aestivum",
    "Hordeum_vulgare",
    "Zea_mays",
    "Oryza_sativa",
}
# SPECIES_ALLOW = None  # uncomment to disable enforcement

from collections import defaultdict
species_to_cols = defaultdict(list)
for c in intensity_cols:
    sp = species_key(c)
    if not sp:
        continue
    if SPECIES_ALLOW is None or sp in SPECIES_ALLOW:
        species_to_cols[sp].append(c)

# sanity check
for sp, cols in species_to_cols.items():
    print(f"[check] {sp}: {len(cols)} columns")

# --- robust column median (ignore NaN and zeros) ---
def col_median(series: pd.Series) -> float:
    vals = pd.to_numeric(series, errors="coerce")
    vals = vals.replace(0, np.nan)  # treat zero as missing
    return float(vals.median(skipna=True))

# --- compute scales and apply: scale = species_global / col_median ---
scales = {}          # (species, col) -> scale (full precision)
species_summary = [] # for logging

for sp, cols in species_to_cols.items():
    col_meds = {c: col_median(df[c]) for c in cols}
    valid = [v for v in col_meds.values() if np.isfinite(v) and v > 0]
    if not valid:
        for c in cols:
            scales[(sp, c)] = 1.0
        species_summary.append((sp, np.nan, {c: np.nan for c in cols}))
        continue

    # per-species global median (median of column medians)
    sp_global = float(np.median(valid))

    per_col_scales = {}
    for c in cols:
        m = col_meds[c]
        if np.isfinite(m) and m > 0:
            s = sp_global / m                 # full precision
        else:
            s = 1.0
        per_col_scales[c] = s
        scales[(sp, c)] = s

    species_summary.append((sp, sp_global, per_col_scales))

# --- apply scales in place (full precision) ---
for (sp, c), s in scales.items():
    if s != 1.0:
        df[c] = pd.to_numeric(df[c], errors="coerce") * s

# --- optional: show a brief summary in console (10 decimals for display only) ---
print("\n=== Per-species median normalization (first few species) ===")
for sp, sp_global, per_col_scales in species_summary[:8]:
    short = ", ".join(f"{k}:×{round10(per_col_scales[k]):.10f}" for k in list(per_col_scales.keys())[:4])
    gdisp = "NA" if not np.isfinite(sp_global) else f"{sp_global:.4g}"
    print(f"[{sp if sp else '(no-prefix)'}] global={gdisp}  |  {short}{' ...' if len(per_col_scales)>4 else ''}")

# --- save normalized table ---
Path(OUT_PATH).parent.mkdir(parents=True, exist_ok=True)
df.to_csv(OUT_PATH, sep="\t", index=False)
print(f"\n[OK] wrote normalized table → {OUT_PATH}")

