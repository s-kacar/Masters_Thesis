#!/usr/bin/env python3
# map_protein_names_minimal.py

import re, sys
import pandas as pd
from pathlib import Path

# ── EDIT THESE THREE PATHS ────────────────────────────────────────────
HOG_PATH      = r"E:\Guido\sibel\Masters_Thesis\data\protein_names_mapped_to_Orthogroups_intensity_MEAN.MEDIAN_NORM.tsv"
MAPPING_PATH  = r"E:\Guido\sibel\import\wheat import (this one for saving the files, tha path pipeline runs from is (impot. then your file name)\storage_protein_inventory_template_OG_24_07_deconvoluted.xlsx"
OUT_PATH      = r"E:\Guido\sibel\Masters_Thesis\data\protein_names_mapped_to_Orthogroups_intensity_mean.tsv"

# If your Excel has multiple sheets, set the correct one here:
SHEET_NAME = "Triticum aestivum"   # or an int like 0

# Column names in the mapping file (exact, after normalization)
ID_COL   = "Entry (UniProt)"
NAME_COL = "Protein Name"

# HOG column prefix to map & output prefix
ENTRY_PREFIX = "Entry_"
OUT_PREFIX   = "ProteinName_"

# ── HELPERS ───────────────────────────────────────────────────────────
def read_any_table(path, sheet=None, dtype=None):
    p = Path(path); suf = p.suffix.lower()
    if suf in {".xlsx", ".xls"}:
        # force a single DataFrame (first sheet if None)
        sheet_name = 0 if sheet is None else sheet
        return pd.read_excel(p, sheet_name=sheet_name, dtype=dtype, engine="openpyxl")
    if suf == ".tsv":
        return pd.read_csv(p, sep="\t", dtype=dtype)
    if suf == ".csv":
        return pd.read_csv(p, dtype=dtype)
    # fallback
    return pd.read_csv(p, sep="\t", dtype=dtype)

def build_id_to_name(mapping_df, id_col=ID_COL, name_col=NAME_COL):
    id_to_name = {}
    for _, row in mapping_df[[id_col, name_col]].dropna(subset=[id_col]).iterrows():
        raw  = str(row[id_col])
        name = str(row.get(name_col, "")).strip()
        ids = [x.strip() for x in re.split(r"[,\s;]+", raw) if x.strip()]
        for acc in ids:
            if acc not in id_to_name or not id_to_name[acc]:
                id_to_name[acc] = name
    return id_to_name

def ordered_unique(seq):
    seen, out = set(), []
    for x in seq:
        if x and x not in seen:
            seen.add(x); out.append(x)
    return out

def map_collapsed_ids_to_names(df, id_to_name,
                               entry_prefix=ENTRY_PREFIX,
                               out_prefix=OUT_PREFIX,
                               split_regex=r"[,\s]+",
                               joiner="; "):
    entry_cols = [c for c in df.columns if c.startswith(entry_prefix)]
    if not entry_cols:
        print(f"[ERROR] No columns starting with '{entry_prefix}' found.", file=sys.stderr)
        return df
    for ec in entry_cols:
        parts = (
            df[ec].fillna("").astype(str)
              .apply(lambda s: [t.strip() for t in re.split(split_regex, s) if t.strip()])
        )
        mapped_lists = []
        for id_list in parts:
            names = [id_to_name.get(uid, "") for uid in id_list]
            names = [n for n in names if n]
            names = ordered_unique(names)
            mapped_lists.append(names)
        df[out_prefix + ec.replace(entry_prefix, "", 1)] = [
            joiner.join(names) if names else "" for names in mapped_lists
        ]
    return df

# ── SCRIPT BODY ───────────────────────────────────────────────────────
# 0) sanity checks
if not Path(HOG_PATH).exists():
    sys.exit(f"[ERROR] HOG file not found: {HOG_PATH}")
if not Path(MAPPING_PATH).exists():
    sys.exit(f"[ERROR] Mapping file not found: {MAPPING_PATH}")

# 1) read inputs
hog = read_any_table(HOG_PATH)
mapping = read_any_table(MAPPING_PATH, sheet=SHEET_NAME)

# normalize mapping headers (kill weird spaces, unify)
mapping.columns = (
    mapping.columns.astype(str)
    .str.replace("\u00A0", " ", regex=False)
    .str.replace(r"\s+", " ", regex=True)
    .str.strip()
)

print("[INFO] HOG columns (first 15):", list(hog.columns[:15]))
print("[INFO] Mapping columns (first 15):", list(mapping.columns[:15]))

# 2) build dictionary
if ID_COL not in mapping.columns or NAME_COL not in mapping.columns:
    sys.exit(f"[ERROR] Expected columns not found in mapping.\n"
             f"  Need: '{ID_COL}' and '{NAME_COL}'\n"
             f"  Got:  {list(mapping.columns)}")

id_to_name = build_id_to_name(mapping, id_col=ID_COL, name_col=NAME_COL)
if not id_to_name:
    sys.exit("[ERROR] Mapping dictionary is empty. Check mapping headers/content.")

# 3) ensure HOG has Entry_* columns
entry_cols = [c for c in hog.columns if c.startswith(ENTRY_PREFIX)]
if not entry_cols:
    sys.exit(f"[ERROR] No columns starting with '{ENTRY_PREFIX}' in HOG.\n"
             f"Available: {list(hog.columns)}")

# 4) map
hog_out = map_collapsed_ids_to_names(
    hog.copy(), id_to_name,
    entry_prefix=ENTRY_PREFIX,
    out_prefix=OUT_PREFIX,
    split_regex=r"[,\s]+",
    joiner="; "
)

# 5) save
Path(OUT_PATH).parent.mkdir(parents=True, exist_ok=True)
hog_out.to_csv(OUT_PATH, sep="\t", index=False)
print(f"[OK] wrote {OUT_PATH}")
