#!/usr/bin/env python3
"""
map_uniprot_into_wide.py
    python map_uniprot_into_wide.py  hog_wide.tsv  all_mappings.tsv  hog_wide_with_uni.tsv
"""

import sys, time, re, pandas as pd
from pathlib import Path

if len(sys.argv) != 4:
    sys.exit("USAGE  python map_uniprot_into_wide.py  wide.tsv  mappings.tsv  out.tsv")

WIDE   = Path(sys.argv[1])
MAP_T  = Path(sys.argv[2])
OUT_T  = Path(sys.argv[3])

t0 = time.perf_counter()

# ── 1. load wide + mapping tables ─────────────────────────────────────
wide = pd.read_csv(WIDE, sep="\t", dtype=str).fillna("")
maps = pd.read_csv(MAP_T, sep="\t", dtype=str)[["query_id", "Entry"]].dropna()

id2entry = dict(zip(maps["query_id"].str.strip(), maps["Entry"].str.strip()))
print(f"Loaded wide  {wide.shape[0]:,} rows × {wide.shape[1]} cols")
print(f"Mappings     {len(id2entry):,} ID pairs   ({time.perf_counter()-t0:.1f}s)")

# ── 2. detect species-only columns (exact species name, no extra '_') ──
import re

# ── detect species-only ID columns ─────────────────────────────────────
species_cols = [
    c for c in wide.columns
    if re.fullmatch(r"[A-Za-z]+_[A-Za-z]+", c)
]
print("Species ID columns found:", len(species_cols))
# e.g. ['Triticum_aestivum', 'Hordeum_vulgare', 'Zea_mays', …]
print("\nDEBUG – first 20 column names:")
print(list(wide.columns[:20]))

print("\nDEBUG – species_cols detected:")
print(species_cols)
# ── 3. helper to map a single cell ────────────────────────────────────
def map_cell(cell: str) -> str:
    if not cell:
        return ""
    uni = [id2entry.get(pid.strip()) for pid in cell.split(",")]
    uni = [u for u in uni if u]      # drop None
    return ",".join(uni)

# ── 4. add / overwrite Entry_<species> columns ────────────────────────
for sp_col in species_cols:
    entry_col = f"Entry_{sp_col}"
    mapped = wide[sp_col].apply(map_cell)

    if entry_col in wide.columns:          # already exists → overwrite
        wide[entry_col] = mapped
    else:                                  # insert next to ProteinIDs column
        wide.insert(
            wide.columns.get_loc(sp_col) + 1,
            entry_col,
            mapped
        )

# ── 5. save ───────────────────────────────────────────────────────────
wide.to_csv(OUT_T, sep="\t", index=False)
print(f"✓ wrote {OUT_T} in {time.perf_counter()-t0:.1f}s")
