#!/usr/bin/env python3
"""
find_inventory_hits.py
    python find_inventory_hits.py  inventory.xlsx  wide.tsv  hits.tsv
• Reads sheet "Wheat Inventory" with column 'Entry (UniProt)'
• Scans every Entry_<species> column in wide.tsv
• Keeps rows where *any* Entry cell contains at least one UniProt ID
  from the inventory list.
"""

import sys, time, pandas as pd
from pathlib import Path

if len(sys.argv) not in (4, 5):
    sys.exit("USAGE  python find_inventory_hits.py  inventory.xlsx  Finale.tsv  hits.tsv  [Entry_column]")

INV_XLS  = Path(sys.argv[1]).resolve()
WIDE_TSV = Path(sys.argv[2]).resolve()
OUT_TSV  = Path(sys.argv[3]).resolve()
ENTRY_COL_ARG = sys.argv[4] if len(sys.argv) == 5 else None   # optional


t0 = time.perf_counter()

# ── 1. load inventory list ─────────────────────────────────────────────
inv_df = pd.read_excel(INV_XLS, sheet_name="Triticum aestivum", dtype=str)
inventory = (
    inv_df['Entry (UniProt)']
    .dropna()
    .str.strip()
    .unique()
    .tolist()
)
inventory_set = set(inventory)
print(f"Inventory IDs: {len(inventory_set):,}")

# ── 2. load wide table ─────────────────────────────────────────────────
wide = pd.read_csv(WIDE_TSV, sep="\t")

# find every Entry_<species> column automatically
entry_cols = [c for c in wide.columns if c.startswith("Entry_")]
print(f"Entry columns detected: {len(entry_cols)}")

# ── 3. mark rows that contain ≥1 inventory ID in ANY Entry column ——––
# helper
def cell_has_hit(cell) -> bool:
    if not isinstance(cell, str) or not cell:
        return False
    return bool(inventory_set.intersection(id.strip() for id in cell.split(",")))

# row mask
hit_mask = (
    wide[entry_cols]
        .apply(lambda col: col.map(cell_has_hit))
        .any(axis=1)
)

# ── build row mask (True = row contains ≥1 inventory UniProt ID) ──────
hit_mask = (
    wide[entry_cols]
        .apply(lambda col: col.map(cell_has_hit))   # column-wise map
        .any(axis=1)
)

# NEW — slice the wide table with that mask
hits = wide[hit_mask]

print(f"Rows with at least one inventory protein: {hits.shape[0]:,}")

# ── save ───────────────────────────────────────────────────────────────
hits.to_csv(OUT_TSV, sep="\t", index=False)
print(f"✓ wrote {OUT_TSV}  in {time.perf_counter()-t0:.1f}s")



