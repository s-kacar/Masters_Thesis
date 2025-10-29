#!/usr/bin/env python3
# mean_per_cell_ignore_multival.py
# Usage: python mean_per_cell_ignore_multival.py  input.tsv  output.tsv

import sys, re, math
import pandas as pd
from pathlib import Path

if len(sys.argv) != 3:
    print("USAGE  python mean_per_cell_ignore_multival.py  input.tsv  output.tsv")
    sys.exit(0)

IN  = Path(sys.argv[1])
OUT = Path(sys.argv[2])

df = pd.read_csv(IN, sep="\t", dtype=str).fillna("")

# 1) pick intensity columns only
intensity_cols = [c for c in df.columns if re.search(r"_P\d", c) or "MaxLFQ" in c]

# 2) helper function
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
