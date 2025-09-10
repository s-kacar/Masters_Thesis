#!/usr/bin/env python3
"""
parse_lfq.py  –  turn FragPipe protein table into tidy long-format
Usage:
  python src/parse_lfq.py \
         --in  data/search_results/proteinGroups_wheat.txt \
         --out tmp/lfq_wheat.tsv \
         --species wheat
"""

import argparse
import pandas as pd
from pathlib import Path
import re

p = argparse.ArgumentParser()
p.add_argument("--in",  required=True, dest="infile",
               help="raw proteinGroups or FragPipe TSV")
p.add_argument("--out", default="tmp/lfq.tsv",
               help="output TSV inside the repo")
p.add_argument("--species", required=True,
               help="short species name (matches OrthoFinder column)")
args = p.parse_args()

IN   = Path(args.infile).expanduser().resolve()
OUT  = Path(args.out).expanduser().resolve()
OUT.parent.mkdir(parents=True, exist_ok=True)

# --------- tweak these three strings if your file differs ---------
# ---------- tweak these strings & function if your file differs ----------
id_col     = "Protein ID"                # column that has protein IDs
lfq_suffix = " MaxLFQ Intensity"         # common tail of every LFQ column

code_regex = re.compile(r"P\d+_(\d+)")   # grabs the digits after the underscore

def tissue_lookup(colname: str) -> str:
    """
    Turn a full LFQ header like
    'P093501_1 MaxLFQ Intensity'  →  'Primary leaf'
    """
    m = code_regex.search(colname)
    if not m:
        return colname                 # fallback: give raw header

    code = m.group(1)                  # "1", "2", …

    code_map = {                       # fill in or edit as needed
        "1":  "Primary leaf",
        "2":  "Secondary leaf",
        "3":  "Young root",
        "4":  "Node",
        "5":  "Internode",
        "6":  "Adult root *",
        "7":  "Adult root **",
        "8":  "Adult root ***",
        "9":  "Anther *",
        "10": "Anther **",
        "11": "Pollen",
        "12": "Stigma/style/ovary *",
        "13": "Stigma/style/ovary **",
        "14": "Immature seed *",
        "15": "Immature seed **",
        "16": "Mature seed",
    }

    return code_map.get(code, f"tissue_{code}")

df = pd.read_csv(IN, sep="\t")

# keep ID col + all LFQ cols
lfq_cols = [c for c in df.columns if c.endswith(lfq_suffix)]
slim     = df[[id_col] + lfq_cols]

# reshape
tidy = (
    slim
    .melt(id_vars=id_col,
          var_name="TissueCol",
          value_name="MaxLFQ")
    .assign(
        Tissue=lambda d: d.TissueCol.apply(tissue_lookup),
        Species=args.species,
        ProteinID=lambda d: d[id_col]
    )
    .drop(columns=["TissueCol", id_col])
)

tidy["MaxLFQ"] = pd.to_numeric(tidy["MaxLFQ"], errors="coerce")

tidy.to_csv(OUT, sep="\t", index=False)
print(f"✅ wrote {len(tidy):,} rows → {OUT.relative_to(Path.cwd())}")