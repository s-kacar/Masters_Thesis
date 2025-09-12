#!/usr/bin/env python3
"""
Turn OrthoFinder's Orthogroups.tsv into a long table:

Orthogroup | Species | ProteinID
OG0000001  | wheat   | TraesCS1A02G012300
OG0000001  | wheat   | TraesCS1A02G012400
OG0000001  | barley  | BAJ97982
…
Explode: for the cell containing many protein IDs(comma seperated) for that species will also be melted with the species columns and its HOG.
Usage:
  python explode_orthogroups.py Orthogroups.tsv og_long.tsv
"""

import sys, re
import pandas as pd
from pathlib import Path
import numpy as np

ortho_tsv, out_tsv = map(Path, sys.argv[1:3])

# read as raw table (tab-separated, first column is "Orthogroup")
df = pd.read_csv(ortho_tsv, sep="\t")

meta_cols = ["HOG"]
drop_cols = ["OG", "Gene Tree Parent Clade"]
species_cols = [c for c in df.columns if c not in meta_cols + drop_cols]
# reshape to long form: one column for species, one for "IDs"
long = (df.melt(id_vars=meta_cols, value_vars=species_cols, var_name="Species", value_name="IDs").assign(IDs=lambda d: d.IDs.fillna("")).assign(ProteinID=lambda d: d["IDs"].str.split(r", \s*")).explode("ProteinID").drop(columns=["IDs"]).query("ProteinID != ''").reset_index(drop=True))

long.to_csv(out_tsv, sep="\t", index=False)
print(f"Wrote {len(long):,} rows → {out_tsv}")
