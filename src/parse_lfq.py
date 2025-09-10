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
id_col        = "Protein ID"          # column that has search IDs
lfq_prefix    = " MaxLFQ Intensity"       # every LFQ column starts with this
tissue_lookup = lambda col: col.split(lfq_prefix)[1]  # pull tissue name
# ------------------------------------------------------------------
code_regex = re.compile(r"P\d+_(\d+)")

def tissue_lookup(col: str) -> str:
    """Extract tissue name from LFQ column header."""
    m = code_regex.search(colname)
    if not m:
        return colname
    code = m.group(1)

    code_map = {'93501_1 MaxLFQ Intensity': 'Primary leaf', '93502_2 MaxLFQ Intensity': 'Secondary leaf', '93503_3 MaxLFQ Intensity': 'Young root',
            '93504_4 MaxLFQ Intensity': 'Node', '93505_5 MaxLFQ Intensity': 'Internode', '93506_6 MaxLFQ Intensity': 'Adult root *',
            '93507_7 MaxLFQ Intensity': 'Adult root **', '93508_8 MaxLFQ Intensity': 'Adult root ***', '93509_9 MaxLFQ Intensity': 'Anther *',
            '93510_10 MaxLFQ Intensity': 'Anther **', '93511_11 MaxLFQ Intensity': 'Pollen', '93512_12 MaxLFQ Intensity': 'Stigma, style and ovary *',
            '93513_13 MaxLFQ Intensity': 'Stigma, style and ovary **', '93514_14 MaxLFQ Intensity': 'Immature seed *', '93515_15 MaxLFQ Intensity': 'Immature seed **',
            '93516_16 MaxLFQ Intensity': 'Mature seed'}

    return code_map.get(code, f"tissue_{code}")
python src/src/parse_lfq.py --in "data/database_search_results/ENB/wheat/combined_protein.tsv" --out "tmp/lfq_wheat.tsv" --species wheat












df = pd.read_csv(IN, sep="\t")

# keep ID col + all LFQ cols
lfq_cols = [c for c in df.columns if c.endswith(lfq_prefix)]
slim     = df[[id_col] + lfq_cols]

# reshape
tidy = (
    slim
    .melt(id_vars=id_col,
          var_name="TissueCol",
          value_name="LFQ")
    .assign(
        Tissue=lambda d: d.TissueCol.apply(tissue_lookup),
        Species=args.species,
        ProteinID=lambda d: d[id_col]
    )
    .drop(columns=["TissueCol", id_col])
)

tidy.to_csv(OUT, sep="\t", index=False)
print(f"✅ wrote {len(tidy):,} rows → {OUT.relative_to(Path.cwd())}")