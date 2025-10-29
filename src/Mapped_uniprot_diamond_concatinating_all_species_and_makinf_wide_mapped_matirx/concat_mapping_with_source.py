# concatinate mapped ids
#!/usr/bin/env python3
"""
concat_mappings_with_source.py
    python concat_mappings_with_source.py  <dir_with_files>  <out.tsv>

Each input file must contain the columns:
    query_id    Entry
The script adds       SourceFile
"""

import sys, glob, pandas as pd
from pathlib import Path

if len(sys.argv) != 3:
    sys.exit("USAGE  python concat_mappings_with_source.py  <dir>  <out.tsv>")

IN_DIR   = Path(sys.argv[1]).resolve()
OUT_FILE = Path(sys.argv[2]).resolve()

# pattern: *final.(tsv|csv|xls|xlsx)
patterns = ["*final.tsv", "*final.csv", "*final.xls", "*final.xlsx"]
files = [p for pat in patterns for p in IN_DIR.glob(pat)]
if not files:
    sys.exit("No mapping files ending in '*final.*' found in " + str(IN_DIR))

print("Found", len(files), "files:")
for f in files: print("  •", f.name)

frames = []
for f in files:
    suffix = f.suffix.lower()
    if suffix in {".xls", ".xlsx"}:
        df = pd.read_excel(f, dtype=str)
    elif suffix == ".csv":
        df = pd.read_csv(f, sep=",", dtype=str)
    else:                            # .tsv (or anything else → treat as tab)
        df = pd.read_csv(f, sep="\t", dtype=str)

    # keep only the two columns (case-insensitive)
    df.columns = [c.lower() for c in df.columns]
    if not {"query_id", "entry"} <= set(df.columns):
        print("⚠  skipped (missing columns):", f.name)
        continue

    df = df[["query_id", "entry"]].rename(
        columns={"query_id": "query_id", "entry": "Entry"}
    )
    df["SourceFile"] = f.stem            # filename without extension
    frames.append(df)

if not frames:
    sys.exit("No usable mapping tables found.")

combined = pd.concat(frames, ignore_index=True)
combined.to_csv(OUT_FILE, sep="\t", index=False)

print(f"\n✓ wrote {len(combined):,} rows → {OUT_FILE}")


