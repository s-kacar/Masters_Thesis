import argparse, glob, pandas as pd
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('--root', default='.', help="project root folder")
p.add_argument("--hog", default="tmp/hog_log.tsv", help="Hog long table relative to --root")
p.add_argument("--lfqdir", default="tmp", help="folder with lfq_*.tsv files (relative to --root)")
p.add_argument("--out", default="tmp/hog_lfq_long.tsv", help="output tsv also relative to --root")
args = p.parse_args()

ROOT = Path(args.root).expanduser().resolve()
hog_tsv = (ROOT / args.hog).resolve()
lfq_dir = (ROOT / args.lfqdir).resolve()
out_tsv = (ROOT / args.out).resolve()
out_tsv.parent.mkdir(parents=True, exist_ok=True)
print(f"Reading HOG data from {hog_tsv}")
# Read HOG file
hog = pd.read_csv(hog_tsv, sep="\t")

# Read all LFQ files
lfq_files = glob.glob(str(lfq_dir / "lfq_*.tsv"))
print(f"Reading LFQ data from {len(lfq_files)} files in {lfq_dir}")

if not lfq_files:
    raise SystemExit(f"No LFQ files found in {lfq_dir}")

# Columns expected in LFQ files
lfq_cols = ['Tissue', 'Species', 'ProteinID', 'MaxLFQ']

# Read and concatenate all LFQ files
lfq_list = []
for f in lfq_files:
    df = pd.read_csv(f, sep="\t", usecols=lfq_cols)
    lfq_list.append(df)

print(df.head())
df.to_csv("testing append on merge.tsv", sep="\t", index=False)

lfq = pd.concat(lfq_list, ignore_index=True)

# Optional: merge LFQ with HOG on ProteinID
merged = pd.merge(lfq, hog, on='ProteinID', how='left')  # or how='inner' if you want only matches

print(merged.head())
merged.to_csv(out_tsv, sep="\t", index=False)

print(f"Wrote {len(merged):,} rows data to {out_tsv.relative_to(ROOT)}")