import glob, re
import pandas as pd
from pathlib import Path

ROOT = Path(".").resolve()
HG_TSV = ROOT / "tmp/hog_long.tsv"
OUt_TSV = ROOT / "tmp/hog_intensity_long.tsv"
ID_COL = "Protein ID"
LFQ_TAIL = " MaxLFQ Intensity"

hog = pd.read_csv(HG_TSV, sep="\t")

folder2hog = {'wheat': 'GCF_018294505.1_Triticum_aestivum',
              'rice': 'GCF_034140825.1_Oryza_sativa', 'maize' : 'GCF_902167145.1_Zea_mays', 'barley': 'GCF_904849725.1_Hordeum_vulgare', 
              'sugarcane': 'Saccharum_officinarum_LA-Purple.protein', 'sorghum': 'GCF_000003195.3_Sorghum_bicolor', 'pearlmillet': 'GCA_963924085.1_Cenchrus_americanus.helixer'}

out_rows = []

for tsv in glob.glob(str(ROOT / "data/database_search_results/ENB/*/combined_protein.tsv")):
    folder = Path(tsv).parent.name
    hog_tag = folder2hog.get(folder)
    if hog_tag is None:
        print(f"!!!! skip {folder} (not mapped)")
        continue

    raw = pd.read_csv(tsv, sep="\t")
    lfq_cols = [c for c in raw.columns if c.endswith(LFQ_TAIL)]
    if not lfq_cols:
        print(f"!!!! no MaxLFQ found in {tsv}")     
        continue
    slim = raw[[ID_COL] + lfq_cols].copy()
    slim = slim.rename(columns={ID_COL: "ProteinID"})
    slim["Species"] = hog_tag

    merged = (hog.merge(slim, on=["ProteinID", "Species"], how="left")).dropna(how="all", subset=lfq_cols)
    out_rows.append(merged)

    final = pd.concat(out_rows, ignore_index=True)
    OUt_TSV.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(OUt_TSV, sep="\t", index=False)
    print(f"Wrote {len(final):,} rows to {OUt_TSV.relative_to(ROOT)}")


              