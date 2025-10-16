# -*- coding: utf-8 -*-
"""
TSV duplicate finder + mapping report + overlap with another TSV

What you get (in OUTDIR):
  1) duplicate_summary.tsv
     - main_id, dup_count, mapped_ids (unique, ';'-joined)
  2) duplicate_pairs.tsv
     - one row per duplicate pair: main_id, mapped_id (exploded), row_index
  3) unique_ids_file1.tsv
     - unique main_id list from file1
  4) overlap_report.tsv
     - id, in_file1, in_file2, where (overlap/only_file1/only_file2)
"""

import pandas as pd
import numpy as np
from pathlib import Path


# ---- CONFIG ----# now Excel
              # <-- set your sheet name
FILE1 = r"E:\Guido\sibel\Masters_Thesis\Uniprot_id_matching\Wheat ReseechacekX - Copy.tsv"   # the file where you want to find duplicates
FILE2 = r"c:\Users\Sibel\Downloads\idmapping_2025_10_08.xlsx"   # the file to compare overlap against
FILE2_SHEET = "Sheet0"
OUTDIR = r"E:\Guido\sibel\Masters_Thesis\Uniprot_id_matching\Inventory macth on RefSeq starting with XP_script.tsv"

# column names in file1:
COL_MAIN = "Entry"        # column to check for duplicates
COL_MAP  = "From"          # column on the same row you want to associate (mapping partner)

# column name in file2 to compare against:
COL_OTHER = "Entry"  

READ_KW_TSV = dict(sep="\t", dtype=str, keep_default_na=False, na_values=["", "NA", "NaN"])
READ_KW_XLSX = dict(dtype=str)  # read everything as str from Excel

def normalize_series(s: pd.Series) -> pd.Series:
    s = s.replace("", np.nan)
    return s.astype("string").str.strip()

def read_file1(path):
    return pd.read_csv(path, **READ_KW_TSV)

def read_file2(path, sheet):
    suffix = Path(path).suffix.lower()
    if suffix in [".xlsx", ".xls"]:
        return pd.read_excel(path, sheet_name=sheet, **READ_KW_XLSX)
    elif suffix in [".tsv", ".txt"]:
        return pd.read_csv(path, **READ_KW_TSV)
    else:
        raise ValueError(f"Unsupported FILE2 extension: {suffix}")

def main():
    Path(OUTDIR).mkdir(parents=True, exist_ok=True)

    # ---- load file1 (TSV)
    df1 = read_file1(FILE1)
    if COL_MAIN not in df1.columns or COL_MAP not in df1.columns:
        raise KeyError(f"file1 needs '{COL_MAIN}' and '{COL_MAP}'. Got: {df1.columns.tolist()}")
    df1[COL_MAIN] = normalize_series(df1[COL_MAIN])
    df1[COL_MAP]  = normalize_series(df1[COL_MAP])
    df1 = df1.dropna(subset=[COL_MAIN]).copy()
    df1["_row_index"] = range(len(df1))

    # ---- duplicate summary on file1
    counts = (df1.groupby(COL_MAIN, as_index=False).size()
                .rename(columns={"size":"dup_count"}))
    dup_only = counts[counts["dup_count"] > 1]

    mapped = (df1[df1[COL_MAIN].isin(dup_only[COL_MAIN])]
                .groupby(COL_MAIN, as_index=False)
                .agg(mapped_ids=(COL_MAP, lambda x: ";".join(pd.Series(x).dropna().unique()))))

    duplicate_summary = (dup_only.merge(mapped, on=COL_MAIN, how="left")
                               .sort_values(["dup_count", COL_MAIN], ascending=[False, True]))

    duplicate_pairs = (df1[df1[COL_MAIN].isin(dup_only[COL_MAIN])]
                        .rename(columns={COL_MAP:"mapped_id"})
                        [[COL_MAIN, "mapped_id", "_row_index"]]
                        .sort_values([COL_MAIN, "_row_index"]))

    unique_file1 = (df1[[COL_MAIN]].drop_duplicates()
                               .rename(columns={COL_MAIN:"id"})
                               .sort_values("id").reset_index(drop=True))

        # ---- load file2 (Excel or TSV)
    df2 = read_file2(FILE2, FILE2_SHEET)
    if COL_OTHER not in df2.columns:
        raise KeyError(f"file2 needs '{COL_OTHER}'. Got: {df2.columns.tolist()}")
    unique_file2 = (pd.DataFrame({"id": normalize_series(df2[COL_OTHER])})
                      .dropna().drop_duplicates()
                      .sort_values("id").reset_index(drop=True))

    # ---- overlap
    set1, set2 = set(unique_file1["id"]), set(unique_file2["id"])
    overlap = sorted(set1 & set2); only1 = sorted(set1 - set2); only2 = sorted(set2 - set1)
    overlap_df = pd.DataFrame({
        "id": overlap + only1 + only2,
        "in_file1": [True]*len(overlap) + [True]*len(only1) + [False]*len(only2),
        "in_file2": [True]*len(overlap) + [False]*len(only1) + [True]*len(only2),
    })
    overlap_df["where"] = overlap_df.apply(
        lambda r: "overlap" if (r.in_file1 and r.in_file2) else ("only_file1" if r.in_file1 else "only_file2"),
        axis=1
    )

    # ---- save TSV reports
    p_summary = Path(OUTDIR) / "duplicate_summary.tsv"
    p_pairs   = Path(OUTDIR) / "duplicate_pairs.tsv"
    p_uniq1   = Path(OUTDIR) / "unique_ids_file1.tsv"
    p_uniq2   = Path(OUTDIR) / "unique_ids_file2.tsv"
    p_overlap = Path(OUTDIR) / "overlap_report.tsv"

    duplicate_summary.to_csv(p_summary, sep="\t", index=False)
    duplicate_pairs.to_csv(p_pairs, sep="\t", index=False)
    unique_file1.to_csv(p_uniq1, sep="\t", index=False)
    unique_file2.to_csv(p_uniq2, sep="\t", index=False)
    overlap_df.to_csv(p_overlap, sep="\t", index=False)

    print("✅ Done")
    print("Saved:")
    for p in [p_summary, p_pairs, p_uniq1, p_uniq2, p_overlap]:
        print(" -", p)

if __name__ == "__main__":
    main()