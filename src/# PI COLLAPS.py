#!/usr/bin/env python3
# collapse_wide_no_sample.py
#
# Input  : TSV or CSV with columns
#            HOG, Species, ProteinID, <many MaxLFQ columns>
# Output : Wide CSV (one row per HOG):
#            HOG, ProteinIDs_<sp>, <sp>_<origCol>, …  (for every species)
#            Empty cells are "" (no NaN)

import sys
from pathlib import Path
import pandas as pd


def tidy(text: str) -> str:
    """Make text safe for column names: spaces → '_', commas removed."""
    return str(text).replace(" ", "_").replace(",", "")


def collapse(input_path: Path, output_path: Path) -> None:
    # 1. read file (auto-detect tab or comma)
    df = pd.read_csv(input_path, sep="\t")
    print("Loaded", df.shape[0], "rows with", df.shape[1], "columns")

    # 2. basic checks
    for col in ("HOG", "Species", "ProteinID"):
        if col not in df.columns:
            sys.exit(f"Missing column '{col}' in {input_path}")

    # 3. tidy species names once
    df["Species"] = df["Species"].apply(tidy)

    # 4. identify all MaxLFQ intensity columns
    intensity_cols = [c for c in df.columns if "MaxLFQ" in c]
    if not intensity_cols:
        sys.exit("No columns containing 'MaxLFQ' were found.")

    # 5. build the wide table row-by-row
    rows = []

    PROGRESS_STEP = 500

    for idx, (hog, hog_block) in enumerate(df.groupby("HOG"), 1):
        if idx % PROGRESS_STEP == 0:
            print(f"Processed {idx} HOGs…", flush=True)
        if idx >= 23300:
            print("Currently on", hog, flush=True)

        row = {"HOG": hog}


        for species, sp_block in hog_block.groupby("Species"):
            # ----------------------------------------------
            # 1) protein list for this HOG × species
            ids = sp_block["ProteinID"].dropna().astype(str).unique()
            id_list = ",".join(ids)
            row[f"ProteinIDs_{species}"] = id_list

            # *** NEW RULE: if no proteins for this species,
            #               skip all its LFQ columns ***
            if not id_list:           # id_list == ""  ➜ leave LFQ cells blank
                continue
            # ----------------------------------------------

            # 2) write every LFQ column for this species
            for col in intensity_cols:
                short = col.split()[0]                 # e.g. "P095141_1"
                header = f"{species}_{short}"
                vals = sp_block[col].dropna().astype(str)
                row[header] = ",".join(vals)

        rows.append(row)

    # 6. rows → DataFrame
    print("Loop finished.  Building DataFrame…", flush=True)
    wide = pd.DataFrame(rows)

    # 7. order columns: ProteinIDs_<sp> then that species’ intensity columns
    species_order = sorted(
        {c.split("_", 1)[1] for c in wide.columns if c.startswith("ProteinIDs_")}
    )
    ordered = ["HOG"]
    for sp in species_order:
        ordered.append(f"ProteinIDs_{sp}")
        ordered.extend(sorted([c for c in wide.columns if c.startswith(f"{sp}_")]))

    print("Reordering columns…", flush=True)
    wide = wide.reindex(columns=ordered)

    # 8. write CSV
    print("Writing TSV to disk…", flush=True)
    wide.to_csv(output_path, sep="\t", index=False, na_rep="")

    print(
        f"✓ Done!  Wrote '{output_path}'  "
        f"[{wide.shape[0]} HOG rows x {wide.shape[1]} columns]"
    )

# ── run from command line ───────────────────────────────────────────────
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(
            "\nUSAGE  python collapse_wide_no_sample.py  input.tsv  output.csv\n"
        )
        sys.exit(0)

    collapse(Path(sys.argv[1]), Path(sys.argv[2]))

