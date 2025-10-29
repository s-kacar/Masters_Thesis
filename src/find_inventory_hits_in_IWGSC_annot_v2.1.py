#!/usr/bin/env python3
"""
find_inventory_hits_v2.py
Usage:
  python find_inventory_hits_in_IWGSC_annot_v2.1.py  inventory.xlsx  second_table.tsv  matched_rows.tsv  audit.tsv

What it does
- Inventory: read Excel 'inventory.xlsx' and build a set of UniProt accessions from
  'Entry (UniProt)' (fallback to 'Entry Name (UniProt)').
- Second table: read TSV with columns 'f.name' and 'f.type'.
- Consider only rows where f.type == 'Blast-Hit-Accession'.
- Parse accessions from f.name:
    * If pipe-delimited like 'sp|P12345|PROT_HUMAN' → take middle token 'P12345'.
    * If no pipes → use token as-is.
    * Multiple items per cell split by [comma|semicolon|whitespace].
- Mark rows with at least one accession present in the inventory.
- Save matching rows to matched_rows.tsv.
- Save audit (totals + per f.type counts) to audit.tsv.
"""

import sys
import time
import re
from pathlib import Path
import pandas as pd

if len(sys.argv) != 5:
    sys.exit("USAGE: python find_inventory_hits_v2.py inventory.xlsx second_table.tsv matched_rows.tsv audit.tsv")

INV_XLS  = Path(sys.argv[1]).resolve()
SECOND   = Path(sys.argv[2]).resolve()
OUT_HITS = Path(sys.argv[3]).resolve()
OUT_AUD  = Path(sys.argv[4]).resolve()
# add after your existing arg parsing
SHEET = sys.argv[5] if len(sys.argv) == 6 else "Sheet0"  # or pick your default



t0 = time.perf_counter()

# ───────────────────────── helpers ─────────────────────────
SPLIT_MULTI = re.compile(r"[,\s;]+")

def normalize_uniprot_token(tok: str) -> str:
    """
    Return a candidate UniProt accession from a token.
    Examples:
      'sp|P12345|PROT_HUMAN' → 'P12345'
      'tr|Q9ABC1|NAME_MOUSE' → 'Q9ABC1'
      'P02768'               → 'P02768'
    """
    if tok is None:
        return ""
    tok = tok.strip()
    if not tok:
        return ""
    if "|" in tok:
        parts = tok.split("|")
        # common UniProt format is db|ACC|NAME  → index 1
        if len(parts) >= 3 and parts[1].strip():
            return parts[1].strip()
        # fallback: if exactly 2 tokens, try the second
        if len(parts) == 2 and parts[1].strip():
            return parts[1].strip()
        # otherwise, take the middle-most non-empty token
        mids = [p.strip() for p in parts if p.strip()]
        if len(mids) >= 2:
            return mids[len(mids)//2]
        return mids[0] if mids else ""
    return tok

def explode_accessions(cell: str):
    """
    Split a cell into candidate tokens (by comma/semicolon/space),
    normalize each to a UniProt-style accession, and return a unique set.
    """
    if not isinstance(cell, str) or not cell.strip():
        return set()
    tokens = [t for t in SPLIT_MULTI.split(cell.strip()) if t]
    normed = {normalize_uniprot_token(t) for t in tokens}
    return {x for x in normed if x}

# ─────────────────── 1) load inventory ────────────────────
# replace your current read_excel line with:
inv_df = pd.read_excel(INV_XLS, sheet_name=SHEET, dtype=str)  # optionally add engine="openpyxl"

inv_col = None
for cand in ["Entry"]:
    if cand in inv_df.columns:
        inv_col = cand
        break
if inv_col is None:
    sys.exit("ERROR: Could not find an inventory column among "
             "['Entry (UniProt)', 'Entry Name (UniProt)', 'Entry', 'UniProt', 'Accession']")

inventory_accessions = (
    inv_df[inv_col]
      .dropna()
      .astype(str)
      .str.strip()
      .map(normalize_uniprot_token)
)
inventory_set = {x for x in inventory_accessions if x}
print(f"Inventory accessions: {len(inventory_set):,}")

# ─────────────────── 2) load second table ─────────────────
# Expecting TSV; change sep if yours is CSV.
df = pd.read_csv(SECOND, dtype=str)
n_total_rows = len(df)

# sanity for required columns
for required in ["f.name", "f.type"]:
    if required not in df.columns:
        sys.exit(f"ERROR: '{required}' column is missing in {SECOND.name}")

# ─────────── 3) restrict to 'Blast-Hit-Accession' rows ────
is_blast = df["f.type"].fillna("").eq("Blast-Hit-Accession")
df["__is_blast"] = is_blast
n_blast_rows = int(is_blast.sum())

# ─────────── 4) parse & compare accessions from f.name ────
# Build a set per row of normalized candidate accessions
df["__acc_set"] = df["f.name"].map(explode_accessions)

# Row has a hit if any accession intersects inventory_set (only check for Blast-Hit-Accession rows)
def has_match(row) -> bool:
    if not row["__is_blast"]:
        return False
    accs = row["__acc_set"]
    if not accs:
        return False
    return not inventory_set.isdisjoint(accs)

df["__has_match"] = df.apply(has_match, axis=1)

# matched rows (subset of 'Blast-Hit-Accession')
matched = df[df["__has_match"]].copy()

# ─────────── 5) write matched rows ────────────────────────
# Drop helper columns in the output table
to_drop = ["__is_blast", "__acc_set", "__has_match"]
matched_out = matched.drop(columns=[c for c in to_drop if c in matched.columns])
matched_out.to_csv(OUT_HITS, sep="\t", index=False)

# ─────────── 6) build audit table ─────────────────────────
audit_rows = []

# global summary
audit_rows.append({
    "metric": "total_rows_in_second_doc",
    "value": n_total_rows
})
audit_rows.append({
    "metric": "total_blast_hit_rows",
    "value": n_blast_rows
})
audit_rows.append({
    "metric": "matched_blast_hit_rows",
    "value": int(matched.shape[0])
})
audit_rows.append({
    "metric": "unmatched_blast_hit_rows",
    "value": int(n_blast_rows - matched.shape[0])
})
# unique accession coverage among BLAST rows
blast_acc_sets = df.loc[is_blast, "__acc_set"]
all_blast_accs = set().union(*blast_acc_sets) if len(blast_acc_sets) else set()
matched_accs = all_blast_accs.intersection(inventory_set)
audit_rows.append({
    "metric": "unique_blast_accessions_total",
    "value": len(all_blast_accs)
})
audit_rows.append({
    "metric": "unique_blast_accessions_matched_inventory",
    "value": len(matched_accs)
})

# per f.type breakdown (counts and matches)
type_counts = df["f.type"].fillna("NA").value_counts().to_dict()
for ftype, cnt in type_counts.items():
    sub = df[df["f.type"].fillna("NA").eq(ftype)]
    audit_rows.append({
        "metric": f"rows_in_type::{ftype}",
        "value": int(cnt)
    })
    audit_rows.append({
        "metric": f"matched_rows_in_type::{ftype}",
        "value": int(sub["__has_match"].sum())
    })

audit_df = pd.DataFrame(audit_rows, columns=["metric", "value"])
audit_df.to_csv(OUT_AUD, sep="\t", index=False)

# ─────────── 7) console summary ───────────────────────────
print(f"Second document rows (total):         {n_total_rows:,}")
print(f"Rows with f.type == Blast-Hit-Accession: {n_blast_rows:,}")
print(f"Matched Blast-Hit rows:                {matched.shape[0]:,}")
print(f"✓ wrote matches → {OUT_HITS}")
print(f"✓ wrote audit   → {OUT_AUD}")
print(f"Done in {time.perf_counter()-t0:.1f}s")
