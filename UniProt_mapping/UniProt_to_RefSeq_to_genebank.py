import os
import time
import json
import re
import requests
from pathlib import Path
from Bio import SeqIO
import pandas as pd

# -----------------------------
# CONFIG
# -----------------------------
FASTA_PATH = r'c:\Users\Sibel\Downloads\idmapping_2025_10_09.fasta'
OUTPUT_TSV = r'e:\Guido\sibel\Masters_Thesis\UniProt_id_matching\refseq_genbank_from_fasta.tsv'
OUTDIR = r'e:\Guido\sibel\Masters_Thesis\UniProt_id_matching\Test_uniprot_to_genebank_refseq'

PAUSE = 0.2  # polite spacing between API calls
RETRY = 3

UNIPROT_ACC_RE = re.compile(
    r"""
    (?<![A-Z0-9])                              # left boundary: not alnum
    (?:                                        # accession core
        [OPQ][0-9][A-Z0-9]{3}[0-9]             # classic 6-char starting O/P/Q
      | [A-NR-Z][0-9][A-Z0-9]{3}[0-9]          # classic 6-char others
      | [A-NR-Z][0-9][A-Z0-9]{8}               # 10-char (e.g., A0A0K2QJB4)
    )
    (?:-\d+)?                                  # optional isoform suffix: -2, -10, ...
    (?![A-Z0-9])                               # right boundary: not alnum
    """,
    re.VERBOSE
)


# --- Extract UniProt accessions from sp|..| / tr|..| or bare ---
def extract_uniprot_acc(header_id: str, header_desc: str) -> str | None:
    # Prefer UniProt pipe format if present: sp|ACC|..., tr|ACC|...
    if '|' in header_id:
        parts = header_id.split('|')
        if len(parts) >= 2 and UNIPROT_ACC_RE.fullmatch(parts[1]):
            return parts[1]
    # Fallback: scan id/description
    m = UNIPROT_ACC_RE.search(header_id) or UNIPROT_ACC_RE.search(header_desc)
    return m.group(0) if m else None

def extract_uniprot_acc(rec_id: str, rec_desc: str):
    # pipe format (most reliable)
    if "|" in rec_id:
        parts = rec_id.split("|")
        if len(parts) >= 2 and UNIPROT_ACC_RE.fullmatch(parts[1]):
            return parts[1]
    # scan description/id
    m = UNIPROT_ACC_RE.search(rec_desc) or UNIPROT_ACC_RE.search(rec_id)
    return m.group(0) if m else None

def iter_unique_accs(fasta_path):
    seen = set()
    for rec in SeqIO.parse(fasta_path, "fasta"):
        acc = extract_uniprot_acc(rec.id, rec.description)
        if acc and acc not in seen:
            seen.add(acc)
            yield acc

# --- UniProt helpers ---
def fetch_uniprot_json(acc: str):
    """
    Try exact accession; if 404 and it's an isoform (ACC-2), try base ACC.
    Return (json_or_None, note)
    """
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
    last_err = None
    for i in range(RETRY):
        try:
            r = requests.get(url, timeout=30)
            if r.status_code == 200:
                return r.json(), "ok"
            if r.status_code == 404 and "-" in acc:
                base = acc.split("-", 1)[0]
                r2 = requests.get(f"https://rest.uniprot.org/uniprotkb/{base}.json", timeout=30)
                if r2.status_code == 200:
                    return r2.json(), f"isoform_fallback:{base}"
                last_err = f"404 for {acc} and {base}"
            else:
                last_err = f"http_{r.status_code}"
        except Exception as e:
            last_err = f"error:{e}"
        time.sleep(PAUSE * (i + 1))
    return None, last_err or "failed"

def parse_xrefs(ujson: dict):
    refseq, genbank_like = set(), set()
    if not ujson:
        return refseq, genbank_like
    for x in (ujson.get("uniProtKBCrossReferences") or []):
        db = (x.get("database") or "").upper()
        xid = x.get("id") or ""
        if not xid:
            continue
        if db == "REFSEQ":
            refseq.add(xid)
        elif db in {"EMBL", "GENBANK", "DDBJ"}:
            genbank_like.add(xid)
    return refseq, genbank_like

def main():
    Path(OUTDIR).mkdir(parents=True, exist_ok=True)
    accs = list(iter_unique_accs(FASTA_PATH))
    print(f"Found {len(accs)} unique UniProt accessions from FASTA.")

    rows = []
    for i, acc in enumerate(accs, 1):
        print(f"[{i}/{len(accs)}] {acc}")
        query_acc = acc.split('-', 1)[0]  # use base accession if isoform like P12345-2
        ujson, note = fetch_uniprot_json(query_acc)
                                              
        time.sleep(PAUSE)
        refseq, genbank = parse_xrefs(ujson)
        rows.append({
            "uniprot": acc,
            "status": note,                  # ok / isoform_fallback:<base> / http_404 / http_5xx / error:...
            "refseq_count": len(refseq),
            "genbank_family_count": len(genbank),
            "refseq_ids": ",".join(sorted(refseq)),
            "genbank_family_ids": ",".join(sorted(genbank)),
        })

    audit = pd.DataFrame(rows)

    # Summary tables
    coverage = (
        audit.assign(has_any=audit["refseq_count"].gt(0) | audit["genbank_family_count"].gt(0))
             .groupby("status", as_index=False)["has_any"].agg(total="count", with_any="sum")
    )
    overall = {
        "total_acc": len(audit),
        "with_any": int((audit["refseq_count"] + audit["genbank_family_count"] > 0).sum()),
        "with_refseq": int(audit["refseq_count"].gt(0).sum()),
        "with_genbank_family": int(audit["genbank_family_count"].gt(0).sum()),
    }

    # Save TSVs
    p_audit = Path(OUTDIR) / "uniprot_xref_audit.tsv"
    p_missing = Path(OUTDIR) / "no_xrefs.tsv"
    p_coverage = Path(OUTDIR) / "coverage_by_status.tsv"
    p_overall = Path(OUTDIR) / "overall_counts.tsv"

    audit.to_csv(p_audit, sep="\t", index=False)
    audit[(audit["refseq_count"] == 0) & (audit["genbank_family_count"] == 0)].to_csv(p_missing, sep="\t", index=False)
    coverage.to_csv(p_coverage, sep="\t", index=False)
    pd.DataFrame([overall]).to_csv(p_overall, sep="\t", index=False)

    print("\n== Overall ==")
    for k, v in overall.items():
        print(f"{k}: {v}")
    print("\nSaved:")
    for p in [p_audit, p_missing, p_coverage, p_overall]:
        print(" -", p)

if __name__ == "__main__":
    main()
