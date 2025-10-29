#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ---------- imports (keep tight & early) ----------
import sys, io, os, tempfile, atexit, traceback
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from decimal import Decimal, ROUND_HALF_UP, getcontext

# ---------- robust tee: file + console, line-buffered ----------
LOG_DIR  = Path(r"E:\Guido\sibel\Masters_Thesis\data")
LOG_DIR.mkdir(parents=True, exist_ok=True)
LOG_PATH = LOG_DIR / f"norm_log_{datetime.now():%Y%m%d_%H%M%S}.txt"

try:
    _log_f = open(LOG_PATH, "w", encoding="utf-8", buffering=1)  # line-buffered
except PermissionError:
    # fallback to temp if the target is locked or not writable
    LOG_PATH = Path(tempfile.gettempdir()) / f"norm_log_{datetime.now():%Y%m%d_%H%M%S}.txt"
    _log_f = open(LOG_PATH, "w", encoding="utf-8", buffering=1)

class Tee(io.TextIOBase):
    def __init__(self, *streams): self.streams = streams
    def write(self, s):
        for st in self.streams:
            try:
                st.write(s)
                st.flush()
            except Exception:
                pass
        return len(s)
    def flush(self):
        for st in self.streams:
            try: st.flush()
            except Exception: pass

_sys_stdout_orig = sys.stdout
_sys_stderr_orig = sys.stderr
sys.stdout = Tee(sys.stdout, _log_f)
sys.stderr = Tee(sys.stderr, _log_f)

@atexit.register
def _close_log():
    # flush + close safely; restore originals
    try:
        sys.stdout.flush(); sys.stderr.flush()
    except Exception:
        pass
    try:
        _log_f.flush(); _log_f.close()
    except Exception:
        pass
    sys.stdout = _sys_stdout_orig
    sys.stderr = _sys_stderr_orig

print(f"[log] writing console + file → {LOG_PATH}", flush=True)

# ---------- config / helpers ----------
getcontext().prec = 28

def round10(x: float) -> float:
    return float(Decimal(str(x)).quantize(Decimal("0.0000000000"), rounding=ROUND_HALF_UP))

def as_num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")

def median_excl_zero(series: pd.Series) -> float:
    vals = as_num(series).replace(0, np.nan)
    return float(vals.median(skipna=True))

# ---------- paths ----------
IN_PATH  = r"E:\Guido\sibel\Masters_Thesis\data\melted_OG_matrix_non_comma_seperated_based_oninventory28.10.tsv"
OUT_PATH = r"E:\Guido\sibel\Masters_Thesis\data\Inventory29102025_Orthogroups_intensity_mean.MEDIAN_NORM.tsv"

# ---------- whitelist ----------
SPECIES_ALLOW = {
    'GCF_904849725.1_Hordeum_vulgare',
    'GCF_902167145.1_Zea_mays',
    'GCA_963924085.1_Cenchrus_americanus.helixer',
    'GCF_034140825.1_Oryza_sativa',
    'GCF_000003195.3_Sorghum_bicolor',
    'GCF_018294505.1_Triticum_aestivum',
}

# ---------- main ----------
def main():
    print("[start] normalization_on_NON_collapsed.py")
    print(f"[env] PYTHONUNBUFFERED={os.environ.get('PYTHONUNBUFFERED')}")
    print(f"[paths] IN_PATH={IN_PATH}")
    print(f"[paths] OUT_PATH={OUT_PATH}")

    # load
    df = pd.read_csv(IN_PATH, sep="\t")
    print(f"[chk] loaded: shape={df.shape}, ncols={len(df.columns)}")

    # species cleanup + whitelist
    df['Species'] = df['Species'].astype(str).str.strip()
    species_values = [sp for sp in df['Species'].unique() if sp in SPECIES_ALLOW]
    unknown = sorted(set(df['Species'].unique()) - SPECIES_ALLOW)
    if unknown:
        print(f"[warn] {len(unknown)} species in file but not whitelisted (skipped): {unknown}")

    print(f"[chk] will normalize across {len(species_values)} species (whitelisted)")

    # intensity columns
    META_EXACT  = {"HOG", "orthogroup", "Species", "ProteinID"}
    META_PREFIX = ("Entry_",)
    intensity_cols = [
        c for c in df.columns
        if ("MaxLFQ" in c) and (c not in META_EXACT) and (not any(c.startswith(p) for p in META_PREFIX))
    ]
    if not intensity_cols:
        raise ValueError("No intensity columns matched '*MaxLFQ*'.")

    print(f"[chk] intensity_cols={len(intensity_cols)}")
    print("[chk] first 10 intensity cols:", intensity_cols[:10])

    species_logs = []

    # per species
    for sp in species_values:
        mask_sp = (df["Species"].astype(str) == sp)
        n_rows_sp = int(mask_sp.sum())
        print(f"\n[chk][species={sp}] rows={n_rows_sp}")
        if n_rows_sp == 0:
            print("[chk] no rows for this species (skipping)")
            species_logs.append((sp, np.nan, {}))
            continue

        # column medians within species
        col_medians = {}
        for c in intensity_cols:
            numeric_mask_c = as_num(df.loc[mask_sp, c]).notna()
            if not numeric_mask_c.any():
                col_medians[c] = np.nan
                continue
            med = median_excl_zero(df.loc[mask_sp, c][numeric_mask_c])
            col_medians[c] = med
            n_num = int(numeric_mask_c.sum())
            preview_vals = as_num(df.loc[mask_sp, c])[numeric_mask_c].head(3).tolist()
            print(f"  [chk][{sp}][col={c}] numeric_cells={n_num}, head3={preview_vals}")
            print(f"  [chk][{sp}][col={c}] median_excl_zero={med}")

        valid_meds = [m for m in col_medians.values() if np.isfinite(m) and m > 0]
        if not valid_meds:
            print(f"  [chk][{sp}] no valid medians; skipping normalization")
            species_logs.append((sp, np.nan, {}))
            continue

        sp_global = float(np.median(valid_meds))
        print(f"  [chk][{sp}] species_global_median={sp_global}")

        per_col_scales = {}
        for c in intensity_cols:
            m = col_medians[c]
            scale = (sp_global / m) if (np.isfinite(m) and m > 0) else 1.0
            per_col_scales[c] = scale
            print(f"  [chk][{sp}][col={c}] scale = {sp_global} / {m} = {scale}")

            if scale != 1.0:
                col_vals = as_num(df.loc[mask_sp, c])
                apply_mask = col_vals.notna() & (col_vals != 0)
                if apply_mask.any():
                    before_vals = col_vals[apply_mask].head(3).tolist()
                    print(f"    [chk][{sp}][col={c}] will apply to {int(apply_mask.sum())} cells; before head3={before_vals}")

                    # apply once
                    df.loc[mask_sp, c] = col_vals.where(~apply_mask, col_vals * scale)

                    # verify
                    post_med = median_excl_zero(df.loc[mask_sp, c][apply_mask])
                    after_vals = as_num(df.loc[mask_sp, c])[apply_mask].head(3).tolist()
                    print(f"    [chk][{sp}][col={c}] after head3={after_vals}")
                    print(f"    [chk][{sp}][col={c}] post-median_excl_zero={post_med}")

        species_logs.append((sp, sp_global, per_col_scales))
        _first = list(per_col_scales.items())[:4]
        print("  [chk] first scales:", ", ".join(f"{k}:×{v:.6g}" for k,v in _first))

    # save
    Path(OUT_PATH).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_PATH, sep="\t", index=False)
    print("\n[chk] quick post-normalization spot-check (first species, first 3 columns)")
    if len(species_values) and len(intensity_cols):
        sp0 = species_values[0]
        mask0 = (df["Species"].astype(str) == sp0)
        for c in intensity_cols[:3]:
            post_med0 = median_excl_zero(df.loc[mask0, c])
            print(f"  [chk][{sp0}][{c}] post-median_excl_zero={post_med0}")

    print(f"[OK] Wrote normalized table → {OUT_PATH}")
    print("\n=== Per-species global medians (first few) ===")
    for sp, sp_global, scales in species_logs[:8]:
        g = "NA" if not np.isfinite(sp_global) else f"{sp_global:.4g}"
        few = list(scales.items())[:4]
        preview = ", ".join(f"{k}:×{v:.6g}" for k, v in few)
        print(f"{sp}: global={g} | {preview}{' ...' if len(scales)>4 else ''}")

# ---------- guarded run (ensures tracebacks go to log) ----------
if __name__ == "__main__":
    try:
        main()
    except Exception:
        print("\n[ERROR] Unhandled exception:", file=sys.stderr)
        traceback.print_exc()  # this goes to both console and log via tee
        raise
