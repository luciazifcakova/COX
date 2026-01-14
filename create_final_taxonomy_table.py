#!/usr/bin/env python3
"""
create_final_taxonomy_table.py

Builds a single final_taxonomy_table.tsv from per-sample outputs:

Inputs per sample:
- 03_taxonomy/<sample>_taxonomy.tsv
- 02_consensus/<sample>_clusters.tsv

Goals:
- Parse BOTH "clean" TSVs (typical NT output) AND "messy" TSVs (some MIDORI outputs can be non-rectangular).
- Attach read_count from clusters.tsv when possible.
- If a row is a collapsed query like:
    consensus_cl_id_0+3_total_supporting_reads_107604+114567
  then derive read_count as sum of supporting reads if clusters join fails.

Output:
- <output_dir>/final_taxonomy_table.tsv

Notes:
- We do NOT collapse rows here (that is handled earlier by collapse_taxonomy_if_same_species.py if desired).
- We keep all rows from each taxonomy.tsv.
"""

import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

# collapsed query format: ...total_supporting_reads_123+456+789
SUPPORT_RE = re.compile(r"total_supporting_reads_([0-9+]+)")


# -------------------------
# Robust TSV parsing
# -------------------------

def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def _try_read_tsv_strict(path: Path) -> pd.DataFrame:
    """
    Fast path: strict tab-separated with a header row.
    """
    return pd.read_csv(path, sep="\t", dtype=str, engine="python")


def _robust_parse_table(path: Path) -> pd.DataFrame:
    """
    Robust fallback:
    - Find the header line containing at least 'query' and 'taxonomy'
    - Split header by tabs OR whitespace
    - Parse following lines, splitting similarly, then pad/truncate to match header length
    This salvages many MIDORI "non-rectangular" tables.
    """
    text = _read_text(path)
    lines = [ln.rstrip("\n") for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        return pd.DataFrame()

    header_idx = None
    header_cols: List[str] = []
    header_sep = None  # "tab" or "ws"

    for i, ln in enumerate(lines):
        # prefer tab header
        if "\t" in ln:
            cols = [c.strip() for c in ln.split("\t")]
            if "query" in cols and "taxonomy" in cols:
                header_idx = i
                header_cols = cols
                header_sep = "tab"
                break

    if header_idx is None:
        # try whitespace header
        for i, ln in enumerate(lines):
            cols = ln.split()
            if "query" in cols and "taxonomy" in cols:
                header_idx = i
                header_cols = cols
                header_sep = "ws"
                break

    if header_idx is None or not header_cols:
        # cannot recover
        return pd.DataFrame()

    def split_row(ln: str) -> List[str]:
        if header_sep == "tab":
            return [c.strip() for c in ln.split("\t")]
        return ln.split()

    rows: List[List[str]] = []
    ncol = len(header_cols)

    for ln in lines[header_idx + 1:]:
        parts = split_row(ln)

        if not parts:
            continue

        # If whitespace-splitting, stitle can blow up into many tokens.
        # We salvage by:
        # - if too many fields, merge extras into the last column
        if len(parts) > ncol:
            merged = parts[:ncol-1] + [" ".join(parts[ncol-1:])]
            parts = merged

        # If too few fields, pad with empty strings
        if len(parts) < ncol:
            parts = parts + [""] * (ncol - len(parts))

        rows.append(parts[:ncol])

    if not rows:
        return pd.DataFrame(columns=header_cols)

    return pd.DataFrame(rows, columns=header_cols)


def safe_read_taxonomy(path: Path) -> Tuple[pd.DataFrame, Optional[str]]:
    """
    Returns (df, error_reason)
    """
    try:
        df = _try_read_tsv_strict(path)
        return df, None
    except Exception:
        try:
            df = _robust_parse_table(path)
            if df.empty:
                return df, "robust_parse_empty"
            return df, None
        except Exception as e:
            return pd.DataFrame(), f"robust_parse_exception:{e}"


# -------------------------
# read_count helpers
# -------------------------

def read_clusters(clusters_path: Path) -> Dict[str, int]:
    """
    clusters.tsv should have:
      consensus_id  read_count
    """
    if not clusters_path.exists():
        return {}

    try:
        df = pd.read_csv(clusters_path, sep="\t", dtype=str)
    except Exception:
        return {}

    if df.empty:
        return {}

    # tolerate missing header
    cols = [c.lower() for c in df.columns.tolist()]
    if "consensus_id" not in cols or "read_count" not in cols:
        # attempt to force first two cols
        if df.shape[1] >= 2:
            df = df.iloc[:, :2].copy()
            df.columns = ["consensus_id", "read_count"]
        else:
            return {}

    if "consensus_id" not in df.columns or "read_count" not in df.columns:
        return {}

    df["read_count"] = pd.to_numeric(df["read_count"], errors="coerce").fillna(0).astype(int)
    return dict(zip(df["consensus_id"].astype(str), df["read_count"].astype(int)))


def readcount_from_collapsed_query(q: str) -> int:
    """
    If query contains total_supporting_reads_123+456 => sum.
    """
    if not isinstance(q, str):
        return 0
    m = SUPPORT_RE.search(q)
    if not m:
        return 0
    nums = m.group(1).split("+")
    total = 0
    for x in nums:
        x = x.strip()
        if x.isdigit():
            total += int(x)
    return total


# -------------------------
# Main merge logic
# -------------------------

def ensure_columns(df: pd.DataFrame, required: List[str]) -> pd.DataFrame:
    for c in required:
        if c not in df.columns:
            df[c] = ""
    return df


def main():
    if len(sys.argv) < 3:
        print("Usage: create_final_taxonomy_table.py <output_dir> <final_output.tsv>", file=sys.stderr)
        sys.exit(1)

    output_dir = Path(sys.argv[1])
    final_output = Path(sys.argv[2])

    taxonomy_dir = output_dir / "03_taxonomy"
    consensus_dir = output_dir / "02_consensus"

    tax_files = sorted(taxonomy_dir.glob("*_taxonomy.tsv"))
    if not tax_files:
        print(f"ERROR: No taxonomy files found in {taxonomy_dir}", file=sys.stderr)
        sys.exit(1)

    required_base = [
        "query", "taxonomy", "confidence", "num_hits",
        "best_pident", "best_evalue", "best_bitscore", "best_accession", "alignment_length",
    ]

    # optional but supported by your newer assign_lca.py
    optional_cols = ["db_used"]

    rank_cols = []
    for r in RANKS:
        rank_cols += [r, f"{r}_taxid"]

    all_rows = []
    skipped: List[Tuple[str, str]] = []

    for tax_file in tax_files:
        sample = tax_file.name.replace("_taxonomy.tsv", "")

        tax_df, err = safe_read_taxonomy(tax_file)
        if tax_df.empty:
            reason = err or "empty"
            skipped.append((tax_file.name, reason))
            continue

        # must have at least query + taxonomy
        if "query" not in tax_df.columns or "taxonomy" not in tax_df.columns:
            skipped.append((tax_file.name, "missing_query_or_taxonomy"))
            continue

        tax_df = ensure_columns(tax_df, required_base)
        tax_df = ensure_columns(tax_df, optional_cols)
        tax_df = ensure_columns(tax_df, rank_cols)

        # attach read_count
        clusters_path = consensus_dir / f"{sample}_clusters.tsv"
        rc_map = read_clusters(clusters_path)

        tax_df["read_count"] = tax_df["query"].map(lambda q: rc_map.get(str(q), None))

        # fill missing read_count via collapsed query name
        missing = tax_df["read_count"].isna()
        if missing.any():
            tax_df.loc[missing, "read_count"] = tax_df.loc[missing, "query"].map(readcount_from_collapsed_query)

        tax_df["read_count"] = pd.to_numeric(tax_df["read_count"], errors="coerce").fillna(0).astype(int)

        # add sample col
        tax_df.insert(0, "sample", sample)

        out_cols = (
            ["sample"]
            + required_base
            + ["read_count"]
            + ["db_used"]   # always included (empty if absent)
            + rank_cols
        )

        tax_df = ensure_columns(tax_df, out_cols)
        tax_df = tax_df[out_cols].copy()

        # numeric conversions (safe)
        for c in ["read_count", "confidence", "best_pident", "best_evalue", "best_bitscore", "alignment_length", "num_hits"]:
            if c in tax_df.columns:
                tax_df[c] = pd.to_numeric(tax_df[c], errors="coerce")

        all_rows.append(tax_df)

    if not all_rows:
        print("ERROR: No rows produced (all taxonomy files failed to parse).", file=sys.stderr)
        if skipped:
            print("Examples of skipped files:", file=sys.stderr)
            for f, r in skipped[:10]:
                print(f"  - {f}: {r}", file=sys.stderr)
        sys.exit(1)

    final_df = pd.concat(all_rows, ignore_index=True)

    # Write output
    final_df.to_csv(final_output, sep="\t", index=False)

    print(f"✓ Created final taxonomy table: {final_output}")
    print(f"  Rows: {len(final_df)}")
    print(f"  Samples: {final_df['sample'].nunique()}")

    if skipped:
        print(f"  WARNING: skipped {len(skipped)} taxonomy.tsv files because they could not be parsed.")
        # show a few reasons
        reason_counts: Dict[str, int] = {}
        for _, r in skipped:
            reason_counts[r] = reason_counts.get(r, 0) + 1
        top = sorted(reason_counts.items(), key=lambda x: (-x[1], x[0]))[:10]
        for r, n in top:
            print(f"    - {r}: {n}")


if __name__ == "__main__":
    main()
