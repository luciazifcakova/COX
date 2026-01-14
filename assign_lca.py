#!/usr/bin/env python3
"""
assign_lca.py  (best-hit taxonomy assignment for BLAST outfmt 6)

Supports BOTH:
1) MIDORI (taxonomy embedded in sseqid after ###)
2) NT (taxonomy from staxids using NCBI taxdump)

Key behavior:
- Reads BLAST outfmt 6 tables (12 cols classic or 15 cols with staxids/sscinames/stitle).
- Filters hits by min_identity, max_evalue, min_bitscore.
- Keeps ALL filtered hits optionally (--hits_output).
- Selects best hit per query by pident desc, bitscore desc, evalue asc.
- If best hit resolves to generic species like "Genus sp." then falls back to next hit
  that resolves to species-specific name (still filtered).

Output taxonomy.tsv columns:
query, db_used, taxonomy, confidence, num_hits, best_pident, best_evalue, best_bitscore,
best_accession, alignment_length, kingdom..species + *_taxid
"""

import argparse
import re
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import pandas as pd

RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


def parse_args():
    p = argparse.ArgumentParser(description="Assign BEST HIT taxonomy from BLAST results")
    p.add_argument("--input", required=True, help="BLAST results TSV (outfmt 6)")
    p.add_argument("--output", required=True, help="Output taxonomy TSV (one row per query)")
    p.add_argument("--hits_output", default="", help="Optional: write per-hit table with resolved taxonomy")
    p.add_argument("--db_used", default="", help="Database label written to output (MIDORI or NT)")
    p.add_argument("--min_hits", type=int, default=1, help="Minimum hits per query after filtering")
    p.add_argument("--min_identity", type=float, default=0.0, help="Minimum percent identity")
    p.add_argument("--max_evalue", type=float, default=1e-5, help="Maximum E-value")
    p.add_argument("--min_bitscore", type=float, default=0.0, help="Minimum bitscore (strictly >)")
    p.add_argument("--taxdump_dir", default="", help="NCBI taxdump dir (nodes.dmp, names.dmp, merged.dmp)")
    return p.parse_args()


# -----------------------------
# BLAST parsing (MIDORI + NT)
# -----------------------------

def _read_blast_table(path: str) -> pd.DataFrame:
    """
    Try strict TSV first. If it has fewer than 12 columns, fall back to whitespace buffering
    (like your old working script).
    """
    try:
        df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str, engine="python")
        if df.shape[1] >= 12:
            return df
    except Exception:
        pass

    # fallback: whitespace-buffered parsing (old script style)
    records = []
    buf: List[str] = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            buf.extend(line.strip().split())
            if len(buf) >= 12:
                records.append(buf)
                buf = []

    if records:
        maxc = max(len(r) for r in records)
        norm = [r + [""] * (maxc - len(r)) for r in records]
        return pd.DataFrame(norm)

    return pd.DataFrame()


def _assign_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Tolerate:
    - 12 cols: classic outfmt 6
    - 13 cols: classic + stitle merged
    - >=15 cols: expected modern outfmt with staxids/sscinames/stitle
      qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames stitle
    """
    if df.empty:
        return df

    base12 = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ]
    n = df.shape[1]

    if n == 12:
        df.columns = base12
        return df

    if n >= 15:
        # if more than 15, merge tail into stitle (last column)
        if n > 15:
            df[14] = df.iloc[:, 14:].astype(str).agg(" ".join, axis=1)
            df = df.iloc[:, :15]
        df.columns = base12 + ["staxids", "sscinames", "stitle"]
        return df

    # 13 or 14: merge tail to stitle
    df[12] = df.iloc[:, 12:].astype(str).agg(" ".join, axis=1)
    df = df.iloc[:, :13]
    df.columns = base12 + ["stitle"]
    return df


def parse_blast_results(path: str) -> pd.DataFrame:
    df = _read_blast_table(path)
    df = _assign_columns(df)
    for c in ["pident", "length", "evalue", "bitscore"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


# -----------------------------
# MIDORI taxonomy parsing
# -----------------------------

def extract_midori_taxonomy_from_sseqid(sseqid: str) -> str:
    if not isinstance(sseqid, str):
        return ""
    if "###" not in sseqid:
        return ""
    return sseqid.split("###", 1)[1].strip()


def split_tax_levels(tax: str) -> List[str]:
    if not isinstance(tax, str):
        return []
    tax = tax.strip()
    if not tax:
        return []
    levels = [t.strip() for t in tax.split(";") if t.strip()]
    # drop root_...
    if levels and levels[0].lower().startswith("root"):
        levels = levels[1:]
    return levels


def split_name_taxid(level: str) -> Tuple[str, str]:
    if not isinstance(level, str):
        return ("", "")
    level = level.strip()
    if not level:
        return ("", "")
    m = re.match(r"^(.*)_(\d+)$", level)
    if m:
        return (m.group(1), m.group(2))
    return (level, "")


def levels_to_rank_columns(levels: List[str]) -> Dict[str, str]:
    out = {r: "" for r in RANKS}
    out.update({f"{r}_taxid": "" for r in RANKS})
    for i, level in enumerate(levels[:len(RANKS)]):
        name, taxid = split_name_taxid(level)
        r = RANKS[i]
        out[r] = name
        out[f"{r}_taxid"] = taxid
    return out


def taxonomy_string_from_rank_cols(rank_cols: Dict[str, str]) -> str:
    parts = []
    for r in RANKS:
        v = (rank_cols.get(r) or "").strip()
        if v:
            parts.append(v)
    return ";".join(parts) if parts else "Unassigned"


def is_generic_species_name(species_name: str) -> bool:
    if not species_name:
        return True
    s = str(species_name).strip()
    # "sp." token anywhere
    return bool(re.search(r"(^|\s)sp\.?($|\s)", s))


# -----------------------------
# Taxdump support for NT
# -----------------------------

class Taxdump:
    def __init__(self, taxdump_dir: str):
        self.taxdump_dir = Path(taxdump_dir) if taxdump_dir else None
        self.parent: Dict[str, str] = {}
        self.rank: Dict[str, str] = {}
        self.name: Dict[str, str] = {}
        self.merged: Dict[str, str] = {}

        if not self.taxdump_dir:
            return

        nodes = self.taxdump_dir / "nodes.dmp"
        names = self.taxdump_dir / "names.dmp"
        merged = self.taxdump_dir / "merged.dmp"

        if nodes.exists():
            self._load_nodes(nodes)
        if names.exists():
            self._load_names(names)
        if merged.exists():
            self._load_merged(merged)

    def ok(self) -> bool:
        return bool(self.parent) and bool(self.rank) and bool(self.name)

    def _load_nodes(self, path: Path):
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                parts = [p.strip() for p in line.split("|")]
                if len(parts) < 3:
                    continue
                tax_id, parent_id, rank = parts[0], parts[1], parts[2]
                self.parent[tax_id] = parent_id
                self.rank[tax_id] = rank

    def _load_names(self, path: Path):
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                parts = [p.strip() for p in line.split("|")]
                if len(parts) < 4:
                    continue
                tax_id, name_txt, _, name_class = parts[0], parts[1], parts[2], parts[3]
                if name_class == "scientific name":
                    self.name[tax_id] = name_txt

    def _load_merged(self, path: Path):
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                parts = [p.strip() for p in line.split("|")]
                if len(parts) < 2:
                    continue
                old_id, new_id = parts[0], parts[1]
                self.merged[old_id] = new_id

    def normalize_taxid(self, taxid: str) -> str:
        if not taxid:
            return ""
        t = taxid.strip()
        for _ in range(10):
            if t in self.merged:
                t = self.merged[t]
            else:
                break
        return t

    @lru_cache(maxsize=200000)
    def lineage_rank_map(self, taxid: str) -> Dict[str, Tuple[str, str]]:
        """
        Map our ranks to (name, taxid) by walking parent links.
        Also maps NCBI 'superkingdom' -> our 'kingdom' if kingdom missing.
        """
        out = {r: ("", "") for r in RANKS}
        if not self.ok():
            return out

        t = self.normalize_taxid(taxid)
        if not t or t not in self.parent:
            return out

        cur = t
        visited = set()
        while cur and cur not in visited:
            visited.add(cur)
            r = self.rank.get(cur, "")
            if r in RANKS and out[r][1] == "":
                out[r] = (self.name.get(cur, ""), cur)

            parent = self.parent.get(cur, "")
            if not parent or parent == cur:
                break
            cur = parent

        # fill kingdom from superkingdom if needed
        if out["kingdom"][1] == "":
            cur = t
            visited = set()
            while cur and cur not in visited:
                visited.add(cur)
                if self.rank.get(cur, "") == "superkingdom":
                    out["kingdom"] = (self.name.get(cur, ""), cur)
                    break
                parent = self.parent.get(cur, "")
                if not parent or parent == cur:
                    break
                cur = parent

        return out


def first_taxid_from_staxids(staxids_val: str) -> str:
    if not isinstance(staxids_val, str):
        return ""
    s = staxids_val.strip()
    if not s:
        return ""
    for tok in re.split(r"[;,]", s):
        tok = tok.strip()
        if tok.isdigit():
            return tok
    return ""


# -----------------------------
# Resolve taxonomy for a hit
# -----------------------------

def resolve_rank_cols_for_row(row: pd.Series, db_used: str, taxdump: Taxdump) -> Tuple[Dict[str, str], str, str]:
    """
    Returns: (rank_cols, taxonomy_str, method)
    method: "midori" | "taxdump" | "title"
    """
    db = (db_used or "").strip().upper()

    # If MIDORI is used, prefer parsing from MIDORI header
    if db == "MIDORI":
        midori = extract_midori_taxonomy_from_sseqid(row.get("sseqid", ""))
        if midori:
            levels = split_tax_levels(midori)
            rank_cols = levels_to_rank_columns(levels)
            return rank_cols, taxonomy_string_from_rank_cols(rank_cols), "midori"
        # If MIDORI db but no ###, fall through to title

    # If NT (or anything else), prefer taxdump if available
    tid = first_taxid_from_staxids(row.get("staxids", ""))
    if tid and taxdump and taxdump.ok():
        rm = taxdump.lineage_rank_map(tid)
        rank_cols = {r: "" for r in RANKS}
        rank_cols.update({f"{r}_taxid": "" for r in RANKS})
        for r in RANKS:
            nm, tx = rm.get(r, ("", ""))
            rank_cols[r] = nm or ""
            rank_cols[f"{r}_taxid"] = tx or ""
        return rank_cols, taxonomy_string_from_rank_cols(rank_cols), "taxdump"

    # Fallback: title/sseqid
    rank_cols = {r: "" for r in RANKS}
    rank_cols.update({f"{r}_taxid": "" for r in RANKS})
    title = str(row.get("stitle", "") or "").strip()
    if title:
        return rank_cols, title, "title"
    sseqid = str(row.get("sseqid", "") or "").strip()
    return rank_cols, (sseqid if sseqid else "Unassigned"), "title"


def pick_best_row_with_sp_skip(group: pd.DataFrame, db_used: str, taxdump: Taxdump):
    """
    Sort by (pident desc, bitscore desc, evalue asc)
    Then choose first candidate with a non-generic species name (when rank cols provide species).
    If none, return top row.
    """
    g = group.copy()

    sort_cols = []
    ascending = []
    if "pident" in g.columns:
        sort_cols.append("pident"); ascending.append(False)
    if "bitscore" in g.columns:
        sort_cols.append("bitscore"); ascending.append(False)
    if "evalue" in g.columns:
        sort_cols.append("evalue"); ascending.append(True)

    if sort_cols:
        g = g.sort_values(sort_cols, ascending=ascending, na_position="last").reset_index(drop=True)

    for i in range(len(g)):
        row = g.iloc[i]
        rank_cols, taxonomy_str, _ = resolve_rank_cols_for_row(row, db_used, taxdump)
        sp_name = (rank_cols.get("species") or "").strip()
        if sp_name and not is_generic_species_name(sp_name):
            return row, rank_cols, taxonomy_str

    top = g.iloc[0]
    rank_cols, taxonomy_str, _ = resolve_rank_cols_for_row(top, db_used, taxdump)
    return top, rank_cols, taxonomy_str


def main():
    args = parse_args()
    df = parse_blast_results(args.input)

    out_cols = [
        "query", "db_used", "taxonomy", "confidence", "num_hits",
        "best_pident", "best_evalue", "best_bitscore", "best_accession", "alignment_length"
    ]
    for r in RANKS:
        out_cols += [r, f"{r}_taxid"]

    # Always write a valid TSV even if empty
    if df.empty or "qseqid" not in df.columns:
        pd.DataFrame(columns=out_cols).to_csv(args.output, sep="\t", index=False)
        if args.hits_output:
            pd.DataFrame().to_csv(args.hits_output, sep="\t", index=False)
        return

    # Filters
    if "pident" in df.columns and args.min_identity > 0:
        df = df[df["pident"].notna() & (df["pident"] >= args.min_identity)]
    if "evalue" in df.columns and args.max_evalue is not None:
        df = df[df["evalue"].notna() & (df["evalue"] <= args.max_evalue)]
    if "bitscore" in df.columns and args.min_bitscore > 0:
        df = df[df["bitscore"].notna() & (df["bitscore"] > args.min_bitscore)]

    taxdump = Taxdump(args.taxdump_dir)

    # Optional: per-hit output
    if args.hits_output:
        hit_rows = []
        for _, row in df.iterrows():
            rank_cols, taxonomy_str, method = resolve_rank_cols_for_row(row, args.db_used, taxdump)
            hit_rows.append({
                "query": row.get("qseqid", ""),
                "db_used": args.db_used,
                "method": method,
                "pident": row.get("pident", ""),
                "evalue": row.get("evalue", ""),
                "bitscore": row.get("bitscore", ""),
                "alignment_length": row.get("length", ""),
                "accession": row.get("sseqid", ""),
                "staxids": row.get("staxids", ""),
                "sscinames": row.get("sscinames", ""),
                "stitle": row.get("stitle", ""),
                "taxonomy": taxonomy_str,
                **{r: rank_cols.get(r, "") for r in RANKS},
                **{f"{r}_taxid": rank_cols.get(f"{r}_taxid", "") for r in RANKS},
            })
        pd.DataFrame(hit_rows).to_csv(args.hits_output, sep="\t", index=False)

    # Per-query best hits
    results = []
    for qseqid, group in df.groupby("qseqid", sort=False):
        if len(group) < args.min_hits:
            continue

        best_row, rank_cols, taxonomy_str = pick_best_row_with_sp_skip(group, args.db_used, taxdump)

        conf = 0.0
        if pd.notna(best_row.get("pident", None)):
            conf = min(float(best_row["pident"]) / 100.0, 1.0)

        row_out = {
            "query": qseqid,
            "db_used": args.db_used,
            "taxonomy": taxonomy_str if taxonomy_str else "Unassigned",
            "confidence": float(conf),
            "num_hits": int(len(group)),
            "best_pident": float(best_row["pident"]) if pd.notna(best_row.get("pident", None)) else "",
            "best_evalue": float(best_row["evalue"]) if pd.notna(best_row.get("evalue", None)) else "",
            "best_bitscore": float(best_row["bitscore"]) if pd.notna(best_row.get("bitscore", None)) else "",
            "best_accession": str(best_row.get("sseqid", "") or ""),
            "alignment_length": int(best_row["length"]) if pd.notna(best_row.get("length", None)) else "",
        }
        for r in RANKS:
            row_out[r] = rank_cols.get(r, "")
            row_out[f"{r}_taxid"] = rank_cols.get(f"{r}_taxid", "")
        results.append(row_out)

    out_df = pd.DataFrame(results)
    for col in out_cols:
        if col not in out_df.columns:
            out_df[col] = ""
    out_df = out_df.reindex(columns=out_cols)
    out_df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
