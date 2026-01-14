#!/usr/bin/env python3
"""
Extract cluster membership (read counts per consensus) from NGSpeciesID output.

Produces TSV:
consensus_id    read_count

consensus_id is FASTA record id (first token) so it matches BLAST qseqid.
"""

import argparse
import gzip
import os
import re
import glob
import pandas as pd


def parse_arguments():
    p = argparse.ArgumentParser(description="Extract cluster membership from NGSpeciesID output")
    p.add_argument("--input_dir", required=True, help="NGSpeciesID output directory")
    p.add_argument("--output", required=True, help="Output TSV file")
    return p.parse_args()


def open_text_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8", errors="replace")


def count_fastq_reads(fq_path) -> int:
    n_lines = 0
    with open_text_maybe_gz(fq_path) as fh:
        for _ in fh:
            n_lines += 1
    return n_lines // 4


def find_consensus_fastas(input_dir):
    candidates = []
    p1 = os.path.join(input_dir, "consensus_references.fasta")
    if os.path.exists(p1):
        candidates.append(p1)

    candidates.extend(sorted(glob.glob(os.path.join(input_dir, "consensus_reference_*.fasta"))))
    candidates.extend(sorted(glob.glob(os.path.join(input_dir, "consensus_reference_*.fa"))))
    candidates.extend(sorted(glob.glob(os.path.join(input_dir, "consensus_reference_*.fna"))))

    # de-dup preserve order
    seen = set()
    uniq = []
    for c in candidates:
        if c not in seen:
            seen.add(c)
            uniq.append(c)
    return uniq


def parse_consensus_headers_for_counts(fasta_path):
    results = []
    pat = re.compile(r"(?:total_supporting_reads_|supporting_reads_)(\d+)", re.IGNORECASE)

    with open_text_maybe_gz(fasta_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            rec_id = header.split()[0]
            m = pat.search(header)
            read_count = int(m.group(1)) if m else None
            results.append({"consensus_id": rec_id, "read_count": read_count})

    return results


def fallback_counts_from_reads_to_consensus(input_dir, consensus_ids):
    fq_files = sorted(glob.glob(os.path.join(input_dir, "reads_to_consensus_*.fastq*")))
    counts_by_n = {}

    for fq in fq_files:
        base = os.path.basename(fq)
        m = re.match(r"reads_to_consensus_(\d+)\.fastq(\.gz)?$", base)
        if not m:
            continue
        n = m.group(1)
        counts_by_n[n] = count_fastq_reads(fq)

    if not counts_by_n:
        return {}

    mapped = {}
    for cid in consensus_ids:
        m = re.search(r"cl_id_(\d+)", cid)
        if m:
            n = m.group(1)
            if n in counts_by_n:
                mapped[cid] = counts_by_n[n]

    return mapped


def main():
    args = parse_arguments()
    input_dir = args.input_dir

    consensus_fastas = find_consensus_fastas(input_dir)
    if not consensus_fastas:
        pd.DataFrame(columns=["consensus_id", "read_count"]).to_csv(args.output, sep="\t", index=False)
        return

    rows = []
    for fasta in consensus_fastas:
        rows.extend(parse_consensus_headers_for_counts(fasta))

    if not rows:
        pd.DataFrame(columns=["consensus_id", "read_count"]).to_csv(args.output, sep="\t", index=False)
        return

    df = pd.DataFrame(rows)

    if df["consensus_id"].duplicated().any():
        df["_has"] = df["read_count"].notna().astype(int)
        df = (
            df.sort_values(["consensus_id", "_has"], ascending=[True, False])
              .drop_duplicates("consensus_id", keep="first")
              .drop(columns=["_has"])
        )

    missing = df["read_count"].isna().sum()
    if missing > 0:
        fallback_map = fallback_counts_from_reads_to_consensus(input_dir, df["consensus_id"].tolist())
        if fallback_map:
            df.loc[df["read_count"].isna(), "read_count"] = df.loc[df["read_count"].isna(), "consensus_id"].map(fallback_map)

    df["read_count"] = df["read_count"].fillna(0).astype(int)
    df[["consensus_id", "read_count"]].to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
