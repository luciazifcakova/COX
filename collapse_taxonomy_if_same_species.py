#!/usr/bin/env python3
"""
Collapse taxonomy.tsv to a single row ONLY IF all assigned consensus sequences
resolve to the same species (species_taxid), and sum read_count across clusters.

Produces query like:
  consensus_cl_id_1+3_total_supporting_reads_456+467

If species_taxid is missing, it will not collapse.
Unassigned rows are never used for collapsing.
"""

import argparse
import re
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--taxonomy_tsv", required=True)
    p.add_argument("--clusters_tsv", required=True)
    p.add_argument("--output_tsv", required=True)
    return p.parse_args()


CL_RE = re.compile(r"cl_id_(\d+)")
SUPPORT_RE = re.compile(r"total_supporting_reads_(\d+)")


def build_collapsed_query(ids, supports):
    return f"consensus_cl_id_{'+'.join(ids)}_total_supporting_reads_{'+'.join(supports)}"


def main():
    args = parse_args()

    tax = pd.read_csv(args.taxonomy_tsv, sep="\t", dtype=str).fillna("")
    if tax.empty:
        tax.to_csv(args.output_tsv, sep="\t", index=False)
        return

    cl = pd.read_csv(args.clusters_tsv, sep="\t", dtype=str).fillna("")
    rc = {}
    if "consensus_id" in cl.columns and "read_count" in cl.columns:
        cl["read_count"] = pd.to_numeric(cl["read_count"], errors="coerce").fillna(0).astype(int)
        rc = dict(zip(cl["consensus_id"], cl["read_count"]))

    if "read_count" not in tax.columns:
        tax["read_count"] = tax["query"].map(lambda q: rc.get(q, 0)).astype(int)

    assigned = tax[(tax["taxonomy"] != "Unassigned") & (tax.get("species_taxid", "") != "")]
    if assigned.empty:
        tax.drop(columns=["read_count"], errors="ignore").to_csv(args.output_tsv, sep="\t", index=False)
        return

    species_taxids = assigned["species_taxid"].unique().tolist()
    if len(species_taxids) != 1:
        tax.drop(columns=["read_count"], errors="ignore").to_csv(args.output_tsv, sep="\t", index=False)
        return

    if assigned["query"].nunique() <= 1:
        tax.drop(columns=["read_count"], errors="ignore").to_csv(args.output_tsv, sep="\t", index=False)
        return

    ids = []
    supports = []
    for q in assigned["query"].tolist():
        m1 = CL_RE.search(q)
        m2 = SUPPORT_RE.search(q)
        ids.append(m1.group(1) if m1 else q)
        supports.append(m2.group(1) if m2 else "0")

    collapsed_query = build_collapsed_query(ids, supports)

    row0 = assigned.iloc[0].copy()
    row0["query"] = collapsed_query
    row0["read_count"] = int(assigned["read_count"].sum())

    # keep best stats
    for col in ["confidence", "best_pident", "best_bitscore", "alignment_length", "num_hits"]:
        if col in assigned.columns:
            row0[col] = pd.to_numeric(assigned[col], errors="coerce").max()
    if "best_evalue" in assigned.columns:
        row0["best_evalue"] = pd.to_numeric(assigned["best_evalue"], errors="coerce").min()

    unassigned = tax[tax["taxonomy"] == "Unassigned"].copy()
    out = pd.concat([pd.DataFrame([row0]), unassigned], ignore_index=True)

    out.drop(columns=["read_count"], errors="ignore").to_csv(args.output_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
