# COX

NGSpeciesID — species-first consensus pipeline

Clusters long reads by similarity
Builds species-level consensus sequences
Polishes them (Racon / Medaka)
Outputs one consensus per species
Key idea  - “Species are the fundamental unit — not reads.”

✔ Designed specifically for Nanopore barcoding
✔ Produces publication-grade consensuses
✔ Handles mixed species samples
✔ Widely used in ONT barcoding literature


Pipeline logic:
run_coi_pipeline.sh (core orchestrator):
Simplified main steps:
1. QC with NanoFilt
2. NGSpeciesID clustering
3. extract_cluster_membership.py - Counts ACTUAL raw reads per consensus
4. BLAST taxonomy assignment (MIDORI2 or fallback to ncbi nt)
5. create finale taxonomy table
