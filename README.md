# COX

Probable purpose of the study: metabarcoding study of an agroecosystem or monitored natural habitat. The researcher likely extracted DNA from a bulk sample (like a trap catch) to get a snapshot of the local community. The goal could have been to monitor biodiversity, track invasive species, or assess agricultural ecosystem health without needing to visually identify each specimen.

Looking at this list of species, they appear to be connected by their association with plants as pests, predators, or pollinators—making them highly relevant to agriculture or ecosystem monitoring.

| Category | Example Species | Role/Concern |
| :--- | :--- | :--- |
| **Major Crop Pests** | `Ceresa bubalus`, `Aphrophora major`, `Palomena prasina`, `Coreus marginatus`| Direct damage to crops (orchards, vineyards). |
| **Generalist Herbivores** | `Graphosoma italicum`,  | Feed on various cultivated plants. |
| **Beneficials / Pollinators** | `Halictus tetrazonianellus`, `Vespula vulgaris` , `Myrmeleon formicarius` | Pollination and/or pest predation. |
| **Disease Vectors** | `Lipoptena fortisetosa` , `Dictyophara europaea` | Potential vector for animal pathogens. |


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


use MIDORI2 database for Cytochrome c oxidase subunit 1 (CO1) 
https://onlinelibrary.wiley.com/doi/10.1002/edn3.303

#download longest representative sequnces from database in blast format to hpc
wget https://www.reference-midori.info/download/Databases/GenBank268_2025-08-14/BLAST/longest/MIDORI2_LONGEST_NUC_GB268_CO1_BLAST.zip
