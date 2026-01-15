# COX

Probable purpose of the study: metabarcoding study of an agroecosystem or monitored natural habitat. The researcher likely extracted DNA from a bulk sample (like a trap catch) to get a snapshot of the local community.

Looking at this list of species, they appear to be connected by their association with plants as pests, predators, or pollinators—making them highly relevant to agriculture or ecosystem monitoring.

Here is a breakdown of the likely connections:

| Category | Example Species | Role/Concern |
| :--- | :--- | :--- |
| **Major Crop Pests** | `Ceresa bubalus`, `Aphrophora major`| Direct damage to crops (orchards, vineyards). |
| **Generalist Herbivores** | `Palomena prasina`, `Coreus marginatus` (Dock bug) | Feed on various cultivated plants. |
| **Beneficials / Pollinators** | `Halictus tetrazonianellus` (Sweat bee), `Vespula vulgaris` (Common wasp) | Pollination and/or pest predation. |
| **Disease Vectors** | `Lipoptena fortisetosa` (Deer ked) | Potential vector for animal pathogens. |
| **Ecosystem Bioindicators** | `Succinea putris` (Amber snail), `Myrmeleon formicarius` (Antlion) | Indicate soil moisture, habitat quality. |

### 🧬 Most Likely Sequencing Scenario
Given these connections, the most probable reason for sequencing them together would be an **environmental DNA (eDNA) metabarcoding study**.

*   **Method**: Scientists collect a single environmental sample (like **soil, leaf washings, or insect trap contents**).
*   **Result**: The DNA of all organisms in that sample—pests, pollinators, predators, decomposers—is sequenced simultaneously.
*   **Goal**: To monitor **biodiversity, track invasive species, or assess agricultural ecosystem health** without needing to visually identify each specimen.

### 🔍 Other Possible Research Contexts
While eDNA is the strongest link, other focused studies could explain this list:

1.  **Agricultural Pest Surveillance**: A study monitoring pest populations and their natural enemies in a specific crop system (e.g., orchards or vineyards).
2.  **Food Web or Trophic Study**: Research on a particular ecosystem (like a meadow or forest edge) examining predator-prey or plant-insect interactions.
3.  **Invasive Species Impact Study**: Investigating how an invasive plant or insect affects the local arthropod community.

**In summary:** This list looks like the output from a **metabarcoding study** of an **agroecosystem or monitored natural habitat**. The researcher likely extracted DNA from a bulk sample (like soil or a malaise trap catch) to get a snapshot of the arthropod community, revealing pests, beneficials, and bioindicators all at once.

Would you like to explore how such a metabarcoding study is typically designed, or how you could analyze this type of species list data further?

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
