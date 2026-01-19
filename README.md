# COX, COI, CO1
Cytochrome c oxidase subunit I is a mitochondrial protein-coding marker widely used for DNA barcoding because it provides high species-level resolution across most metazoans due to its balance of conserved priming sites and rapidly evolving regions; however, its use is limited in some metazoa by incomplete reference databases, mitochondrial introgression and nuclear pseudogenes. Although COI is primarily used as an animal mitochondrial marker, homologs of cytochrome c oxidase subunit I are also present in intracellular bacteria such as Rickettsia because mitochondria originated from an alphaproteobacterial ancestor, and these bacteria retain a functional respiratory chain that includes cytochrome c oxidase for oxidative phosphorylation within the host cell.

## Probable purpose of the study: 

This study seems like metabarcoding study of an agroecosystem or monitored natural habitat. DNA was likely extracted from a bulk sample (like a trap catch) to get a snapshot of the local community. The goal could have been to monitor biodiversity, track invasive species, or assess agricultural ecosystem health without needing to visually identify each specimen.

Looking at this list of species, they appear to be connected by their association with plants as pests, predators, or pollinators—making them highly relevant to agricultural ecosystem monitoring.

| Category | Example Species | Role/Concern |
| :--- | :--- | :--- |
| **Major Crop Pests** | `Ceresa bubalus`, `Aphrophora major`, `Palomena prasina`, `Coreus marginatus` | Direct damage to crops (orchards, vineyards). |
| **Generalist Herbivores** | `Graphosoma italicum`, `Tropidothorax leucopterus`  | Feed on various cultivated plants. |
| **Beneficials / Pollinators** | `Halictus tetrazonianellus`, `Vespula vulgaris` , `Myrmeleon formicarius` | Pollination and/or pest predation. |
| **Disease Vectors** | `Lipoptena fortisetosa` , `Dictyophara europaea` | Potential vector for pathogens. |

## Pipeline logic:

run_coi_pipeline.sh (core orchestrator):
Simplified main steps:
1. QC with NanoFilt, NanoPlot, MultiQC
2. NGSpeciesID clustering
3. BLAST taxonomy assignment (MIDORI2 or fallback to ncbi238 nt), BLAST filters: pident>=95, evalue<=1e-25, max_targets=10, bitscore>400, if no species identified either way, fall back to lowest common ancestor (LCA)
4. extract_cluster_membership.py - extract NGSpeciesID cluster membership (read counts per consensus) 
5. assign_lca.py - reads BLAST output form Midori2 or ncbi/nt238, filters hits by min_identity, max_evalue, min_bitscore, selects best hit per query by pident desc, bitscore desc, evalue.
4. 


The workflow is fully containerized with pinned software versions, enabling deterministic reruns across HPC and local environments. All analytical steps are executed via Slurm-compatible scripts, supporting large scale processing while ensuring reproducibility and traceability of results. There is a Slurm-native execution model with array jobs for per-sample processing and dedicated merge job for cohort-level summaries (MultiQC, final tables). 
Because COI datasets are typically PCR-amplified and derived from bulk mixed-organism samples, read counts supporting each NGSpeciesID consensus are interpreted as a sequencing/PCR signal rather than direct organism abundance. Taxonomic calls were assigned by filtered BLAST hits (pident ≥95, evalue ≤1e−25, bitscore >400) with LCA reporting when multiple high-scoring hits were not same, and species-level labels are treated as high-confidence only when identity and hit specificity support unambiguous assignment.

NGSpeciesID https://github.com/ksahlin/NGSpeciesID species-first consensus pipeline was used as it clusters long reads by similarity, builds species-level consensus sequences, polishes them (Medaka), outputs one consensus per cluster. Designed specifically for Nanopore barcoding, handles mixed species samples. NGSpeciesID reduces random sequencing errors and mitigates Nanopore-specific noise before taxonomic assignment.

I have used curated, full 16S lenght MIDORI2 database for Cytochrome c oxidase subunit 1 (CO1) (https://onlinelibrary.wiley.com/doi/10.1002/edn3.303). Database of longest representative sequnces was downloaded from https://www.reference-midori.info/download/Databases/GenBank268_2025-08-14/BLAST/longest/MIDORI2_LONGEST_NUC_GB268_CO1_BLAST.zip, wich reduces partial-hit ambiguity typical of short COI fragments. As a fallback was used taxified ncbi-blast nt version 238, already rpesent on hpc. 


## Limitations:
COI cannot fully resolve cases of mitochondrial introgression or cryptic species without nuclear markers. NUMTs (Nuclear Mitochondrial DNA Segments - fragments of mitochondrial DNA inserted into the cell's nuclear genome) and chimeras were not explicitly filtered via translation-based screening in this run (recommended as a future enhancement).


## Results:

Read length distributions showed a tight peak corresponding to the expected COI amplicon length, indicating successful primer amplification and minimal off-target products. Quality score distributions were consistent across samples and within acceptable ranges for Nanopore amplicon sequencing. Cluster structure (NGSpeciesID) showed strong dominance of a single consensus cluster per sample, as expected for samples dominated by one organism, with minor secondary clusters likely representing natural within species diversity. Read support per consensus was used as a confidence measure (supporting evidence), not as a proxy for organism abundance. QC outputs were preserved per sample and aggregated via MultiQC, providing both sample-level diagnostics and cohort-level overview suitable for production reporting. Fallback strategy via blast nt was implemented but rarely needed, indicating good coverage of the target taxa in MIDORI2.

See  final_taxonomy_table file for full results or short table here:

| Sample ID | Percent Identity | Read count supporting a consensus cluster| Class | Order | Family | Genus | Species | Notes |
|:---|:---:|:---:|:---|:---|:---|:---|:---|:---|
| **S3155_001_1_3_Filtered** | 98.9% | 215,139 | Insecta | Orthoptera | Tettigoniidae | *Phaneroptera* | *Phaneroptera* sp. MAA-2007 | <span style="color: #d9534f;">Pest</span> (ID as *P. nana* via MIDORI2) |
| **S3155_002_2_4_Filtered** | 99.4% | 225,818 | Insecta | Hemiptera | Pentatomidae | *Graphosoma* | *Graphosoma italicum* | Sap feeder |
| **S3155_003_3_6_Filtered** | 99.4% | 248,268 | Insecta | Hemiptera | Dictyopharidae | *Dictyophara* | *Dictyophara europaea* | <span style="color: #d9534f;">Potential vector of phytoplasmas (e.g., Flavescence dorée)</span> |
| **S3155_004_4_9_Filtered** | 99.1% | 162,551 | Insecta | Hemiptera | Aphrophoridae | *Aphrophora* | *Aphrophora major* | <span style="color: #d9534f;">Vector of fungus Xylella fastidiosa</span> |
| **S3155_005_5_10_Filtered** | 99.8% | 164,673 | Insecta | Diptera | Hippoboscidae | *Lipoptena* | *Lipoptena fortisetosa* | <span style="color: #d9534f;">Disease vector</span> |
| **S3155_006_6_11_Filtered** | 99.7% | 138,044 | Insecta | Diptera | Tachinidae | *Ectophasia* | *Ectophasia crassipennis* | <span style="color: #5cb85c;">Beneficial</span> (pest control) |
| **S3155_007_7_13_Filtered** | 99.8% | 161,892 | Insecta | Hemiptera | Nabidae | *Nabis* | *Nabis pseudoferus* | <span style="color: #5cb85c;">Beneficial</span> (pest control) |
| **S3155_008_8_14_Filtered** | 100% | 1,980 | Insecta | Orthoptera | Trigonidiidae | *Nemobius* | *Nemobius sylvestris* | <span style="color: #5cb85c;">Beneficial</span> |
| **S3155_009_9_15_Filtered** | 99.8% | 238,813 | Arachnida | Trombidiformes | Trombidiidae | *Allothrombium* | *Allothrombium fuliginosum* | <span style="color: #5cb85c;">Beneficial</span> (pest control) |
| **S3155_010_10_16_Filtered** | 98.1% | 259,004 | Insecta | Hemiptera | Pyrrhocoridae | *Pyrrhocoris* | *Pyrrhocoris apterus* | Not a pest |
| **S3155_011_11_17_Filtered** | 99.4% | 148,121 | Insecta | Hemiptera | Membracidae | *Ceresa* | *Ceresa bubalus* | <span style="color: #d9534f;">Pest</span> |
| **S3155_012_12_18_Filtered** | 99.8% | 222,171 | Insecta | Neuroptera | Myrmeleontidae | *Myrmeleon* | *Myrmeleon formicarius* | <span style="color: #5cb85c;">Beneficial</span> (pest control) |
| **S3155_013_13_19_Filtered** | 98.6% | 244,797 | Insecta | Hemiptera | Lygaeidae | *Tropidothorax* | *Tropidothorax leucopterus* | Sap feeding |
| **S3155_014_14_21_Filtered** | 100% | 97,348 | Insecta | Hemiptera | Coreidae | *Coreus* | *Coreus marginatus* | <span style="color: #d9534f;">Pest</span> |
| **S3155_015_15_22_Filtered** | 98.8% | 270,061 | Insecta | Hymenoptera | Vespidae | *Vespula* | *Vespula vulgaris* | Not a pest |
| **S3155_016_16_23_Filtered** | 98.9% | 215,679 | Insecta | Hemiptera | Pentatomidae | *Palomena* | *Palomena prasina* | <span style="color: #d9534f;">Pest</span> |
| **S3155_017_17_25_Filtered** | 98.4% | 105,132 | Insecta | Hymenoptera | Halictidae | *Halictus* | *Halictus tetrazonianellus* | <span style="color: #5cb85c;">Beneficial pollinator</span> |
| **S3155_018_18_26_Filtered** | 98.3% | 259,546 | Insecta | Hymenoptera | Vespidae | *Vespa* | *Vespa crabro* | <span style="color: #d9534f;">Pest</span> |
| **S3155_019_19_27_Filtered** | 98.7% | 215,166 | Arachnida | Araneae | Araneidae | *Araneus* | *Araneus quadratus* | <span style="color: #5cb85c;">Beneficial</span> (pest control) |
| **S3155_020_20_28_Filtered** | 99.0% | 234,397 | Insecta | Hymenoptera | Vespidae | *Polistes* | *Polistes dominula* | <span style="color: #d9534f;">Pest</span> |
| **S3155_021_21_29_Filtered** | 98.0% | 235,690 | Arachnida | Opiliones | Phalangiidae | *Phalangium* | *Phalangium opilio* | <span style="color: #5cb85c;">Beneficial</span> |
| **S3155_022_22_30_Filtered** | 97.0% | 232,867 | Insecta | Orthoptera | Acrididae | *Oedipoda* | *Oedipoda caerulescens* | <span style="color: #d9534f;">Pest</span> |
| **S3155_023_23_41_Filtered** | 99.7% | 251,326 | Diplopoda | Julida | Julidae | *Megaphyllum* | *Megaphyllum unilineatum* | Not a pest |
| **S3155_024_24_42_Filtered** | 100% | 116,645 | Gastropoda | Stylommatophora | Succineidae | *Succinea* | *Succinea putris* | Possible pest |
| **S3155_PK_COX_Filtered** | 97.7% | 256,275 | Aves | Charadriiformes | Charadriidae | *Charadrius* | *Charadrius hiaticula* | Positive control |

## Suggestions for further analyses:

By using additional analyses, can we turn nanopore metabarcoding into recurring, high-margin revenue for the company?
We can convert species detections into a Pest Risk Index, where weight species by economic damage, outbreak likelihood, regulatory relevance, vector status, invasivness. Beneficial vs pest balance metrics - using thing like ratios  of predator + parasitoid / herbivore, pollinator presence index, biocontrol capacity score. 
Multiple sampling time points for trend and change detection can turn one-off sequencing into subscription monitoring. Pathogen and symbiont screening using the same samples to get extra information that can be presented as different dataset, host–plant interaction inference using archive DNA approach to answer "what are pests actually feeding on?". 
Finally, these analyses enable the creation of long-term monitoring programs, client-specific data-viewing dashboards, predictive risk frameworks (If pest X appears at abundance Y what is the risk Z in T weeks of it damaging the crop W?), transforming single sequencing projects into durable, high-value service contracts.
