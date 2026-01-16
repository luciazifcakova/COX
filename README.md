# COX

Probable purpose of the study: metabarcoding study of an agroecosystem or monitored natural habitat. DNA was likely extracted from a bulk sample (like a trap catch) to get a snapshot of the local community. The goal could have been to monitor biodiversity, track invasive species, or assess agricultural ecosystem health without needing to visually identify each specimen.

Looking at this list of species, they appear to be connected by their association with plants as pests, predators, or pollinators—making them highly relevant to agriculture or ecosystem monitoring.

| Category | Example Species | Role/Concern |
| :--- | :--- | :--- |
| **Major Crop Pests** | `Ceresa bubalus`, `Aphrophora major`, `Palomena prasina`, `Coreus marginatus` | Direct damage to crops (orchards, vineyards). |
| **Generalist Herbivores** | `Graphosoma italicum`, `Tropidothorax leucopterus`  | Feed on various cultivated plants. |
| **Beneficials / Pollinators** | `Halictus tetrazonianellus`, `Vespula vulgaris` , `Myrmeleon formicarius` | Pollination and/or pest predation. |
| **Disease Vectors** | `Lipoptena fortisetosa` , `Dictyophara europaea` | Potential vector for pathogens. |

Pipeline logic:
run_coi_pipeline.sh (core orchestrator):
Simplified main steps:
1. QC with NanoFilt
2. NGSpeciesID clustering
3. extract_cluster_membership.py - Counts raw reads per consensus
4. BLAST taxonomy assignment (MIDORI2 or fallback to ncbi238 nt)
5. create finale taxonomy table
   
NGSpeciesID — species-first consensus pipeline:
Clusters long reads by similarity, builds species-level consensus sequences, polishes them (Racon / Medaka), outputs one consensus per cluster. Designed specifically for Nanopore barcoding, handles mixed species samples.

use curated, full 16S lenght MIDORI2 database for Cytochrome c oxidase subunit 1 (CO1) (https://onlinelibrary.wiley.com/doi/10.1002/edn3.303), downloaded longest representative sequnces from database in blast format to hpc wget https://www.reference-midori.info/download/Databases/GenBank268_2025-08-14/BLAST/longest/MIDORI2_LONGEST_NUC_GB268_CO1_BLAST.zip

Results - see  final_taxonomy_table file or short table here:

| Sample ID | Percent Identity | Read Count | Class | Order | Family | Genus | Species | Notes |
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
