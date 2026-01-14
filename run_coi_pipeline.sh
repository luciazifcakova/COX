#!/bin/bash
# COI taxonomy-only pipeline: QC → optional primer trim → NGSpeciesID → BLAST → assign_best_hit → collapse_if_same_species → final table
set -euo pipefail

START_TIME=$(date +%s)

# -------------------------
# Configuration
# -------------------------
INPUT_DIR="${INPUT_DIR:-/data/input}"
OUTPUT_DIR="${OUTPUT_DIR:-/data/output}"
DB_DIR="${DB_DIR:-/data/databases}"
THREADS="${THREADS:-4}"

SAMPLE_FILE="${SAMPLE_FILE:-}"                   # optional: process only this fastq(.gz)
SAMPLE_NAME_OVERRIDE="${SAMPLE_NAME_OVERRIDE:-}" # optional fixed sample name
MERGE_ONLY="${MERGE_ONLY:-false}"                # skip per-sample; only final merge
DO_FINAL_TABLE="${DO_FINAL_TABLE:-true}"

REMOVE_PRIMERS="${REMOVE_PRIMERS:-false}"
PRIMER_FILE="${INPUT_DIR}/primers.fasta"

MIN_IDENT="${MIN_IDENT:-80.0}"
MAX_EVALUE="${MAX_EVALUE:-1e-10}"
MAX_TARGETS="${MAX_TARGETS:-10}"
MIN_BITSCORE="${MIN_BITSCORE:-400}"              # taxonomy filter

# nt fallback DB prefix (must be visible inside container if you want fallback)
# Example: /data/nt/nt
NT_PREFIX="${NT_PREFIX:-}"

# NCBI taxdump dir (must be visible inside container for nt fallback taxonomy)
# Example: /data/taxdump  (nodes.dmp, names.dmp, merged.dmp required)
TAXDUMP_DIR="${TAXDUMP_DIR:-}"

echo "========================================="
echo "COI TAXONOMY PIPELINE"
echo "Input dir:     ${INPUT_DIR}"
echo "Output dir:    ${OUTPUT_DIR}"
echo "DB dir:        ${DB_DIR}"
echo "Threads:       ${THREADS}"
echo "Merge only:    ${MERGE_ONLY}"
echo "Final table:   ${DO_FINAL_TABLE}"
echo "Primer trim:   ${REMOVE_PRIMERS}"
echo "BLAST filters: pident>=${MIN_IDENT}, evalue<=${MAX_EVALUE}, max_targets=${MAX_TARGETS}, bitscore>${MIN_BITSCORE}"
echo "nt fallback:   ${NT_PREFIX:-<disabled>}"
echo "taxdump dir:   ${TAXDUMP_DIR:-<disabled>}"
echo "========================================="

mkdir -p "${OUTPUT_DIR}"/{00_logs,01_qc,02_consensus,03_taxonomy}

pick_db_prefix() {
  local prefix=""
  if [ -f "${DB_DIR}/MIDORI2_LONGEST_NUC_GB268_CO1_BLAST2.nsq" ] || [ -f "${DB_DIR}/MIDORI2_LONGEST_NUC_GB268_CO1_BLAST2.nal" ]; then
    prefix="${DB_DIR}/MIDORI2_LONGEST_NUC_GB268_CO1_BLAST2"
  elif [ -f "${DB_DIR}/MIDORI2_UNIQ_NUC_GB268_CO1_BLAST2.nsq" ] || [ -f "${DB_DIR}/MIDORI2_UNIQ_NUC_GB268_CO1_BLAST2.nal" ]; then
    prefix="${DB_DIR}/MIDORI2_UNIQ_NUC_GB268_CO1_BLAST2"
  else
    shopt -s nullglob
    local nsq_files=( "${DB_DIR}"/*.nsq )
    local nal_files=( "${DB_DIR}"/*.nal )
    shopt -u nullglob
    if [ "${#nsq_files[@]}" -gt 0 ]; then
      prefix="${nsq_files[0]%.nsq}"
    elif [ "${#nal_files[@]}" -gt 0 ]; then
      prefix="${nal_files[0]%.nal}"
    fi
  fi
  echo "${prefix}"
}

blastdb_exists() {
  local p="$1"
  [ -n "${p}" ] && ( [ -f "${p}.nsq" ] || [ -f "${p}.nal" ] )
}

build_query_fasta() {
  local out_fa="$1"; shift
  local files=( "$@" )

  : > "${out_fa}"
  for f in "${files[@]}"; do
    if [ -s "${f}" ] && grep -q '^>' "${f}"; then
      cat "${f}" >> "${out_fa}"
      echo >> "${out_fa}"
    fi
  done

  grep -q '^>' "${out_fa}"
}

if [ "${MERGE_ONLY}" != "true" ]; then
  INPUT_FILES=()

  if [ -n "${SAMPLE_FILE}" ]; then
    [ -f "${SAMPLE_FILE}" ] || { echo "ERROR: SAMPLE_FILE not found: ${SAMPLE_FILE}"; exit 1; }
    INPUT_FILES+=("${SAMPLE_FILE}")
  else
    shopt -s nullglob
    INPUT_FILES+=( "${INPUT_DIR}"/*.fastq "${INPUT_DIR}"/*.fastq.gz )
    shopt -u nullglob
  fi

  [ "${#INPUT_FILES[@]}" -gt 0 ] || { echo "ERROR: No FASTQ/FASTQ.GZ found in ${INPUT_DIR}"; exit 1; }

  DB_PREFIX="$(pick_db_prefix)"
  if ! blastdb_exists "${DB_PREFIX}"; then
    echo "ERROR: No MIDORI BLAST database (*.nsq or *.nal) found in ${DB_DIR}"
    ls -lah "${DB_DIR}" || true
    exit 1
  fi
  echo "Using MIDORI BLAST DB prefix: ${DB_PREFIX}"

  for READS in "${INPUT_FILES[@]}"; do
    if [ -n "${SAMPLE_NAME_OVERRIDE}" ]; then
      SAMPLE="${SAMPLE_NAME_OVERRIDE}"
    else
      SAMPLE=$(basename "$READS")
      SAMPLE="${SAMPLE%.fastq.gz}"
      SAMPLE="${SAMPLE%.fastq}"
    fi

    echo ""
    echo "Processing sample: ${SAMPLE}"

    mkdir -p "${OUTPUT_DIR}/01_qc/${SAMPLE}"
    mkdir -p "${OUTPUT_DIR}/02_consensus/${SAMPLE}"

    # STEP 1: QC
    if [[ "$READS" == *.gz ]]; then
      READ_CMD=(pigz -dc "$READS")
    else
      READ_CMD=(cat "$READS")
    fi

    "${READ_CMD[@]}" | NanoFilt -q 15 -l 600 --maxlength 2000 \
      > "${OUTPUT_DIR}/01_qc/${SAMPLE}/filtered.fastq"

    LINES=$(wc -l < "${OUTPUT_DIR}/01_qc/${SAMPLE}/filtered.fastq" || echo 0)
    READS_N=$(( LINES / 4 ))
    if [ "${READS_N}" -le 0 ]; then
      echo "WARNING: 0 reads after QC for ${SAMPLE}; writing empty cluster+taxonomy and continue."
      echo -e "consensus_id\tread_count" > "${OUTPUT_DIR}/02_consensus/${SAMPLE}_clusters.tsv"
      echo -e "query\tdb_used\ttaxonomy\tconfidence\tnum_hits\tbest_pident\tbest_evalue\tbest_bitscore\tbest_accession\talignment_length\tkingdom\tkingdom_taxid\tphylum\tphylum_taxid\tclass\tclass_taxid\torder\torder_taxid\tfamily\tfamily_taxid\tgenus\tgenus_taxid\tspecies\tspecies_taxid" \
        > "${OUTPUT_DIR}/03_taxonomy/${SAMPLE}_taxonomy.tsv"
      continue
    fi

    # STEP 2: Primer trimming (optional)
    FINAL_READS="${OUTPUT_DIR}/01_qc/${SAMPLE}/filtered.fastq"
    if [ "${REMOVE_PRIMERS}" = "true" ] && [ -f "${PRIMER_FILE}" ]; then
      cutadapt -g file:"${PRIMER_FILE}" \
               -a file:"${PRIMER_FILE}" \
               -e 0.2 \
               --discard-untrimmed \
               -j "${THREADS}" \
               -o "${OUTPUT_DIR}/01_qc/${SAMPLE}/trimmed.fastq" \
               "${FINAL_READS}" \
               &> "${OUTPUT_DIR}/00_logs/cutadapt_${SAMPLE}.log"
      FINAL_READS="${OUTPUT_DIR}/01_qc/${SAMPLE}/trimmed.fastq"
    fi

    # STEP 3: NGSpeciesID
    NGSpeciesID --ont \
      --fastq "${FINAL_READS}" \
      --outfolder "${OUTPUT_DIR}/02_consensus/${SAMPLE}" \
      --consensus \
      --medaka \
      --t "${THREADS}" \
      &> "${OUTPUT_DIR}/00_logs/ngspeciesid_${SAMPLE}.log"

    # STEP 4: Cluster membership
    python /app/scripts/extract_cluster_membership.py \
      --input_dir "${OUTPUT_DIR}/02_consensus/${SAMPLE}" \
      --output "${OUTPUT_DIR}/02_consensus/${SAMPLE}_clusters.tsv"

    # STEP 5: Taxonomy
    CONS_DIR="${OUTPUT_DIR}/02_consensus/${SAMPLE}"
    CONS_CANDIDATES=()
    [ -f "${CONS_DIR}/consensus_references.fasta" ] && CONS_CANDIDATES+=("${CONS_DIR}/consensus_references.fasta")
    shopt -s nullglob
    CONS_CANDIDATES+=( "${CONS_DIR}"/consensus_reference*.fasta )
    shopt -u nullglob

    if [ "${#CONS_CANDIDATES[@]}" -eq 0 ]; then
      echo "WARNING: No consensus FASTA found in: ${CONS_DIR}"
      ls -lah "${CONS_DIR}" || true
      continue
    fi

    # If more than one consensus FASTA exists, build one query file
    if [ "${#CONS_CANDIDATES[@]}" -eq 1 ]; then
      CONSENSUS_FILE="${CONS_CANDIDATES[0]}"
      if ! grep -q '^>' "${CONSENSUS_FILE}"; then
        echo "WARNING: Consensus FASTA has no sequences for ${SAMPLE}: ${CONSENSUS_FILE}"
        continue
      fi
    else
      CONSENSUS_FILE="${CONS_DIR}/consensus_for_blast.fasta"
      if ! build_query_fasta "${CONSENSUS_FILE}" "${CONS_CANDIDATES[@]}"; then
        echo "WARNING: Built consensus_for_blast.fasta but it contains no sequences for ${SAMPLE}"
        continue
      fi
    fi

    # IMPORTANT: include staxids/sscinames so nt fallback can be resolved via taxdump
    BLAST_OUTFMT='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames stitle'

    BLAST_OUT="${OUTPUT_DIR}/03_taxonomy/${SAMPLE}_blast_hits.tsv"
    DB_USED_OUT="${OUTPUT_DIR}/03_taxonomy/${SAMPLE}_blast_db_used.txt"

    echo "Taxonomy assignment (MIDORI first)..."

    # 1) MIDORI
    blastn -query "${CONSENSUS_FILE}" \
      -db "${DB_PREFIX}" \
      -strand both \
      -outfmt "${BLAST_OUTFMT}" \
      -max_target_seqs "${MAX_TARGETS}" \
      -perc_identity "${MIN_IDENT}" \
      -evalue "${MAX_EVALUE}" \
      -num_threads "${THREADS}" \
      > "${BLAST_OUT}"

    DB_USED="MIDORI"

    # 2) fallback to nt if MIDORI has zero hits
    if [ ! -s "${BLAST_OUT}" ]; then
      if [ -n "${NT_PREFIX}" ] && blastdb_exists "${NT_PREFIX}"; then
        echo "WARNING: No MIDORI hits for ${SAMPLE}. Falling back to nt: ${NT_PREFIX}"
        blastn -query "${CONSENSUS_FILE}" \
          -db "${NT_PREFIX}" \
          -strand both \
          -outfmt "${BLAST_OUTFMT}" \
          -max_target_seqs "${MAX_TARGETS}" \
          -perc_identity "${MIN_IDENT}" \
          -evalue "${MAX_EVALUE}" \
          -num_threads "${THREADS}" \
          > "${BLAST_OUT}"
        DB_USED="NT"
      else
        echo "WARNING: No MIDORI hits and nt fallback not available."
      fi
    fi

    echo "${DB_USED}" > "${DB_USED_OUT}"

    # Assign taxonomy (best hit) with:
    # - bitscore filter
    # - taxdump-based rank mapping for nt
    # - skip generic first hit like "... sp."
    python /app/scripts/assign_lca.py \
      --input "${BLAST_OUT}" \
      --output "${OUTPUT_DIR}/03_taxonomy/${SAMPLE}_taxonomy.tsv" \
      --hits_output "${OUTPUT_DIR}/03_taxonomy/${SAMPLE}_all_hits.tsv" \
      --db_used "${DB_USED}" \
      --min_hits 1 \
      --min_identity "${MIN_IDENT}" \
      --max_evalue "${MAX_EVALUE}" \
      --min_bitscore "${MIN_BITSCORE}" \
      --taxdump_dir "${TAXDUMP_DIR}" \
      &> "${OUTPUT_DIR}/00_logs/assign_lca_${SAMPLE}.log"

    # Collapse only if ALL consensus resolve to the SAME species
    python /app/scripts/collapse_taxonomy_if_same_species.py \
      --taxonomy_tsv "${OUTPUT_DIR}/03_taxonomy/${SAMPLE}_taxonomy.tsv" \
      --clusters_tsv "${OUTPUT_DIR}/02_consensus/${SAMPLE}_clusters.tsv" \
      --output_tsv "${OUTPUT_DIR}/03_taxonomy/${SAMPLE}_taxonomy.tsv" \
      &> "${OUTPUT_DIR}/00_logs/collapse_taxonomy_${SAMPLE}.log"

    echo "Finished: ${SAMPLE}"
  done
else
  echo "MERGE_ONLY=true -> skipping per-sample steps."
fi

echo ""
if [ "${DO_FINAL_TABLE}" = "true" ]; then
  echo "Creating final taxonomy table..."
  python /app/scripts/create_final_taxonomy_table.py \
    "${OUTPUT_DIR}" \
    "${OUTPUT_DIR}/final_taxonomy_table.tsv"
else
  echo "Skipping final taxonomy table (DO_FINAL_TABLE=false)"
fi

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

echo "========================================="
echo "PIPELINE FINISHED"
echo "Final table: ${OUTPUT_DIR}/final_taxonomy_table.tsv"
echo "Runtime: $((DURATION/3600))h $((DURATION%3600/60))m $((DURATION%60))s"
echo "=========================================""
