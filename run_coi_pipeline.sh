#!/bin/bash
# COI taxonomy-only pipeline: QC → optional primer trim → NGSpeciesID → BLAST → assign_lca → collapse_if_same_species → final table
#
# HPC-safe design for Slurm arrays:
# - In array mode (one sample per task): writes ONLY per-sample outputs + per-sample qc summary line
# - In merge mode (MERGE_ONLY=true): builds combined qc_summary.tsv, runs MultiQC once, and creates final taxonomy table
#
# Controls:
#   MERGE_ONLY=false (default): process samples (use SAMPLE_FILE or SAMPLE_LIST or auto-discover)
#   MERGE_ONLY=true:  do only merge steps (MultiQC + qc_summary + final taxonomy table)

set -euo pipefail

START_TIME=$(date +%s)

# -------------------------
# Configuration
# -------------------------
INPUT_DIR="${INPUT_DIR:-/data/input}"
OUTPUT_DIR="${OUTPUT_DIR:-/data/output}"
DB_DIR="${DB_DIR:-/data/databases}"
THREADS="${THREADS:-4}"

# For Slurm array usage:
# - Provide SAMPLE_FILE directly, OR
# - Provide SAMPLE_LIST (file containing one fastq(.gz)/fq(.gz) path per line) + SAMPLE_INDEX (0-based)
SAMPLE_FILE="${SAMPLE_FILE:-}"                   # optional: process only this fastq(.gz)/fq(.gz)
SAMPLE_LIST="${SAMPLE_LIST:-}"                   # optional: list of input files (one per line)
SAMPLE_INDEX="${SAMPLE_INDEX:-}"                 # optional: 0-based index into SAMPLE_LIST
SAMPLE_NAME_OVERRIDE="${SAMPLE_NAME_OVERRIDE:-}" # optional fixed sample name

MERGE_ONLY="${MERGE_ONLY:-false}"                # skip per-sample; only merge
DO_FINAL_TABLE="${DO_FINAL_TABLE:-true}"         # used in merge step

REMOVE_PRIMERS="${REMOVE_PRIMERS:-false}"
PRIMER_FILE="${INPUT_DIR}/primers.fasta"

MIN_IDENT="${MIN_IDENT:-80.0}"
MAX_EVALUE="${MAX_EVALUE:-1e-10}"
MAX_TARGETS="${MAX_TARGETS:-10}"
MIN_BITSCORE="${MIN_BITSCORE:-400}"

NT_PREFIX="${NT_PREFIX:-}"
TAXDUMP_DIR="${TAXDUMP_DIR:-}"

# -------------------------
# Logging / dirs
# -------------------------
echo "========================================="
echo "COI TAXONOMY PIPELINE"
echo "Input dir:     ${INPUT_DIR}"
echo "Output dir:    ${OUTPUT_DIR}"
echo "DB dir:        ${DB_DIR}"
echo "Threads:       ${THREADS}"
echo "Merge only:    ${MERGE_ONLY}"
echo "Final table:   ${DO_FINAL_TABLE}"
echo "Primer trim:   ${REMOVE_PRIMERS}"
echo "BLAST filters: pident>=${MIN_IDENT}, evalue<=${MAX_EVALUE}, max_targets=${MAX_TARGETS}, bitscore>=${MIN_BITSCORE}"
echo "nt fallback:   ${NT_PREFIX:-<disabled>}"
echo "taxdump dir:   ${TAXDUMP_DIR:-<disabled>}"
echo "========================================="

mkdir -p "${OUTPUT_DIR}"/{00_logs,01_qc,02_consensus,03_taxonomy}

# -------------------------
# Helpers
# -------------------------
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

# Return a single input file in $1, based on SAMPLE_FILE / SAMPLE_LIST+SAMPLE_INDEX / auto-discovery
pick_one_input() {
  local out_var="$1"

  if [ -n "${SAMPLE_FILE}" ]; then
    [ -f "${SAMPLE_FILE}" ] || { echo "ERROR: SAMPLE_FILE not found: ${SAMPLE_FILE}"; exit 1; }
    printf -v "${out_var}" "%s" "${SAMPLE_FILE}"
    return 0
  fi

  if [ -n "${SAMPLE_LIST}" ] && [ -n "${SAMPLE_INDEX}" ]; then
    [ -f "${SAMPLE_LIST}" ] || { echo "ERROR: SAMPLE_LIST not found: ${SAMPLE_LIST}"; exit 1; }
    local line
    line="$(sed -n "$((SAMPLE_INDEX+1))p" "${SAMPLE_LIST}" | tr -d '\r')"
    [ -n "${line}" ] || { echo "ERROR: No line for SAMPLE_INDEX=${SAMPLE_INDEX} in ${SAMPLE_LIST}"; exit 1; }
    [ -f "${line}" ] || { echo "ERROR: File from SAMPLE_LIST not found: ${line}"; exit 1; }
    printf -v "${out_var}" "%s" "${line}"
    return 0
  fi

  # Fallback: auto-discovery only if not in array mode
  shopt -s nullglob
  local files=( "${INPUT_DIR}"/*.fastq "${INPUT_DIR}"/*.fastq.gz "${INPUT_DIR}"/*.fq "${INPUT_DIR}"/*.fq.gz )
  shopt -u nullglob
  [ "${#files[@]}" -gt 0 ] || { echo "ERROR: No FASTQ/FASTQ.GZ/FQ/FQ.GZ found in ${INPUT_DIR}"; exit 1; }
  # If multiple files are found and you didn't specify which one, this script will process ALL (non-array local mode)
  printf -v "${out_var}" "%s" "__ALL__"
  return 0
}

# Derive sample name from reads filename unless overridden
derive_sample_name() {
  local reads="$1"
  if [ -n "${SAMPLE_NAME_OVERRIDE}" ]; then
    echo "${SAMPLE_NAME_OVERRIDE}"
    return
  fi
  local s
  s="$(basename "${reads}")"
  s="${s%.fastq.gz}"
  s="${s%.fastq}"
  s="${s%.fq.gz}"
  s="${s%.fq}"
  echo "${s}"
}

# Write per-sample QC summary line safely (one file per sample, no races)
write_qc_line() {
  local sample="$1"
  local reads_n="$2"
  local out="${OUTPUT_DIR}/01_qc/${sample}/qc_reads.tsv"
  echo -e "${sample}\t${reads_n}" > "${out}"
}

# Merge per-sample qc_reads.tsv into qc_summary.tsv (merge job)
merge_qc_summary() {
  local qc_summary="${OUTPUT_DIR}/01_qc/qc_summary.tsv"
  echo -e "sample\treads_after_nanofilt" > "${qc_summary}"
  find "${OUTPUT_DIR}/01_qc" -type f -name "qc_reads.tsv" -print0 \
    | xargs -0 cat \
    | sort -k1,1 \
    >> "${qc_summary}" || true
}

run_multiqc() {
  echo "Running MultiQC over QC outputs (all samples)..."

  local outdir="${OUTPUT_DIR}/01_qc/multiqc"
  rm -rf "${outdir}"
  mkdir -p "${outdir}"

  # IMPORTANT:
  # 1) write into a clean directory
  # 2) explicitly set report filename
  # 3) ignore the multiqc output directory itself (prevents recursion)
  # 4) --force overwrites if something exists
  multiqc "${OUTPUT_DIR}/01_qc" \
    --outdir "${outdir}" \
    --filename "multiqc_report.html" \
    --force \
    --ignore "multiqc" \
    &> "${OUTPUT_DIR}/00_logs/multiqc.log"
}

# -------------------------
# Main
# -------------------------
if [ "${MERGE_ONLY}" != "true" ]; then
  # Per-sample processing (array-safe)
  DB_PREFIX="$(pick_db_prefix)"
  if ! blastdb_exists "${DB_PREFIX}"; then
    echo "ERROR: No MIDORI BLAST database (*.nsq or *.nal) found in ${DB_DIR}"
    ls -lah "${DB_DIR}" || true
    exit 1
  fi
  echo "Using MIDORI BLAST DB prefix: ${DB_PREFIX}"

  ONE_INPUT=""
  pick_one_input ONE_INPUT

  INPUT_FILES=()
  if [ "${ONE_INPUT}" = "__ALL__" ]; then
    shopt -s nullglob
    INPUT_FILES+=( "${INPUT_DIR}"/*.fastq "${INPUT_DIR}"/*.fastq.gz "${INPUT_DIR}"/*.fq "${INPUT_DIR}"/*.fq.gz )
    shopt -u nullglob
  else
    INPUT_FILES+=( "${ONE_INPUT}" )
  fi

  for READS in "${INPUT_FILES[@]}"; do
    SAMPLE="$(derive_sample_name "${READS}")"

    echo ""
    echo "Processing sample: ${SAMPLE}"

    mkdir -p "${OUTPUT_DIR}/01_qc/${SAMPLE}"
    mkdir -p "${OUTPUT_DIR}/02_consensus/${SAMPLE}"

    # STEP 1: QC (NanoFilt)
    if [[ "${READS}" == *.gz ]]; then
      READ_CMD=(pigz -dc "${READS}")
    else
      READ_CMD=(cat "${READS}")
    fi

    FILTERED_FASTQ="${OUTPUT_DIR}/01_qc/${SAMPLE}/filtered.fastq"

    "${READ_CMD[@]}" | NanoFilt -q 15 -l 600 --maxlength 2000 \
      > "${FILTERED_FASTQ}"

    LINES=$(wc -l < "${FILTERED_FASTQ}" || echo 0)
    READS_N=$(( LINES / 4 ))
    write_qc_line "${SAMPLE}" "${READS_N}"

    if [ "${READS_N}" -le 0 ]; then
      echo "WARNING: 0 reads after QC for ${SAMPLE}; writing empty cluster+taxonomy and continue."
      echo -e "consensus_id\tread_count" > "${OUTPUT_DIR}/02_consensus/${SAMPLE}_clusters.tsv"
      echo -e "query\tdb_used\ttaxonomy\tconfidence\tnum_hits\tbest_pident\tbest_evalue\tbest_bitscore\tbest_accession\talignment_length\tkingdom\tkingdom_taxid\tphylum\tphylum_taxid\tclass\tclass_taxid\torder\torder_taxid\tfamily\tfamily_taxid\tgenus\tgenus_taxid\tspecies\tspecies_taxid" \
        > "${OUTPUT_DIR}/03_taxonomy/${SAMPLE}_taxonomy.tsv"
      continue
    fi

    # NanoPlot per sample
    NanoPlot --fastq "${FILTERED_FASTQ}" \
      -o "${OUTPUT_DIR}/01_qc/${SAMPLE}/nanoplot" \
      --threads "${THREADS}" \
      --prefix "${SAMPLE}_"  \
      --title "${SAMPLE}" \
      &> "${OUTPUT_DIR}/00_logs/nanoplot_${SAMPLE}.log"

    # STEP 2: Primer trimming (optional)
    FINAL_READS="${FILTERED_FASTQ}"
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

    BLAST_OUTFMT='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames stitle'
    BLAST_OUT="${OUTPUT_DIR}/03_taxonomy/${SAMPLE}_blast_hits.tsv"
    DB_USED_OUT="${OUTPUT_DIR}/03_taxonomy/${SAMPLE}_blast_db_used.txt"

    echo "Taxonomy assignment (MIDORI first)..."

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

    python /app/scripts/collapse_taxonomy_if_same_species.py \
      --taxonomy_tsv "${OUTPUT_DIR}/03_taxonomy/${SAMPLE}_taxonomy.tsv" \
      --clusters_tsv "${OUTPUT_DIR}/02_consensus/${SAMPLE}_clusters.tsv" \
      --output_tsv "${OUTPUT_DIR}/03_taxonomy/${SAMPLE}_taxonomy.tsv" \
      &> "${OUTPUT_DIR}/00_logs/collapse_taxonomy_${SAMPLE}.log"

    echo "Finished: ${SAMPLE}"
  done

  echo ""
  echo "Per-sample run complete. (MultiQC and qc_summary.tsv are generated in MERGE_ONLY=true mode.)"

else
  # Merge-only steps (single Slurm job after array finishes)
  echo "MERGE_ONLY=true -> generating combined QC summary, MultiQC report, and final table (optional)."

  merge_qc_summary
  run_multiqc

  echo ""
  if [ "${DO_FINAL_TABLE}" = "true" ]; then
    echo "Creating final taxonomy table..."
    python /app/scripts/create_final_taxonomy_table.py \
      "${OUTPUT_DIR}" \
      "${OUTPUT_DIR}/final_taxonomy_table.tsv"
  else
    echo "Skipping final taxonomy table (DO_FINAL_TABLE=false)"
  fi
fi

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

echo "========================================="
echo "PIPELINE FINISHED"
echo "QC summary:   ${OUTPUT_DIR}/01_qc/qc_summary.tsv"
echo "MultiQC:      ${OUTPUT_DIR}/01_qc/multiqc/multiqc_report.html"
echo "Final table:  ${OUTPUT_DIR}/final_taxonomy_table.tsv"
echo "Runtime: $((DURATION/3600))h $((DURATION%3600/60))m $((DURATION%60))s"
echo "========================================="
