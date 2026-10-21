#!/usr/bin/env bash
# =============================================================================
# scn2a_haplotype_pipeline.sh
#
# SCN2A haplotype analysis pipeline for ASO target identification.
# Reproduces Tables S1, S2, and ASOG input table from the manuscript,
# and supports reanalysis in subsets (e.g. unrelated 1kGP samples).
#
# Usage:
#   bash scn2a_haplotype_pipeline.sh [config_file]
#
# If no config file is provided, defaults to scn2a_haplotype_config.sh
# in the same directory as this script.
#
# Key outputs (in OUTDIR):
#   haplotype_blocks/          Per-block haplotype TSVs (Table S1)
#   snp_summary_table.tsv      SNP-level summary (Table S2)
#   asog_input_table.tsv       ASOG input table with ranking
#   run.log                    Automated run log with cohort stats
#
# Requirements: bcftools, plink (v1.9), tabix, bgzip, curl, python3, samtools/faidx
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Load config
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${1:-${SCRIPT_DIR}/scn2a_haplotype_config.sh}"

if [[ ! -f "${CONFIG}" ]]; then
    echo "ERROR: Config file not found: ${CONFIG}"
    exit 1
fi
source "${CONFIG}"

# ---------------------------------------------------------------------------
# Setup output directory and log
# ---------------------------------------------------------------------------
mkdir -p "${OUTDIR}" "${OUTDIR}/haplotype_blocks" "${OUTDIR}/plink" "${OUTDIR}/hapsample"

RUN_ID="$(date '+%Y%m%d_%H%M%S')"
LOG="${OUTDIR}/run_${RUN_ID}.log"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG}"; }

log "=== SCN2A Haplotype Pipeline ==="
log "Config: ${CONFIG}"
log "Run ID: ${RUN_ID}"

# ---------------------------------------------------------------------------
# Validate required tools
# ---------------------------------------------------------------------------
log "--- Checking dependencies ---"
for tool in "${bcftools}" "${plink}" "${samtools}" tabix bgzip curl python3; do
    if ! command -v "${tool}" &>/dev/null; then
        log "ERROR: Required tool not found: ${tool}"
        exit 1
    fi
    log "  ${tool}: $(command -v ${tool})"
done

# Check faidx-indexed reference
if [[ ! -f "${REF_FASTA}" ]]; then
    log "ERROR: Reference FASTA not found: ${REF_FASTA}"
    exit 1
fi
if [[ ! -f "${REF_FASTA}.fai" ]]; then
    log "Indexing reference FASTA..."
    "${samtools}" faidx "${REF_FASTA}"
fi

# ---------------------------------------------------------------------------
# Step 1: Download chr2 VCF if not present
# ---------------------------------------------------------------------------
log "--- Step 1: Data acquisition ---"

# The panel is a single multi-gigabyte file and this transfer has been seen to die partway
# with a connection reset. --retry on its own restarts from byte zero, which turns one failed
# transfer into several and can spend a long time leaving nothing behind; --continue-at is
# the flag that resumes, and the retry flags ride on top of it.
fetch_resumable() {
    local url="$1" dest="$2"
    curl --location --fail \
         --continue-at - \
         --retry 5 --retry-delay 10 --retry-all-errors \
         --output "${dest}" "${url}"
}

# The panel and its index are checked separately, so a missing index is fetched on its own
# instead of forcing the whole panel to be downloaded again. A transfer that stopped short
# still leaves a readable gzip header, so presence alone does not mean complete: bgzip -t
# reads every block's checksum, and a panel that fails it is resumed rather than trusted.
# That test reads the whole file, so it costs a pass over the panel on every run.
if [[ -f "${CHR2_VCF}" ]] && bgzip -t "${CHR2_VCF}" 2>/dev/null; then
    log "Chr2 VCF already present and complete, skipping download: ${CHR2_VCF}"
else
    if [[ -f "${CHR2_VCF}" ]]; then
        log "Chr2 VCF is present but incomplete; resuming the transfer"
    else
        log "Downloading chr2 VCF from 1kGP FTP..."
    fi
    # A file that is corrupt rather than short is already at its full length, so the server
    # has nothing left to send and refuses the resume. Say so, because the fix is different.
    if ! fetch_resumable "${VCF_URL}" "${CHR2_VCF}"; then
        log "ERROR: the panel transfer failed. If ${CHR2_VCF} is present and already the"
        log "       full size, the server cannot resume it: delete it and run again."
        exit 1
    fi

    log "Verifying the downloaded panel is complete..."
    if ! bgzip -t "${CHR2_VCF}"; then
        log "ERROR: ${CHR2_VCF} is still truncated or corrupt after the transfer."
        log "       Re-run to resume; if it keeps failing, delete the file and start over."
        exit 1
    fi
    log "Download complete and verified: ${CHR2_VCF}"
fi

if [[ -f "${CHR2_VCF}.tbi" ]]; then
    log "Chr2 VCF index already exists, skipping download: ${CHR2_VCF}.tbi"
else
    log "Downloading chr2 VCF index..."
    if ! fetch_resumable "${VCF_URL}.tbi" "${CHR2_VCF}.tbi"; then
        log "ERROR: the index transfer failed. If ${CHR2_VCF}.tbi is present and already the"
        log "       full size, the server cannot resume it: delete it and run again."
        exit 1
    fi
    log "Index downloaded: ${CHR2_VCF}.tbi"
fi

# ---------------------------------------------------------------------------
# Step 2: Subset to SCN2A region if not present
# ---------------------------------------------------------------------------
if [[ -f "${REGION_VCF}" && -f "${REGION_VCF}.tbi" ]]; then
    log "Region VCF already exists, skipping subsetting: ${REGION_VCF}"
else
    log "Subsetting to SCN2A region: ${REGION}"
    ${bcftools} view \
        --regions "${REGION}" \
        --output-type z \
        --output-file "${REGION_VCF}" \
        "${CHR2_VCF}"
    tabix -p vcf "${REGION_VCF}"
    log "Region VCF created: ${REGION_VCF}"
fi

# ---------------------------------------------------------------------------
# Step 3: Optional sample subsetting
# ---------------------------------------------------------------------------
log "--- Step 3: Sample subsetting ---"

ANALYSIS_VCF="${OUTDIR}/analysis_input.vcf.gz"

if [[ -n "${SAMPLE_SUBSET_FILE}" && -f "${SAMPLE_SUBSET_FILE}" ]]; then
    log "Subsetting to samples in: ${SAMPLE_SUBSET_FILE}"
    N_SUBSET=$(wc -l < "${SAMPLE_SUBSET_FILE}")
    log "  Requested samples: ${N_SUBSET}"

    ${bcftools} view \
        --samples-file "${SAMPLE_SUBSET_FILE}" \
        --output-type z \
        --output-file "${ANALYSIS_VCF}" \
        "${REGION_VCF}"
    tabix -p vcf "${ANALYSIS_VCF}"

    N_SAMPLES=$(${bcftools} query -l "${ANALYSIS_VCF}" | wc -l)
    log "  Samples in output VCF: ${N_SAMPLES}"
else
    log "No sample subset file specified; using all samples"
    ln -sf "$(realpath "${REGION_VCF}")" "${ANALYSIS_VCF}"
    ln -sf "$(realpath "${REGION_VCF}.tbi")" "${ANALYSIS_VCF}.tbi"
    N_SAMPLES=$(${bcftools} query -l "${ANALYSIS_VCF}" | wc -l)
    log "  Total samples: ${N_SAMPLES}"
fi

COHORT_SIZE=${N_SAMPLES}
log "  Cohort size for this run: ${COHORT_SIZE}"

# ---------------------------------------------------------------------------
# Step 4: PLINK haplotype block detection
# ---------------------------------------------------------------------------
log "--- Step 4: PLINK haplotype block detection ---"
log "  Parameters: --blocks-min-maf ${PLINK_MIN_MAF} --blocks-max-kb ${PLINK_MAX_KB}"

PLINK_PREFIX="${OUTDIR}/plink/scn2a"

# Convert VCF to PLINK format
${bcftools} view "${ANALYSIS_VCF}" | \
    ${plink} \
        --vcf /dev/stdin \
        --vcf-half-call m \
        --make-bed \
        --out "${PLINK_PREFIX}" \
        --allow-extra-chr \
        --chr-set 95 \
        2>>"${LOG}"

# Detect haplotype blocks
${plink} \
    --bfile "${PLINK_PREFIX}" \
    --blocks no-pheno-req \
    --blocks-max-kb "${PLINK_MAX_KB}" \
    --blocks-min-maf "${PLINK_MIN_MAF}" \
    --out "${PLINK_PREFIX}" \
    --allow-extra-chr \
    2>>"${LOG}"

# Compute LD (r2 and D')
${plink} \
    --bfile "${PLINK_PREFIX}" \
    --r2 dprime \
    --ld-window-kb "${PLINK_LD_WINDOW_KB}" \
    --ld-window "${PLINK_LD_WINDOW_SNP}" \
    --ld-window-r2 0 \
    --out "${PLINK_PREFIX}" \
    --allow-extra-chr \
    2>>"${LOG}"

N_BLOCKS=$(wc -l < "${PLINK_PREFIX}.blocks.det" 2>/dev/null || echo 0)
log "  Haplotype blocks detected: $((N_BLOCKS - 1))"  # subtract header

# ---------------------------------------------------------------------------
# Step 5: Export hapsample format
# ---------------------------------------------------------------------------
log "--- Step 5: Export hapsample format ---"

HAPSAMPLE_PREFIX="${OUTDIR}/hapsample/scn2a"
${bcftools} convert \
    --hapsample "${HAPSAMPLE_PREFIX}" \
    "${ANALYSIS_VCF}" \
    2>>"${LOG}"

log "  Hapsample files: ${HAPSAMPLE_PREFIX}.hap.gz, ${HAPSAMPLE_PREFIX}.sample"

# ---------------------------------------------------------------------------
# Step 6: Parse haplotypes and generate per-block TSVs + figures
# ---------------------------------------------------------------------------
log "--- Step 6: Parse haplotypes (R) ---"

Rscript "${SCRIPT_DIR}/scn2a_parse_haplotypes.R" \
    "${PLINK_PREFIX}" \
    "${HAPSAMPLE_PREFIX}" \
    "${OUTDIR}/haplotype_blocks" \
    "${MIN_BLOCK_SIZE_BP}" \
    "${VIZ_MIN_HET}" \
    "${METADATA_FILE}" \
    "${LOG}" \
    "${BASELINE_ASOS}" \
    2>>"${LOG}"

log "  Per-block TSVs written to: ${OUTDIR}/haplotype_blocks/"

# ---------------------------------------------------------------------------
# Step 7: Generate ASOG input table
# ---------------------------------------------------------------------------
log "--- Step 7: Generate ASOG input table (Python) ---"

python3 "${SCRIPT_DIR}/scn2a_generate_asog_table.py" \
    --vcf "${ANALYSIS_VCF}" \
    --blocks_dir "${OUTDIR}/haplotype_blocks" \
    --ref_fasta "${REF_FASTA}" \
    --samtools "${samtools}" \
    --region_chrom "${REGION_CHROM}" \
    --region_start "${REGION_START}" \
    --region_end "${REGION_END}" \
    --window_size "${WINDOW_SIZE}" \
    --asog_min_het "${ASOG_MIN_HET}" \
    --cohort_size "${COHORT_SIZE}" \
    --baseline_asos "${BASELINE_ASOS}" \
    --metadata_file "${METADATA_FILE}" \
    --gc_min "${GC_MIN}" \
    --gc_max "${GC_MAX}" \
    --outdir "${OUTDIR}" \
    --log "${LOG}" \
    2>>"${LOG}"

log "  ASOG input table: ${OUTDIR}/asog_input_table.tsv"
log "  SNP summary table: ${OUTDIR}/snp_summary_table.tsv"

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
log ""
log "=== Run Complete ==="
log "  Run ID:          ${RUN_ID}"
log "  Cohort size:     ${COHORT_SIZE}"
log "  Blocks detected: $((N_BLOCKS - 1))"
log "  Output dir:      ${OUTDIR}"
log "  Full log:        ${LOG}"
