# =============================================================================
# scn2a_haplotype_config.sh
# Configuration file for SCN2A haplotype analysis pipeline
# =============================================================================

# ---------------------------------------------------------------------------
# Tools
# ---------------------------------------------------------------------------

plink="/data0/software/plink/plink_1.9.0/plink"
bcftools="/data0/software/bcftools/bcftools-1.3/bcftools"
samtools="/data0/software/samtools/samtools-1.3/samtools"

# ---------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------
OUTDIR="results/pipeline_SCN2A_coding_unrel2504"
REGION_VCF="results/pipeline_SCN2A_coding_unrel2504.vcf.gz"

# Minimum het count for ASOG table inclusion
ASOG_MIN_HET=626

# Optional: path to file of sample IDs to subset to (one per line)
# Leave empty to use all samples
SAMPLE_SUBSET_FILE="data/1kG_30x_2504unrelated.txt"

# ---------------------------------------------------------------------------
# REFERENCES
# ---------------------------------------------------------------------------
# Remote VCF (NYGC 30x phased, chr2 only)
VCF_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr2.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

# Local paths (checked before downloading)
CHR2_VCF="data/1kGP_high_coverage_Illumina.chr2.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

# Reference genome
REF_FASTA="data/Homo_sapiens_assembly38.fasta"

# Sample Metadata
METADATA_FILE="data/igsr-1000genomes_sampleMap_30x_grch38-samples.tsv"


# ---------------------------------------------------------------------------
# Genomic region
# ---------------------------------------------------------------------------
REGION="chr2:165295823-165389822"
REGION_CHROM="chr2"
REGION_START=165295823
REGION_END=165389822

# ---------------------------------------------------------------------------
# PLINK parameters
# ---------------------------------------------------------------------------
PLINK_MIN_MAF=0.15
PLINK_MAX_KB=1000
PLINK_LD_WINDOW_KB=1000
PLINK_LD_WINDOW_SNP=99999

# ---------------------------------------------------------------------------
# Haplotype analysis thresholds
# ---------------------------------------------------------------------------
# Minimum block size (bp) for sequence-level analysis
MIN_BLOCK_SIZE_BP=5000

# Minimum het count for haplotype visualization in figures
VIZ_MIN_HET=100

# ---------------------------------------------------------------------------
# ASOG table parameters
# ---------------------------------------------------------------------------
# Window size (bp) centered on each variant for RefSeq/AltSeq extraction
WINDOW_SIZE=30

# GC content range for filtering (informational only; not used for ranking)
GC_MIN=0.40
GC_MAX=0.60

# Baseline ASO target SNPs used to compute UniqueSamplesAdded
# Format: "rsID:CHROM:POS" (one per line, space-separated list)
BASELINE_ASOS="rs1368238:chr2:165375124 rs72874313:chr2:165322227"

