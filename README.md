# SCN2A Haplospecific ASO Pipeline
Haplotype analysis pipeline supporting the identification of population-common SCN2A variants for haplospecific antisense oligonucleotide (ASO) therapy targeting. Implements the 1000 Genomes Project (1kGP) haplotype analysis described in:

Kim-McManus et al. From N-of-1 to N-of-Many: A haplospecific ASO approach for SCN2A-associated childhood neurological disease. American Journal of Human Genetics (in revision).

### Overview
The pipeline takes phased short-read genotype data from the 1kGP and identifies haplotype blocks across the SCN2A coding region, characterises common haplotypes by population, and generates ranked candidate ASO target variants for in silico screening.

Repository structure
```
project_dir/
├── data/
├── results/
├── scn2a_haplotype_pipeline.sh       # Main pipeline (Steps 1–7)
├── scn2a_haplotype_config.sh         # Configuration (edit before running)
├── scn2a_parse_haplotypes.R          # Step 6: haplotype parsing, Table S1, figures
├── scn2a_generate_asog_table.py      # Step 7: Table S2 and ASOG input table
└── README.md
```

## Dependencies
Tool | Version tested
---|---
bcftools | 1.3
PLINK | 1.9
tabix | any
samtools | any
python3 | 3.9+
R | 4.0+

R packages: `data.table`, `dplyr`, `tidyr`, `ggplot2`, `forcats`, `viridis`, `patchwork`, `gridExtra`

Python packages: none beyond the standard library

## Input data
The following files must be present before running. Set paths in scn2a_haplotype_config.sh.

File | Source
-----|------
`data/Homo_sapiens_assembly38.fasta` (`.fai`)	| GRCh38 reference genome
`data/igsr-1000genomes_sampleMap_30x_grch38-samples.tsv`	| [IGSR sample portal](https://www.internationalgenome.org/data-portal/sample)
`data/1kG_30x_2504unrelated.txt`	| From the [1000 Genomes 30x](https://www.internationalgenome.org/data-portal/data-collection/1000genomes_30x) resources: [Tab delimited file list for the 2504 panel]([url](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index))

The chr2 phased VCF and SCN2A region VCF are downloaded and extracted automatically by the pipeline if not already present

## Configuration
Edit scn2a_haplotype_config.sh before running. Key parameters:

```
# Tool paths
bcftools="/path/to/bcftools"
plink="/path/to/plink"
samtools="/path/to/samtools"

# Optional: restrict to unrelated samples only (reviewer reanalysis)
SAMPLE_SUBSET_FILE="data/1kG_30x_2504unrelated.txt"

# PLINK haplotype block parameters
PLINK_MIN_MAF=0.15          # Minor allele frequency threshold
PLINK_MAX_KB=1000           # Maximum block size

# ASOG candidate thresholds
ASOG_MIN_HET=800            # Minimum het carriers for ASOG table (~25% of cohort)
VIZ_MIN_HET=100             # Minimum het carriers for haplotype plot

# ASO sequence window
WINDOW_SIZE=30              # bp window centred on each variant

# Baseline ASOs for UniqueSamplesAdded calculation
BASELINE_ASOS="rs1368238:chr2:165375124 rs72874313:chr2:165322227"
```

## Usage
```
# Full cohort (all 3,202 samples)
bash scn2a_haplotype_pipeline.sh scn2a_haplotype_config.sh

# Unrelated samples only (reviewer reanalysis, n=2,504)
# Set SAMPLE_SUBSET_FILE in config, then:
bash scn2a_haplotype_pipeline.sh scn2a_haplotype_config.sh
Each run produces a timestamped log at ${OUTDIR}/run_YYYYMMDD_HHMMSS.log.
```

## Outputs
```
${OUTDIR}/
├── run_YYYYMMDD_HHMMSS.log
├── analysis_input.vcf.gz             # Optionally subsetted input VCF
├── plink/                            # PLINK intermediate files
├── hapsample/                        # bcftools hapsample output
├── haplotype_blocks/
│   ├── block_NN_XXXXXX_XXXXXX.tsv   # Per-block haplotype table (Table S1 input)
│   ├── block_NN_haplotype_plot.pdf   # Haplotype pattern figure
│   ├── block_NN_haplotype_plot.svg
│   ├── block_NN_ld_heatmap.pdf       # D' LD matrix figure
│   ├── block_NN_ld_heatmap.svg
│   ├── block_NN_combined.png         # Combined haplotype + LD figure (panels A/B)
│   └── all_blocks_haplotype_summary.tsv  # Table S1
├── snp_summary_table.tsv             # Table S2 (block SNPs only)
├── all_variants_summary.tsv          # All region variants (ASOG context)
└── asog_input_table.tsv              # Ranked ASOG candidate input table
```

### Table descriptions
**Table S1** (all_blocks_haplotype_summary.tsv): One row per haplotype per block. Columns include haplotype ID, block coordinates, SNP count, number of defining alt alleles, ASO target SNP flag, and Hom/Het/Total counts for the overall cohort and each of five 1kGP superpopulations (AFR, AMR, EAS, EUR, SAS).

**Table S2** (snp_summary_table.tsv): One row per informative SNP (MAF ≥15%, within a block ≥5kb). Includes total and superpopulation-stratified Het/HomAlt/HomRef counts, haplotype block assignment, number and identity of haplotypes carrying each allele, and the full list of heterozygous sample IDs.

**ASOG input table** (asog_input_table.tsv): Ranked candidates for in silico ASO screening. Ranked by UniqueSamplesAdded (het samples not already covered by baseline ASOs), with AlleleSpecificScore and HetSamples as tiebreakers. Includes 30bp reference and alternate sequences centred on each variant for ASOG input.

---

## Citation
If you use this pipeline, please cite:

> Kim-McManus et al. From N-of-1 to N-of-Many: A haplospecific ASO approach for SCN2A-associated childhood neurological disease. American Journal of Human Genetics (in revision).

## License
MIT
