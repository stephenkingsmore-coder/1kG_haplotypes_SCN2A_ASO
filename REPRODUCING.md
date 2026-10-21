# Setting up and reproducing the pipeline

A newcomer's guide to running `scn2a_haplotype_pipeline.sh` on a machine that is not
the authors' server. It covers the tools each step needs, the files that have to be
staged by hand, the folders the pipeline expects, a small-subset run that finishes in
minutes, and the handful of things that currently stop a clean run elsewhere.

This guide was written by reading the seven steps rather than by completing a full
run — no 1kGP panel or GRCh38 reference was available. Statements marked **verified**
were checked by executing the code; everything else is read from the scripts, so
treat the first real run as the test and correct this file from it.

---

## 1. Tools

| Tool | Used by | Notes |
|---|---|---|
| `bcftools` | Steps 2, 3, 4, 5 | `view --regions`, `view --samples-file`, `query -l`, `convert --hapsample` |
| `tabix` | Steps 2, 3 | indexes each VCF it writes; ships with htslib |
| `plink` | Step 4 | **must be PLINK 1.9** — `--blocks` does not exist in PLINK 2 |
| `Rscript` | Step 6 | not checked by the pipeline's dependency loop; a missing R fails at Step 6 |
| `python3` | Step 7 | standard library only, 3.7+ in practice (`subprocess.run(capture_output=)`) |
| `samtools` | reference indexing, Step 7 | `faidx` only |
| `curl` | Step 1 | only if you let the pipeline download the panel — see §6.2 |
| `bgzip` | Step 1 | checks the downloaded panel is complete; ships with htslib |

R packages: `data.table`, `dplyr`, `tidyr`, `ggplot2`, `forcats`, `viridis`,
`patchwork`, `gridExtra`, and two the README omits — **`rlang`** (used as
`rlang::sym`) and **`svglite`** (every block writes an `.svg`, and `ggsave()` cannot
write SVG without it). `ragg` is worth adding for the combined `.png` figures.

`environment.yml` in this repository pins a conda environment with all of the above:

```bash
conda env create -f environment.yml     # or mamba env create -f environment.yml
conda activate scn2a-aso
```

The pipeline reads tool paths from the config, not from `PATH`, so after creating the
environment point the three variables at it:

```bash
# in scn2a_haplotype_config.sh
plink="$CONDA_PREFIX/bin/plink"
bcftools="$CONDA_PREFIX/bin/bcftools"
samtools="$CONDA_PREFIX/bin/samtools"
```

`tabix`, `bgzip`, `curl`, `python3` and `Rscript` are called bare and must be on `PATH`. Every use of
samtools now goes through the config variable: the dependency check, the reference
indexing, and Step 7, which receives it as `--samtools`. Inside the conda environment
all of this holds.

---

## 2. Folders

The pipeline creates everything under `OUTDIR` itself. You create `data/`, and the
parent of `OUTDIR`, before the first run:

```bash
mkdir -p data results
```

```
project_dir/
├── data/                  # you stage these; nothing here is in git
│   ├── Homo_sapiens_assembly38.fasta            (+ .fai)
│   ├── 1kGP_high_coverage_Illumina.chr2...vcf.gz (+ .tbi)
│   ├── igsr-1000genomes_sampleMap_30x_grch38-samples.tsv
│   └── 1kG_30x_2504unrelated.txt
└── results/               # REGION_VCF lands here, OUTDIR is a subdirectory
    ├── pipeline_SCN2A_coding_unrel2504.vcf.gz   (+ .tbi)
    └── pipeline_SCN2A_coding_unrel2504/         = OUTDIR
        ├── run_YYYYMMDD_HHMMSS.log
        ├── analysis_input.vcf.gz                (+ .tbi; a symlink if no subsetting)
        ├── plink/        scn2a.{bed,bim,fam,blocks,blocks.det,ld}
        ├── hapsample/    scn2a.{hap.gz,sample}
        ├── haplotype_blocks/
        │   ├── block_NN_<bp1>_<bp2>.tsv
        │   ├── block_NN_haplotype_plot.{pdf,svg}
        │   ├── block_NN_ld_heatmap.{pdf,svg}
        │   ├── block_NN_combined.png
        │   └── all_blocks_haplotype_summary.tsv   = Table S1
        ├── snp_summary_table.tsv                  = Table S2
        ├── all_variants_summary.tsv
        └── asog_input_table.tsv
```

Both `REGION_VCF` and `OUTDIR` are reused if they already exist: Step 2 is skipped
when `REGION_VCF` and its `.tbi` are both present. Change `OUTDIR` for each distinct
run, or you will mix outputs from two cohorts in one folder.

---

## 3. Input files

Four files must be in place before the first run. None is in git, and only the panel
VCF is downloadable by the pipeline.

| File | Where it comes from | Approx. size |
|---|---|---|
| `Homo_sapiens_assembly38.fasta` | GRCh38 reference (Broad bundle); must be the `chr`-prefixed GRCh38 the panel is called against | ~3 GB |
| `1kGP_high_coverage_Illumina.chr2.filtered.SNV_INDEL_SV_phased_panel.vcf.gz` + `.tbi` | [1kGP 30x phased panel](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/) | 2.5 GB |
| `igsr-1000genomes_sampleMap_30x_grch38-samples.tsv` | [IGSR data portal](https://www.internationalgenome.org/data-portal/sample) | small |
| `1kG_30x_2504unrelated.txt` | the 2,504-sample list, one ID per line | small |

**The panel is one 2.5 GB file** (`Content-Length: 2506332497`, read by range request on
2026-10-20), and a transfer that size can drop partway: one attempt from this project's cloud
container died with `curl: (56) Recv failure: Connection reset by peer`. Step 1 therefore
fetches it with `curl --continue-at -`, which **resumes** rather than restarts — retrying
without resuming reads like the fix and is not, since each attempt would start again from byte
zero. It then runs `bgzip -t` over the finished file, because a transfer that stopped short
still leaves a readable gzip header and would otherwise surface as a confusing failure in
Step 2. A panel that is present but fails that test is resumed on the next run rather than
trusted, which costs one pass over the file per run; one that is corrupt rather than short is
already the full length, so there is nothing for the server to resume and the step tells you to
delete it and start again.

Know what that was tested against. **The resume is measured against the 1kGP route itself**, on
2026-10-21: a transfer of the panel killed mid-stream resumed with `Range: bytes=5726208-` and
`206` in a single request, and matched an independently fetched prefix byte for byte; the index,
truncated and resumed, reproduced a clean download's sha256 exactly. Three limits come with that.
No transfer of the 2.5 GB panel was carried through to the end, the largest being about 62 MB, so
`bgzip -t` has never run over a real panel; the machine that measured it had no `bgzip` at all;
and the no-validator path described below is unexercised, that host serving both headers. Two
further limits are worth knowing before you lean on any of it: `bgzip -t` validates each block's
checksum and not that the blocks are the right blocks in the right order, so if the provider
publishes a checksum, verify against that; and a resume of a file that has been re-released since
you were interrupted would splice the new bytes onto your old prefix, so Step 1 records the file's
`ETag` and `Last-Modified` beside the partial download and re-reads them before resuming,
starting over instead if either has changed. Do not reach for `If-Range` here, which is the
obvious idiom and the wrong one: a conforming server answers a stale `If-Range` with the whole
body, and curl, still resuming, would write that body at the resume offset. Where a host serves
neither header the step says so in the log and resumes anyway, and there `bgzip -t` is the only
remaining guard.

The metadata TSV must keep the IGSR portal's own header: both Step 6 and Step 7 look
up the columns **`Sample name`** and **`Superpopulation code`** by name, and Step 7
raises a `ValueError` naming the header it found if either is missing. The 2,504 list
is a bare list of sample IDs, one per line — not the IGSR sequence index, which has a
header and many columns. Extract it with `cut`/`awk` if you downloaded the index.

If the reference has no `.fai`, the pipeline creates one with `samtools faidx` before
Step 1; indexing a 3 GB FASTA takes a minute or two. Step 7 needs that index for more
than sequence fetching — it reads contig lengths from the `.fai` to tell a variant
near a contig end from a broken reference — and it checks for both samtools and the
index upfront rather than part way through scoring.

---

## 4. Choosing which run you are reproducing

The committed config is the **reviewer reanalysis** (2,504 unrelated samples), not the
primary analysis the README documents. Three values differ:

| | Full cohort (README §Usage) | Unrelated, as committed |
|---|---|---|
| `SAMPLE_SUBSET_FILE` | `""` | `data/1kG_30x_2504unrelated.txt` |
| `ASOG_MIN_HET` | `800` | `626` |
| `OUTDIR` / `REGION_VCF` | your own names | `results/pipeline_SCN2A_coding_unrel2504*` |

Both thresholds are ~25% of their cohort (800/3202, 626/2504), so scale
`ASOG_MIN_HET` with the cohort rather than carrying a number across runs. Nothing
else needs changing to switch between the two.

Everything else in the config is shared: the region
`chr2:165,295,823-165,389,822`, `PLINK_MIN_MAF=0.15`, `PLINK_MAX_KB=1000`,
`MIN_BLOCK_SIZE_BP=5000`, `VIZ_MIN_HET=100`, `WINDOW_SIZE=30`, and the two baseline
ASO SNPs.

> **On the MANE Select region:** this repository stays on the manuscript's
> coding-domain region; the MANE Select transcript run is being done in the
> `N-of-Many` repository, so look there before widening anything here. If you do
> widen it, the region is the only lever — `REGION`, `REGION_CHROM`, `REGION_START`
> and `REGION_END` must agree with each other, since Steps 2 and 7 read them
> separately — and it will raise the SNP and block counts while the baseline ASO
> positions stay as they are.

---

## 5. Running it end to end

```bash
conda activate scn2a-aso
bash scn2a_haplotype_pipeline.sh scn2a_haplotype_config.sh
```

The config path is optional and defaults to `scn2a_haplotype_config.sh` beside the
script. The driver is `set -euo pipefail`, so it stops at the first failing step;
each run appends to its own `${OUTDIR}/run_<timestamp>.log`, and the R and Python
steps write into that same log.

What each step consumes and produces:

| Step | Reads | Writes | Fails if |
|---|---|---|---|
| 1 Acquire | `VCF_URL` | `data/…chr2….vcf.gz` + `.tbi` | panel and index are checked separately; a complete panel is skipped, an incomplete one resumed, and the panel must pass `bgzip -t` |
| 2 Region subset | chr2 panel | `REGION_VCF` + `.tbi` | skipped when both already exist |
| 3 Sample subset | `REGION_VCF`, `SAMPLE_SUBSET_FILE` | `OUTDIR/analysis_input.vcf.gz` | symlinks the region VCF when no subset file is set; cohort size is counted here and passed to Step 7 |
| 4 Blocks + LD | `analysis_input.vcf.gz` | `plink/scn2a.{bed,bim,fam,blocks,blocks.det,ld}` | PLINK 2 instead of 1.9 |
| 5 Hapsample | `analysis_input.vcf.gz` | `hapsample/scn2a.{hap.gz,sample}` | |
| 6 Parse (R) | `.blocks.det`, `.ld`, hapsample, metadata | per-block TSVs, figures, Table S1 | `Rscript` or any R package missing |
| 7 ASOG table (Python) | `analysis_input.vcf.gz`, block TSVs, reference + `.fai` | Table S2, `all_variants_summary.tsv`, `asog_input_table.tsv`, `flagged_window_variants.tsv` | samtools is not runnable at the configured path, or the reference is not indexed — both checked upfront |

Step 7 re-counts genotypes from the VCF itself rather than reusing PLINK's, so Table
S2 and Table S1 are produced by two independent passes over the same data. Comparing
their het counts for a shared SNP is a useful consistency check.

---

## 6. A small-subset run

Use this to prove the toolchain works before committing to the full cohort. It keeps
the real region and reference but cuts the cohort to 100 samples, which brings each
step down to seconds or minutes.

**6.1 — Build a 100-sample list, 20 per superpopulation**, from the metadata you
already staged:

```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="Sample name")n=i; if($i=="Superpopulation code")s=i}; next}
            $s!~/,/ && c[$s]++ < 20 {print $n}' \
    data/igsr-1000genomes_sampleMap_30x_grch38-samples.tsv > data/smoke_100.txt
wc -l data/smoke_100.txt      # expect 100
```

The `$s!~/,/` test drops the handful of samples the portal assigns to more than one
superpopulation. Cross-check the IDs against the panel if a run reports fewer samples
than you asked for — `bcftools view --samples-file` is quiet about IDs it cannot find
unless you add `--force-samples`.

**6.2 — Write a config for it.** Copy the real one and override five values, so the
subset run cannot overwrite a full-cohort result:

```bash
cp scn2a_haplotype_config.sh config_smoke.sh
cat >> config_smoke.sh <<'EOF'

# --- small-subset overrides (appended; these win, the file is sourced top to bottom) ---
OUTDIR="results/smoke_100"
REGION_VCF="results/pipeline_SCN2A_region.vcf.gz"   # shared, cohort-independent
SAMPLE_SUBSET_FILE="data/smoke_100.txt"
ASOG_MIN_HET=25        # ~25% of 100, matching 800/3202 and 626/2504
VIZ_MIN_HET=5          # 100 would leave every figure empty at this cohort size
EOF
```

`VIZ_MIN_HET` matters: it is an absolute het count, so leaving it at 100 with a
100-sample cohort silently produces no haplotype figures at all.

**6.3 — Run it:**

```bash
bash scn2a_haplotype_pipeline.sh config_smoke.sh
```

**6.4 — Check it worked.** The log ends with a cohort size and block count; then:

```bash
ls results/smoke_100/haplotype_blocks/           # per-block TSVs, one set per block >=5kb
wc -l results/smoke_100/haplotype_blocks/all_blocks_haplotype_summary.tsv
wc -l results/smoke_100/asog_input_table.tsv
cut -f1,12 results/smoke_100/asog_input_table.tsv | head   # Rank and RefSeq
```

The things worth looking at, in order:

- **Blocks were found.** The log's `Blocks >= 5000 bp:` line being `0` means Step 6
  wrote nothing; with only 100 samples that can be real (MAF 0.15 is harder to clear),
  so retry with a larger subset before suspecting the code.
- **`flagged_window_variants.tsv` is empty or absent.** Variants whose 31 bp window
  would run past the end of the contig are listed there rather than scored, with the
  reason. A long list means the reference and the panel disagree about the contig,
  not that those variants are near the end of chromosome 2.
- **Block SNP IDs matched.** Step 6 normalises hapsample IDs (`chr` prefix stripped,
  `_` to `:`) to match PLINK's `.blocks.det`. If a panel's ID column is formatted
  differently, blocks come out with no SNPs and the per-block TSVs are empty while
  nothing errors.

Subset numbers are not the published numbers and are not meant to be: blocks are
called from the subset's own LD, so both the blocks and the ranking will differ from
the 3,202- and 2,504-sample runs.

---

## 7. Running Step 7 alone, with no external tools

`tests/smoke_step7.py` runs the Step 7 script over a tiny synthetic fixture
(six samples, five variants, one block, a 500 bp contig named `chr2`) and checks the
tables it writes. It needs nothing but `python3` — no panel, no reference, no conda
environment:

```bash
python3 tests/smoke_step7.py
```

**Verified:** this passes on the current `main`. It is the quickest way to tell
whether a change to the Python step altered the genotype counts, the superpopulation
stratification, the block join, the window arithmetic or the ranking.

What it does **not** reach matters as much, because it needs no external tools only by
standing in for them. The test replaces `fetch_sequence` with a pure-Python `faidx` of
the fixture FASTA and stubs `check_samtools`, so nothing in it runs samtools at all. A
green run therefore says nothing about whether your configured `samtools` path works,
whether your reference is faidx-indexed, or whether real sequence windows come back
correct — the first two fail fast in a real run, and §6 is what exercises the third.
It does read a real `.fai` for contig lengths, since the fixture ships one.

---

## 8. Known rough edges

None of these was ever a problem on the authors' server; each one was in the way of a
run anywhere else. Three have since been fixed on `main` and their entries say so,
because a checkout or a run that predates the fix still behaves the old way, and the
failure is easier to recognise than to diagnose. The rest are open.

**8.1 — Step 7 used to die looking for samtools. Fixed on `main`.** The Python step
called `/data0/prod_archive/software/samtools/samtools-1.9/samtools` directly,
ignoring the config's `samtools` variable. **Verified** before the fix: off that
server the call raised `FileNotFoundError` and Step 7 exited non-zero, so under
`set -e` the pipeline aborted there — Steps 1–6 and Table S1 complete, Table S2, the
all-variant summary and the ASOG table never written. (The `returncode != 0` branch
below the call only caught a samtools that existed but errored; a missing binary never
reached it.) The step now takes `--samtools` from the config, checks that samtools runs
and that the reference is indexed before it scores anything, and records a window it
cannot read in `flagged_window_variants.tsv` rather than scoring it zero. If you see
that absolute path in a traceback, you are running a copy from before the fix.

**8.2 — The configured download URL could not be fetched. Fixed on `main`.** `VCF_URL`
held a `urldefense.com` rewrite left by a mail gateway, with a trailing space, so both
of Step 1's downloads failed — the panel and `"${VCF_URL}.tbi"`; both used `wget` then,
and Step 1 uses `curl` now (§3). It now holds the
plain EBI address under the 1kGP FTP directory linked in §3. Staging the panel into
`data/` by hand is still the quicker route if you already have it: Step 1 checks the panel and
its index separately and skips either one that is already there and complete, so a staged panel
with no index costs you the index alone (§3).

**8.3 — `Rscript` is not in the dependency check.** The loop at
`scn2a_haplotype_pipeline.sh:56` validates `bcftools`, `plink`, `samtools`, `tabix`,
`bgzip`, `curl` and `python3`. A machine without R, or missing one R package, still gets through
Steps 1–5 and fails at Step 6, after the expensive subsetting work. This one is open.

**8.4 — `samtools` used to be named three different ways; it is named once now. Fixed
on `main`.** The config's `samtools` variable, a bare `samtools` on `PATH` in the
pipeline, and the absolute 1.9 path inside the Python step all had to resolve, and the
three could be different builds. Every use now goes through the config variable: the
dependency check and the reference indexing read `"${samtools}"`, and Step 7 receives
it as `--samtools` (§1).

**8.5 — The baseline ASO SNPs used to be defined twice; `BASELINE_ASOS` is now the
only place.** The R step previously hardcoded the pair as
`opt_ref_snps <- c("2:165375124:A:G", "2:165322227:A:G")` for Table S1's
`ASOtargetIsAlt` column and again inside the figure shape scale, so changing the
config moved Step 7's `UniqueSamplesAdded` while Table S1 and the figures kept
describing the old pair. The pipeline now passes `BASELINE_ASOS` to
`scn2a_parse_haplotypes.R` as an eighth argument, and both steps match a baseline to
a variant on `CHROM:POS` rather than on the full `CHROM:POS:REF:ALT` ID. Two
consequences worth knowing: the R step called by hand with only its first seven
arguments still falls back to the manuscript's pair, so old invocations behave as
before; and because matching is now positional in both steps, a baseline entry whose
alleles differ from the panel's still matches, where the R step would previously have
missed it.

**8.6 — Reading the output columns.** Three names do not mean what they look like.
**Verified** against the fixture in §7:

- `WINDOW_SIZE=30` yields **31 bp** of reference sequence, not 30: the window is
  `pos ± 15`, with the variant at the centre (`RefSeq[15]`), extended by
  `len(ref) - 1` for indels.
- `DeletionSize` is `len(alt) - len(ref)`, so it is **positive for insertions** and
  zero for every SNV.
- `Sequential` in Table S2 is always written empty.

`asog_input_table.tsv` is ranked by `UniqueSamplesAdded` descending, then
`AlleleSpecificScore`, then `HetSamples` — so a variant reaching fewer new patients
never outranks one reaching more, whatever its sequence properties. `block_id` is
written as the literal `None` for variants outside any block.

**8.7 — `FILTER` is not actually filtered.** Step 7's log line and docstrings say
"PASS variants", but the parser skips only multiallelic records and never reads the
`FILTER` column. With the 1kGP phased panel, whose records are already filtered, the
two are the same thing; with any other input they are not.

---

## 9. What a finished run should look like

From the log: the cohort size (3,202 or 2,504, or your subset), the block count, and
no `WARNING` lines about a missing LD file, a missing `D'` column, or a baseline ASO
not found in the VCF. That last one means `BASELINE_ASOS` and the region or panel
disagree, and it makes every `UniqueSamplesAdded` value wrong while the run still
reports success.

From the outputs: `all_blocks_haplotype_summary.tsv` has one row per haplotype per
block with per-superpopulation Hom/Het/Total; `snp_summary_table.tsv` covers only
SNPs inside blocks ≥ 5 kb, so it is much shorter than `all_variants_summary.tsv`;
`asog_input_table.tsv` is ranked from 1 with non-empty `RefSeq` and `AltSeq`;
`flagged_window_variants.tsv` is absent or holds only a header, since it collects the
variants whose window could not be built (§6.4); and each block ≥ 5 kb has a TSV, two
PDFs, two SVGs and a PNG, except where no haplotype cleared `VIZ_MIN_HET`.
