#!/usr/bin/env python3
"""
Smoke test for Step 7 (scn2a_generate_asog_table.py).

Runs the Step 7 script over the tiny fixture in tests/fixtures/ and checks that
it produces the three expected tables with the expected ranking. The fixture is
synthetic: six samples, five variants on a 500 bp contig named chr2, one
haplotype block. It exercises the genotype counter, the superpopulation
stratification, the block/haplotype join, the sequence window and the ranking,
without needing the 1kGP panel or the GRCh38 reference.

The Step 7 script shells out to samtools faidx at a path hardcoded for the
authors' server (see REPRODUCING.md, "Known rough edges"), so this test
substitutes a pure-Python faidx of the fixture FASTA. That keeps the test
runnable anywhere with python3 alone, and keeps it independent of whether that
path is ever made configurable.

Usage:
    python3 tests/smoke_step7.py

Exits 0 on success, 1 with a message on the first failed check.
"""

import csv
import importlib.util
import pathlib
import sys
import tempfile

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent
FIXTURES = HERE / "fixtures"


def load_step7():
    """Import scn2a_generate_asog_table.py as a module."""
    path = REPO / "scn2a_generate_asog_table.py"
    spec = importlib.util.spec_from_file_location("scn2a_asog", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def python_faidx(ref_fasta, chrom, start, end):
    """Stand in for `samtools faidx <fasta> <chrom>:<start>-<end>` (1-based, inclusive)."""
    seqs, current = {}, None
    with open(ref_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                current = line[1:].strip().split()[0]
                seqs[current] = []
            elif current is not None:
                seqs[current].append(line.strip())
    seq = "".join(seqs.get(chrom, []))
    return seq[max(start - 1, 0):end].upper()


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def check(condition, message):
    if not condition:
        print(f"FAIL: {message}")
        sys.exit(1)
    print(f"  ok: {message}")


def main():
    step7 = load_step7()
    step7.fetch_sequence = python_faidx

    with tempfile.TemporaryDirectory() as outdir:
        sys.argv = [
            "scn2a_generate_asog_table.py",
            "--vcf", str(FIXTURES / "mini.vcf"),
            "--blocks_dir", str(FIXTURES / "blocks"),
            "--ref_fasta", str(FIXTURES / "mini_chr2.fasta"),
            "--region_chrom", "chr2",
            "--region_start", "1",
            "--region_end", "500",
            "--window_size", "30",
            "--asog_min_het", "2",
            "--cohort_size", "6",
            "--baseline_asos", "rsBASE1:chr2:200 rsBASE2:chr2:300",
            "--metadata_file", str(FIXTURES / "metadata.tsv"),
            "--outdir", outdir,
        ]
        step7.main()

        out = pathlib.Path(outdir)
        for name in ("asog_input_table.tsv", "snp_summary_table.tsv",
                     "all_variants_summary.tsv"):
            check((out / name).is_file(), f"{name} written")

        asog = read_tsv(out / "asog_input_table.tsv")
        s2 = read_tsv(out / "snp_summary_table.tsv")

        check(len(asog) == 5, f"all 5 fixture variants clear asog_min_het=2 (got {len(asog)})")
        check(len(s2) == 3, f"Table S2 holds only the 3 block SNPs (got {len(s2)})")

        # Ranking: UniqueSamplesAdded desc, then AlleleSpecificScore desc, then HetSamples desc.
        variants = [r["Variant"] for r in asog]
        expected = ["chr2:250:G:A", "chr2:150:A:G", "chr2:350:A:AGG",
                    "chr2:200:C:T", "chr2:300:T:C"]
        check(variants == expected, f"ranking order is {expected}")
        check([int(r["Rank"]) for r in asog] == [1, 2, 3, 4, 5], "Rank is 1..5 in table order")

        # The two baseline variants add nobody beyond themselves.
        baselines = {r["Variant"]: r for r in asog}
        check(int(baselines["chr2:200:C:T"]["UniqueSamplesAdded"]) == 0,
              "baseline variant adds 0 unique samples")
        check(int(baselines["chr2:250:G:A"]["UniqueSamplesAdded"]) == 1,
              "chr2:250 adds the 1 sample the baselines miss")

        # WINDOW_SIZE=30 means pos +/- 15, i.e. 31 bp of reference.
        snv = baselines["chr2:250:G:A"]
        check(len(snv["RefSeq"]) == 31, f"RefSeq is 31 bp for WINDOW_SIZE=30 (got {len(snv['RefSeq'])})")
        check(len(snv["AltSeq"]) == 31, "AltSeq is 31 bp for an SNV")
        check(snv["RefSeq"][15] == "G" and snv["AltSeq"][15] == "A",
              "the variant sits at the centre of the window")

        # A 2 bp insertion lengthens AltSeq and reports DeletionSize as +2.
        ins = baselines["chr2:350:A:AGG"]
        check(len(ins["AltSeq"]) == len(ins["RefSeq"]) + 2, "insertion lengthens AltSeq by 2 bp")
        check(int(ins["DeletionSize"]) == 2, "DeletionSize is +2 for a 2 bp insertion")

        # Superpopulation counts should add up to the total het count.
        row = {r["SNP"]: r for r in s2}["2:250:G:A"]
        sp_total = sum(int(row[f"{sp}_Het"]) for sp in ("AFR", "AMR", "EAS", "EUR", "SAS"))
        check(sp_total == int(row["Total_Het"]),
              f"superpopulation het counts sum to Total_Het ({sp_total})")

        # Block membership comes from the per-block TSV written by Step 6.
        check(row["HapBlock"] == "1", "block SNP is assigned to block 1")
        # Fixture patterns: 1.1 = 1-0-1, 1.2 = 0-1-0 over (150, 250, 350),
        # so 1.2 is the haplotype carrying alt at chr2:250 and 1.1 carries ref.
        check(row["altHaplotypes"] == "1.2",
              "haplotype 1.2 carries the alt allele at chr2:250")
        check(row["refHaplotypes"] == "1.1",
              "haplotype 1.1 carries the ref allele at chr2:250")

    print("\nsmoke_step7: all checks passed")


if __name__ == "__main__":
    main()
