#!/usr/bin/env python3
"""
scn2a_generate_asog_table.py

Generates the ASOG input table and SNP summary table from per-block
haplotype TSVs and the analysis VCF.

Outputs:
    asog_input_table.tsv    - Ranked table of ASO candidate SNPs
    snp_summary_table.tsv   - Full SNP-level summary (Table S2 equivalent)
"""

import os
import sys
import gzip
import glob
import math
import argparse
import subprocess
from collections import defaultdict

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--vcf",           required=True)
    p.add_argument("--blocks_dir",    required=True)
    p.add_argument("--ref_fasta",     required=True)
    p.add_argument("--samtools",      default="samtools",
                   help="samtools executable (path or name on PATH). "
                        "The pipeline passes the value set in the config file.")
    p.add_argument("--region_chrom",  required=True)
    p.add_argument("--region_start",  type=int, required=True)
    p.add_argument("--region_end",    type=int, required=True)
    p.add_argument("--window_size",   type=int, default=30)
    p.add_argument("--asog_min_het",  type=int, default=800)
    p.add_argument("--cohort_size",   type=int, required=True)
    p.add_argument("--baseline_asos", type=str, required=True,
                   help="Space-separated list of rsID:CHROM:POS for baseline ASOs")
    p.add_argument("--metadata_file", required=True,
                   help="IGSR sample map TSV: columns 'Sample name' and 'Superpopulation code'")
    p.add_argument("--gc_min",        type=float, default=0.40)
    p.add_argument("--gc_max",        type=float, default=0.60)
    p.add_argument("--outdir",        required=True)
    p.add_argument("--log",           default="")
    return p.parse_args()


def log_msg(msg, logfile=""):
    from datetime import datetime
    ts = datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
    print(f"{ts} {msg}", flush=True)
    if logfile:
        with open(logfile, "a") as f:
            f.write(f"{ts} {msg}\n")


# ---------------------------------------------------------------------------
# VCF parsing
# ---------------------------------------------------------------------------
def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_vcf_variants(vcf_path, chrom, start, end):
    """
    Parse all PASS biallelic variants in the region.
    Returns list of dicts with keys:
        chrom, pos, ref, alt, samples_het (set), samples_homalt (set),
        samples_homref (set), n_het, n_homalt, n_homref, maf
    """
    variants = []
    sample_ids = []

    with open_vcf(vcf_path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                sample_ids = line.strip().split("\t")[9:]
                continue

            fields = line.strip().split("\t")
            if len(fields) < 10:
                continue

            vcf_chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            fmt = fields[8]

            if vcf_chrom != chrom:
                continue
            if pos < start or pos > end:
                continue
            if "," in alt:  # skip multiallelic for ASOG table
                continue

            fmt_keys = fmt.split(":")
            gt_idx = fmt_keys.index("GT") if "GT" in fmt_keys else 0

            samples_het = set()
            samples_homalt = set()
            samples_homref = set()

            for i, samp_str in enumerate(fields[9:]):
                gt_raw = samp_str.split(":")[gt_idx]
                gt_clean = gt_raw.replace("|", "/")
                alleles = gt_clean.split("/")
                if "." in alleles:
                    continue
                try:
                    a = [int(x) for x in alleles]
                except ValueError:
                    continue

                sid = sample_ids[i]
                if a[0] == 0 and a[1] == 0:
                    samples_homref.add(sid)
                elif a[0] == 1 and a[1] == 1:
                    samples_homalt.add(sid)
                elif (a[0] == 0 and a[1] == 1) or (a[0] == 1 and a[1] == 0):
                    samples_het.add(sid)

            n_called = len(samples_het) + len(samples_homalt) + len(samples_homref)
            if n_called == 0:
                continue

            alt_count = len(samples_het) + 2 * len(samples_homalt)
            total_alleles = 2 * n_called
            af = alt_count / total_alleles
            maf = min(af, 1 - af)

            variants.append({
                "chrom":         vcf_chrom,
                "pos":           pos,
                "ref":           ref,
                "alt":           alt,
                "samples_het":   samples_het,
                "samples_homalt": samples_homalt,
                "samples_homref": samples_homref,
                "n_het":         len(samples_het),
                "n_homalt":      len(samples_homalt),
                "n_homref":      len(samples_homref),
                "af":            af,
                "maf":           maf,
            })

    return variants, sample_ids


# ---------------------------------------------------------------------------
# Reference sequence extraction
# ---------------------------------------------------------------------------
class SequenceFetchError(RuntimeError):
    """Raised when a reference sequence window cannot be read."""


class WindowOutOfBounds(Exception):
    """
    Raised when a variant's window would extend past either end of its contig.

    This is not a configuration error — the contig is present and readable, the
    variant simply sits too close to its end for a full window. samtools faidx
    would return a short sequence (or none) and exit 0, so these variants are
    detected from the contig length and flagged rather than scored.
    """
    def __init__(self, chrom, contig_len, win_start, win_end):
        self.chrom = chrom
        self.contig_len = contig_len
        self.win_start = win_start
        self.win_end = win_end
        super().__init__(
            f"window {chrom}:{win_start}-{win_end} extends beyond {chrom} "
            f"(length {contig_len})"
        )


def load_contig_lengths(ref_fasta):
    """
    Read contig name -> length from the faidx index (column 2 of the .fai).
    Used to detect windows that would overrun a contig end before asking
    samtools for them.
    """
    lengths = {}
    with open(ref_fasta + ".fai") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                lengths[fields[0]] = int(fields[1])
    return lengths


def check_samtools(samtools, ref_fasta):
    """
    Verify samtools is runnable and the reference is faidx-indexed, so the run
    fails at startup rather than after scoring thousands of empty windows.
    """
    try:
        result = subprocess.run([samtools, "--version"],
                                capture_output=True, text=True)
    except OSError as exc:
        raise SequenceFetchError(
            f"cannot run samtools at '{samtools}': {exc}. "
            f"Set the samtools path in the config file (or pass --samtools)."
        ) from exc
    if result.returncode != 0:
        raise SequenceFetchError(
            f"'{samtools} --version' exited {result.returncode}: "
            f"{result.stderr.strip()}"
        )
    if not os.path.exists(ref_fasta):
        raise SequenceFetchError(f"reference FASTA not found: {ref_fasta}")
    if not os.path.exists(ref_fasta + ".fai"):
        raise SequenceFetchError(
            f"reference FASTA is not faidx-indexed: {ref_fasta}.fai is missing. "
            f"Run '{samtools} faidx {ref_fasta}'."
        )
    return result.stdout.strip().split("\n")[0]


def fetch_sequence(samtools, ref_fasta, chrom, start, end):
    """
    Fetch reference sequence using samtools faidx.
    start/end are 1-based inclusive.
    Returns uppercase string of exactly (end - start + 1) bases.

    Raises SequenceFetchError if samtools cannot be run, if it exits non-zero,
    or if the window comes back empty or short. A window that cannot be read
    is a hard error: silently returning "" scores the variant as zero and the
    bad score is indistinguishable from a real one in the output table.
    """
    region = f"{chrom}:{start}-{end}"
    try:
        result = subprocess.run(
            [samtools, "faidx", ref_fasta, region],
            capture_output=True, text=True
        )
    except OSError as exc:
        raise SequenceFetchError(
            f"cannot run samtools at '{samtools}': {exc}"
        ) from exc

    if result.returncode != 0:
        raise SequenceFetchError(
            f"samtools faidx failed for {region} (exit {result.returncode}): "
            f"{result.stderr.strip()}"
        )

    lines = result.stdout.strip().split("\n")
    seq = "".join(lines[1:]).upper()
    expected = end - start + 1

    if not seq:
        raise SequenceFetchError(
            f"samtools faidx returned no sequence for {region} in {ref_fasta}. "
            f"Check that the contig name matches the reference "
            f"(e.g. 'chr2' vs '2')."
        )
    if len(seq) != expected:
        raise SequenceFetchError(
            f"samtools faidx returned {len(seq)} bases for {region}, "
            f"expected {expected}. The window may run off the end of {chrom}."
        )
    return seq


def gc_content(seq):
    if not seq:
        return None
    gc = sum(1 for b in seq.upper() if b in "GC")
    return gc / len(seq)


# ---------------------------------------------------------------------------
# Metadata loading
# ---------------------------------------------------------------------------
SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]

def load_metadata(metadata_file):
    """
    Load IGSR sample map. Returns dict: sample_id -> superpopulation code.
    Expects tab-delimited file with headers 'Sample name' and 'Superpopulation code'.
    """
    sample_superpop = {}
    with open(metadata_file) as fh:
        header = [h.strip() for h in fh.readline().split("\t")]
        try:
            name_idx = header.index("Sample name")
            sp_idx   = header.index("Superpopulation code")
        except ValueError:
            raise ValueError(
                f"metadata_file must contain 'Sample name' and "
                f"'Superpopulation code' columns. Found: {header}"
            )
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(name_idx, sp_idx):
                continue
            sid = fields[name_idx].strip()
            sp  = fields[sp_idx].strip()
            if sid and sp:
                sample_superpop[sid] = sp
    return sample_superpop


# ---------------------------------------------------------------------------
# Build ref and alt sequences for a variant
# ---------------------------------------------------------------------------
def build_sequences(samtools, ref_fasta, chrom, pos, ref, alt, window_size,
                    contig_len=None):
    """
    Build RefSeq and AltSeq centered on the variant.
    pos is 1-based.
    For SNPs: sequences are equal length (window_size).
    For indels: AltSeq may differ in length.

    Returns: (win_start, win_end, ref_seq, alt_seq)

    Raises WindowOutOfBounds if contig_len is given and the window would run off
    either end of the contig, so the caller can flag the variant instead of
    scoring a truncated window.
    """
    half = window_size // 2
    ref_len = len(ref)
    alt_len = len(alt)

    # Window centered on variant start position
    win_start = pos - half
    win_end   = pos + half + (ref_len - 1)  # extend for ref length

    # A variant too close to either contig end cannot have a full window.
    # Detect it here: samtools faidx would return a short sequence, or none,
    # and still exit 0.
    if contig_len is not None and (win_start < 1 or win_end > contig_len):
        raise WindowOutOfBounds(chrom, contig_len, win_start, win_end)

    # Fetch reference sequence for the full window.
    # Raises SequenceFetchError rather than returning an empty window.
    ref_seq_full = fetch_sequence(samtools, ref_fasta, chrom, win_start, win_end)

    # Build alt sequence by substituting alt allele at the center
    # The variant occupies positions [half, half + ref_len) in the window (0-based)
    center = half  # 0-based index of variant in window
    prefix = ref_seq_full[:center]
    suffix = ref_seq_full[center + ref_len:]
    alt_seq = prefix + alt.upper() + suffix

    return win_start, win_end, ref_seq_full, alt_seq


# ---------------------------------------------------------------------------
# AlleleSpecificScore
# ---------------------------------------------------------------------------
def allele_specific_score(ref_seq, alt_seq, ref_allele, alt_allele, window_size):
    """
    Component 1: compare center base of ref_seq vs alt_seq
        0 if same
        1 if transition (both purines or both pyrimidines)
        2 if transversion
    Component 2: +1 if insertion AND 3-base window around center differs
    """
    # Sequences are built by build_sequences(), which raises rather than
    # returning an empty window, so empty input here is a bug, not a variant
    # that happens to score zero.
    if not ref_seq or not alt_seq:
        raise SequenceFetchError(
            f"empty sequence window for {ref_allele}>{alt_allele}: "
            f"cannot compute an allele-specific score"
        )

    center_0 = window_size // 2  # 0-based center index

    # Safe access
    def safe_char(seq, idx):
        if 0 <= idx < len(seq):
            return seq[idx].upper()
        return ""

    ref_center = safe_char(ref_seq, center_0)
    alt_center = safe_char(alt_seq, center_0)

    if not ref_center or not alt_center:
        raise SequenceFetchError(
            f"sequence window is too short to read the center base "
            f"(index {center_0}; ref {len(ref_seq)} bases, "
            f"alt {len(alt_seq)} bases) for {ref_allele}>{alt_allele}"
        )

    # Component 1
    if ref_center == alt_center:
        comp1 = 0
    else:
        purines = set("AG")
        pyrimidines = set("CT")
        if (ref_center in purines and alt_center in purines) or \
           (ref_center in pyrimidines and alt_center in pyrimidines):
            comp1 = 1
        else:
            comp1 = 2

    # Component 2: insertion bonus
    is_insertion = len(alt_allele) > len(ref_allele)
    if is_insertion:
        # 3-base window around position center_0 - 1 (0-based: positions center-1, center, center+1)
        ref_tri = ref_seq[max(0, center_0-1): center_0+2]
        alt_tri = alt_seq[max(0, center_0-1): center_0+2]
        comp2 = 1 if ref_tri != alt_tri else 0
    else:
        comp2 = 0

    return comp1 + comp2


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    logfile = args.log

    log_msg("=== scn2a_generate_asog_table.py ===", logfile)
    log_msg(f"  VCF: {args.vcf}", logfile)
    log_msg(f"  Region: {args.region_chrom}:{args.region_start}-{args.region_end}", logfile)
    log_msg(f"  Cohort size: {args.cohort_size}", logfile)
    log_msg(f"  Window size: {args.window_size}", logfile)
    log_msg(f"  ASOG min het: {args.asog_min_het}", logfile)

    # -----------------------------------------------------------------------
    # Check samtools before any scoring, so a bad path fails here rather
    # than turning every allele-specific score into a zero.
    # -----------------------------------------------------------------------
    version = check_samtools(args.samtools, args.ref_fasta)
    log_msg(f"  samtools: {args.samtools} ({version})", logfile)

    contig_lengths = load_contig_lengths(args.ref_fasta)
    region_contig_len = contig_lengths.get(args.region_chrom)
    if region_contig_len is None:
        raise SequenceFetchError(
            f"contig '{args.region_chrom}' is not in {args.ref_fasta}.fai. "
            f"Check that the contig name matches the reference "
            f"(e.g. 'chr2' vs '2')."
        )
    log_msg(f"  {args.region_chrom} length: {region_contig_len:,} bp", logfile)

    os.makedirs(args.outdir, exist_ok=True)

    # -----------------------------------------------------------------------
    # Parse baseline ASO target variants
    # -----------------------------------------------------------------------
    log_msg("Parsing baseline ASO targets...", logfile)
    baseline_targets = {}  # rsid -> (chrom, pos)
    for entry in args.baseline_asos.strip().split():
        parts = entry.split(":")
        if len(parts) == 3:
            rsid, chrom, pos = parts
            baseline_targets[rsid] = (chrom, int(pos))
            log_msg(f"  Baseline: {rsid} at {chrom}:{pos}", logfile)

    # -----------------------------------------------------------------------
    # Parse all PASS variants from VCF
    # -----------------------------------------------------------------------
    log_msg("Parsing VCF variants...", logfile)
    variants, sample_ids = parse_vcf_variants(
        args.vcf, args.region_chrom, args.region_start, args.region_end
    )
    log_msg(f"  PASS biallelic variants: {len(variants)}", logfile)

    # Index variants by position
    var_by_pos = {v["pos"]: v for v in variants}

    # -----------------------------------------------------------------------
    # Get het sample sets for baseline ASO variants
    # -----------------------------------------------------------------------
    baseline_het_samples = set()
    baseline_coverage = {}
    for rsid, (chrom, pos) in baseline_targets.items():
        if pos in var_by_pos:
            het_set = var_by_pos[pos]["samples_het"]
            baseline_coverage[rsid] = het_set
            baseline_het_samples |= het_set
            log_msg(f"  {rsid}: {len(het_set)} het samples", logfile)
        else:
            log_msg(f"  WARNING: {rsid} at pos {pos} not found in VCF", logfile)
            baseline_coverage[rsid] = set()

    log_msg(f"  Baseline ASOs cover {len(baseline_het_samples)} unique het samples", logfile)

    # -----------------------------------------------------------------------
    # Load metadata
    # -----------------------------------------------------------------------
    log_msg("Loading sample metadata...", logfile)
    sample_superpop = load_metadata(args.metadata_file)
    # Build per-superpop sample sets
    sp_samples = defaultdict(set)
    for sid, sp in sample_superpop.items():
        if sp in SUPERPOPS:
            sp_samples[sp].add(sid)
    log_msg(f"  Metadata loaded: {len(sample_superpop)} samples", logfile)
    for sp in SUPERPOPS:
        log_msg(f"    {sp}: {len(sp_samples[sp])} samples", logfile)

    # -----------------------------------------------------------------------
    # Load per-block TSVs to get block membership and haplotype info per SNP
    # -----------------------------------------------------------------------
    log_msg("Loading per-block TSVs...", logfile)
    block_files = sorted(glob.glob(os.path.join(args.blocks_dir, "block_*.tsv")))

    # snp_to_block: pos -> block_id
    snp_to_block = {}
    # snp_to_haplotypes: pos -> {"alt": [hapIDs], "ref": [hapIDs]}
    snp_to_haplotypes = defaultdict(lambda: {"alt": [], "ref": []})

    for bf in block_files:
        with open(bf) as fh:
            header = fh.readline().strip().split("\t")
            rows = []
            for line in fh:
                row = dict(zip(header, line.strip().split("\t")))
                rows.append(row)

        if not rows:
            continue

        # SNP positions and IDs for this block
        snp_positions = rows[0]["snp_positions"].split(";")
        snp_ids       = rows[0]["snp_ids"].split(";")
        bid = int(rows[0]["block_id"])

        for p in snp_positions:
            if p:
                snp_to_block[int(p)] = bid

        # For each haplotype row, check each SNP position for alt (1) vs ref (0)
        for row in rows:
            hap_id  = row.get("haplotypeID", "")
            pattern = row.get("RefAltAllelePattern", row.get("pattern", ""))
            if not pattern or not hap_id:
                continue
            alleles = pattern.split("-")
            for snp_id, allele in zip(snp_ids, alleles):
                # Extract position from SNP ID (format: CHROM:POS:REF:ALT)
                parts = snp_id.split(":")
                if len(parts) < 2:
                    continue
                try:
                    pos_i = int(parts[1])
                except ValueError:
                    continue
                if allele == "1":
                    snp_to_haplotypes[pos_i]["alt"].append(hap_id)
                else:
                    snp_to_haplotypes[pos_i]["ref"].append(hap_id)

    # -----------------------------------------------------------------------
    # Build ASOG table rows for all PASS het variants >= asog_min_het
    # -----------------------------------------------------------------------
    log_msg("Building ASOG table...", logfile)
    asog_rows = []
    snp_rows = []
    flagged_variants = []

    for v in variants:
        pos = v["pos"]
        ref = v["ref"]
        alt = v["alt"]
        n_het = v["n_het"]
        het_samples = v["samples_het"]

        # Build sequences. A variant too close to a contig end cannot have a
        # full window: flag it and leave its score blank, never 0, since 0 is a
        # legitimate score. Flagged variants are excluded from the ASOG table,
        # which ranks on that score.
        window_flag = ""
        try:
            win_start, win_end, ref_seq, alt_seq = build_sequences(
                args.samtools, args.ref_fasta, args.region_chrom,
                pos, ref, alt, args.window_size,
                contig_len=contig_lengths.get(args.region_chrom)
            )
        except WindowOutOfBounds as oob:
            window_flag = "window_beyond_contig_end"
            win_start, win_end = oob.win_start, oob.win_end
            ref_seq = alt_seq = ""
            aso_score = ""
            ref_gc = alt_gc = None
            flagged_variants.append({
                "Chrom": v["chrom"], "Pos": pos, "Ref": ref, "Alt": alt,
                "Total_Het": n_het,
                "win_start": win_start, "win_end": win_end,
                "contig_length": oob.contig_len,
                "reason": str(oob),
            })
        else:
            # AlleleSpecificScore
            aso_score = allele_specific_score(ref_seq, alt_seq, ref, alt, args.window_size)

            # GC content
            ref_gc = gc_content(ref_seq)
            alt_gc = gc_content(alt_seq)

        # Unique samples added beyond baseline
        unique_added = len(het_samples - baseline_het_samples)

        # Total covered if this variant were added to the ASO panel
        samples_covered = baseline_het_samples | het_samples
        pct_covered = len(samples_covered) / args.cohort_size

        # Deletion size
        del_size = len(alt) - len(ref)

        # Block membership
        block_id = snp_to_block.get(pos, None)

        # Haplotype assignments for this SNP (only meaningful for block SNPs)
        hap_info     = snp_to_haplotypes.get(pos, {"alt": [], "ref": []})
        alt_haps     = sorted(hap_info["alt"])
        ref_haps     = sorted(hap_info["ref"])
        n_alt_haps   = len(alt_haps)
        n_ref_haps   = len(ref_haps)

        # Superpopulation-stratified counts
        sp_het    = {sp: len(v["samples_het"]    & sp_samples[sp]) for sp in SUPERPOPS}
        sp_homalt = {sp: len(v["samples_homalt"] & sp_samples[sp]) for sp in SUPERPOPS}
        sp_homref = {sp: len(v["samples_homref"] & sp_samples[sp]) for sp in SUPERPOPS}

        # SNP ID in CHROM:POS:REF:ALT format (strip chr prefix to match hap IDs)
        chrom_bare = v["chrom"].replace("chr", "")
        snp_id = f"{chrom_bare}:{pos}:{ref}:{alt}"

        # SNP summary row — Table S2 structure
        snp_row = {
            "Chrom":           v["chrom"],
            "Pos":             pos,
            "Ref":             ref,
            "Alt":             alt,
            "Sequential":      "",
            "SNP":             snp_id,
            # Total counts
            "Total_Het":       n_het,
            "Total_HomAlt":    v["n_homalt"],
            "Total_HomRef":    v["n_homref"],
            # Per-superpop counts
            "AFR_Het":         sp_het["AFR"],
            "AFR_HomAlt":      sp_homalt["AFR"],
            "AFR_HomRef":      sp_homref["AFR"],
            "AMR_Het":         sp_het["AMR"],
            "AMR_HomAlt":      sp_homalt["AMR"],
            "AMR_HomRef":      sp_homref["AMR"],
            "EAS_Het":         sp_het["EAS"],
            "EAS_HomAlt":      sp_homalt["EAS"],
            "EAS_HomRef":      sp_homref["EAS"],
            "EUR_Het":         sp_het["EUR"],
            "EUR_HomAlt":      sp_homalt["EUR"],
            "EUR_HomRef":      sp_homref["EUR"],
            "SAS_Het":         sp_het["SAS"],
            "SAS_HomAlt":      sp_homalt["SAS"],
            "SAS_HomRef":      sp_homref["SAS"],
            # Block and haplotype info
            "HapBlock":        block_id,
            "n_altHaplotypes": n_alt_haps,
            "n_refHaplotypes": n_ref_haps,
            "altHaplotypes":   ", ".join(alt_haps),
            "refHaplotypes":   ", ".join(ref_haps),
            "hetSamples":      ",".join(sorted(v["samples_het"])),
            # ASOG-relevant fields (kept for reference, not in S2)
            "maf":             round(v["maf"], 4),
            "af":              round(v["af"], 4),
            "allele_specific_score": aso_score,
            "unique_samples_added":  unique_added,
            "ref_gc":          round(ref_gc, 3) if ref_gc is not None else "",
            "alt_gc":          round(alt_gc, 3) if alt_gc is not None else "",
            "win_start":       win_start,
            "win_end":         win_end,
            "ref_seq":         ref_seq,
            "alt_seq":         alt_seq,
            "del_size":        del_size,
            "window_flag":     window_flag,
        }
        snp_rows.append(snp_row)

        # ASOG row: only variants with n_het >= asog_min_het, and only those
        # with a full sequence window — the table ranks on AlleleSpecificScore.
        if n_het >= args.asog_min_het and not window_flag:
            asog_rows.append({
                "Variant":            f"{v['chrom']}:{pos}:{ref}:{alt}",
                "AlleleSpecificScore": aso_score,
                "HetSamples":         n_het,
                "UniqueSamplesAdded": unique_added,
                "SamplesCovered":     len(samples_covered),
                "pctCohortCovered":   round(pct_covered, 4),
                "DeletionSize":       del_size,
                "win_start":          win_start,
                "win_end":            win_end,
                "RefSeq":             ref_seq,
                "AltSeq":             alt_seq,
                "RefGC":              round(ref_gc, 3) if ref_gc is not None else "",
                "AltGC":              round(alt_gc, 3) if alt_gc is not None else "",
                "block_id":           block_id,
            })

    log_msg(f"  Total PASS variants: {len(snp_rows)}", logfile)
    log_msg(f"  Variants in ASOG table (n_het >= {args.asog_min_het}): {len(asog_rows)}", logfile)

    # -----------------------------------------------------------------------
    # Report variants flagged for an incomplete window
    # -----------------------------------------------------------------------
    if flagged_variants:
        log_msg(
            f"  FLAGGED (window beyond contig end, not scored): "
            f"{len(flagged_variants)} of {len(snp_rows)} variants", logfile
        )
        for f in flagged_variants[:10]:
            log_msg(f"    {f['Chrom']}:{f['Pos']} {f['Ref']}>{f['Alt']}: {f['reason']}",
                    logfile)
        if len(flagged_variants) > 10:
            log_msg(f"    ... and {len(flagged_variants) - 10} more", logfile)

        flagged_cols = ["Chrom", "Pos", "Ref", "Alt", "Total_Het",
                        "win_start", "win_end", "contig_length", "reason"]
        flagged_out = os.path.join(args.outdir, "flagged_window_variants.tsv")
        with open(flagged_out, "w") as f:
            f.write("\t".join(flagged_cols) + "\n")
            for row in sorted(flagged_variants, key=lambda r: r["Pos"]):
                f.write("\t".join(str(row.get(c, "")) for c in flagged_cols) + "\n")
        log_msg(f"  Flagged variants written: {flagged_out}", logfile)
    else:
        log_msg("  Flagged (window beyond contig end): 0", logfile)

    # -----------------------------------------------------------------------
    # Rank ASOG rows
    # -----------------------------------------------------------------------
    # Primary: UniqueSamplesAdded descending
    # Secondary: AlleleSpecificScore descending
    # Tertiary: HetSamples descending
    asog_rows_sorted = sorted(
        asog_rows,
        key=lambda r: (-r["UniqueSamplesAdded"], -r["AlleleSpecificScore"], -r["HetSamples"])
    )
    for rank, row in enumerate(asog_rows_sorted, 1):
        row["Rank"] = rank

    # -----------------------------------------------------------------------
    # Write ASOG table
    # -----------------------------------------------------------------------
    asog_cols = [
        "Rank", "Variant", "AlleleSpecificScore", "HetSamples",
        "UniqueSamplesAdded", "SamplesCovered", "pctCohortCovered",
        "DeletionSize", "win_start", "win_end", "RefSeq", "AltSeq",
        "RefGC", "AltGC", "block_id"
    ]

    asog_out = os.path.join(args.outdir, "asog_input_table.tsv")
    with open(asog_out, "w") as f:
        f.write("\t".join(asog_cols) + "\n")
        for row in asog_rows_sorted:
            f.write("\t".join(str(row.get(c, "")) for c in asog_cols) + "\n")

    log_msg(f"  ASOG table written: {asog_out}", logfile)

    # -----------------------------------------------------------------------
    # Write Table S2 — block SNPs only, exact column structure
    # -----------------------------------------------------------------------
    s2_cols = [
        "Chrom", "Pos", "Ref", "Alt", "Sequential", "SNP",
        "Total_Het", "Total_HomAlt", "Total_HomRef",
        "AFR_Het", "AFR_HomAlt", "AFR_HomRef",
        "AMR_Het", "AMR_HomAlt", "AMR_HomRef",
        "EAS_Het", "EAS_HomAlt", "EAS_HomRef",
        "EUR_Het", "EUR_HomAlt", "EUR_HomRef",
        "SAS_Het", "SAS_HomAlt", "SAS_HomRef",
        "HapBlock", "n_altHaplotypes", "n_refHaplotypes",
        "altHaplotypes", "refHaplotypes", "hetSamples",
    ]

    # Block SNPs only: those with a block_id assigned
    block_snp_rows = [r for r in snp_rows if r["HapBlock"] is not None]
    log_msg(f"  Block SNPs for Table S2: {len(block_snp_rows)}", logfile)

    snp_out = os.path.join(args.outdir, "snp_summary_table.tsv")
    with open(snp_out, "w") as f:
        f.write("\t".join(s2_cols) + "\n")
        for row in sorted(block_snp_rows, key=lambda r: r["Pos"]):
            f.write("\t".join(str(row.get(c, "")) for c in s2_cols) + "\n")

    log_msg(f"  SNP summary table (Table S2) written: {snp_out}", logfile)

    # -----------------------------------------------------------------------
    # Write full SNP table (all variants, for internal use / ASOG context)
    # -----------------------------------------------------------------------
    full_snp_cols = [
        "Chrom", "Pos", "Ref", "Alt", "SNP", "maf", "af",
        "Total_Het", "Total_HomAlt", "Total_HomRef",
        "HapBlock", "allele_specific_score", "unique_samples_added",
        "ref_gc", "alt_gc", "win_start", "win_end",
        "ref_seq", "alt_seq", "del_size", "window_flag",
    ]

    full_snp_out = os.path.join(args.outdir, "all_variants_summary.tsv")
    with open(full_snp_out, "w") as f:
        f.write("\t".join(full_snp_cols) + "\n")
        for row in sorted(snp_rows, key=lambda r: r["Pos"]):
            f.write("\t".join(str(row.get(c, "")) for c in full_snp_cols) + "\n")

    log_msg(f"  Full variant table written: {full_snp_out}", logfile)
    log_msg("=== scn2a_generate_asog_table.py complete ===", logfile)


if __name__ == "__main__":
    try:
        main()
    except SequenceFetchError as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        sys.exit(1)
