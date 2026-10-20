#!/usr/bin/env Rscript
# =============================================================================
# scn2a_parse_haplotypes.R
#
# Loads PLINK haplotype block output and hapsample genotype matrix.
# For each block >= min_block_bp:
#   - Identifies informative SNPs (MAF >= threshold, already filtered by PLINK)
#   - Collapses per-sample allelic patterns into haplotype strings
#   - Computes hom/het/total counts by superpopulation
#   - Saves per-block TSV
#   - Generates haplotype pattern plot and D' LD heatmap
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(forcats)
    library(viridis)
    library(gridExtra) 
    library(patchwork)
})

# ---------------------------------------------------------------------------
# Arguments — read from command line positional args or environment variables
# Usage: Rscript scn2a_parse_haplotypes.R \
#            <plink_prefix> <hapsample_prefix> <outdir> \
#            <min_block_bp> <viz_min_het> <metadata_file> <log> \
#            [baseline_asos]
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)

# Baseline ASO target SNPs, in the config's own "rsID:CHROM:POS" format, space
# separated. Defaults to the manuscript's pair, so calling this script with the
# first seven arguments behaves exactly as it did before.
DEFAULT_BASELINE_ASOS <- "rs1368238:chr2:165375124 rs72874313:chr2:165322227"

opt <- list(
    plink_prefix     = if (length(args) >= 1) args[1] else Sys.getenv("PLINK_PREFIX"),
    hapsample_prefix = if (length(args) >= 2) args[2] else Sys.getenv("HAPSAMPLE_PREFIX"),
    outdir           = if (length(args) >= 3) args[3] else Sys.getenv("OUTDIR"),
    min_block_bp     = as.integer(if (length(args) >= 4) args[4] else Sys.getenv("MIN_BLOCK_BP", "5000")),
    viz_min_het      = as.integer(if (length(args) >= 5) args[5] else Sys.getenv("VIZ_MIN_HET", "100")),
    metadata_file    = if (length(args) >= 6) args[6] else Sys.getenv("METADATA_FILE"),
    log              = if (length(args) >= 7) args[7] else Sys.getenv("LOG_FILE", ""),
    baseline_asos    = if (length(args) >= 8 && nchar(args[8]) > 0) args[8]
                       else Sys.getenv("BASELINE_ASOS", DEFAULT_BASELINE_ASOS)
)

log_msg <- function(msg) {
    ts <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
    cat(ts, msg, "\n", sep=" ")
    if (nchar(opt$log) > 0) cat(ts, msg, "\n", sep=" ", file=opt$log, append=TRUE)
}

# ---------------------------------------------------------------------------
# Baseline ASO target SNPs
# ---------------------------------------------------------------------------
# The config gives these as "rsID:CHROM:POS", which carries no alleles, while hap
# SNP IDs are CHR:POS:REF:ALT. Both are reduced to a "CHR:POS" key so that the two
# can be matched on position.

# "2:165375124:A:G" or "chr2:165375124:A:G" -> "2:165375124"
snp_pos_key <- function(snp_id) {
    parts <- strsplit(sub("^chr", "", snp_id), ":", fixed=TRUE)
    vapply(parts, function(p)
               if (length(p) >= 2) paste(p[1], p[2], sep=":") else NA_character_,
           character(1), USE.NAMES=FALSE)
}

# "rs1368238:chr2:165375124 rs72874313:chr2:165322227" -> c("2:165375124", "2:165322227")
parse_baseline_asos <- function(spec) {
    entries <- strsplit(trimws(spec), "[[:space:]]+")[[1]]
    entries <- entries[nchar(entries) > 0]
    if (length(entries) == 0) return(character(0))
    keys <- vapply(entries, function(entry) {
        parts <- strsplit(entry, ":", fixed=TRUE)[[1]]
        if (length(parts) < 3) return(NA_character_)
        paste(sub("^chr", "", parts[2]), parts[3], sep=":")
    }, character(1), USE.NAMES=FALSE)
    keys[!is.na(keys)]
}

baseline_keys <- parse_baseline_asos(opt$baseline_asos)
log_msg(paste("Baseline ASO targets:",
              if (length(baseline_keys) > 0) paste(baseline_keys, collapse=", ")
              else "none"))

dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

# ---------------------------------------------------------------------------
# Load PLINK block definitions
# ---------------------------------------------------------------------------
log_msg("Loading PLINK block definitions...")
blocks <- fread(paste0(opt$plink_prefix, ".blocks.det"))
colnames(blocks) <- c("CHR", "BP1", "BP2", "KB", "NSNPS", "SNPS")
blocks[, BLOCK_ID := .I]
blocks[, BLOCK_SIZE_BP := (BP2 - BP1)]

log_msg(paste("  Total blocks:", nrow(blocks)))
blocks_large <- blocks[BLOCK_SIZE_BP >= opt$min_block_bp]
log_msg(paste("  Blocks >=", opt$min_block_bp, "bp:", nrow(blocks_large)))

# ---------------------------------------------------------------------------
# Load hapsample matrix
# ---------------------------------------------------------------------------
log_msg("Loading hapsample matrix...")
hap_file <- paste0(opt$hapsample_prefix, ".hap.gz")
sam_file <- paste0(opt$hapsample_prefix, ".sample")

# Load hap matrix — space-delimited, no header
# Columns: CHROM ID POS REF ALT then two consecutive columns per sample (hap1, hap2)
hap <- read.table(hap_file, header=FALSE)

# .sample file: first header row is column names (ID_1 ID_2 missing),
# second row is type placeholder (0 0 0), then sample data.
# read.table with header=TRUE reads row 1 as header, row 2 (0 0 0) as first data row.
# Drop the placeholder row where ID_1 == "0".
genotype_samples <- read.table(sam_file, header=TRUE)
genotype_samples <- genotype_samples[genotype_samples$ID_1 != "0", ]
sample_ids <- genotype_samples$ID_1
n_samples <- length(sample_ids)

# Name columns: CHROM SNP POS REF ALT then SAMPLEID.1 SAMPLEID.2 per sample
# (matching original format_genotypes convention)
sample_cols <- paste(rep(sample_ids, each=2), rep(1:2, times=n_samples), sep=".")
colnames(hap) <- c("CHROM", "SNP", "POS", "REF", "ALT", sample_cols)

# Clean SNP ID to match PLINK format: strip chr prefix, replace _ with :
hap$SNP <- gsub("^chr", "", hap$SNP)
hap$SNP <- gsub("_", ":", hap$SNP)

log_msg(paste("  Variants in hapsample:", nrow(hap)))
log_msg(paste("  Samples:", n_samples))

# ---------------------------------------------------------------------------
# Load LD matrix
# ---------------------------------------------------------------------------
log_msg("Loading LD (D') data...")
ld_file <- paste0(opt$plink_prefix, ".ld")
if (file.exists(ld_file)) {
    ld <- fread(ld_file)
} else {
    log_msg("  WARNING: LD file not found, skipping LD heatmaps")
    ld <- NULL
}

# ---------------------------------------------------------------------------
# Load sample metadata (superpopulations) from IGSR sample map
# ---------------------------------------------------------------------------
# Columns: Sample name, Sex, Biosample ID, Population code, Population name,
#          Superpopulation code, Superpopulation name, Population elastic ID,
#          Data collections
log_msg("Loading sample metadata...")
meta_raw <- fread(opt$metadata_file, sep="\t", header=TRUE, quote="")
# Standardise column names by stripping spaces
setnames(meta_raw, trimws(colnames(meta_raw)))
meta <- meta_raw[, .(
    ID       = get("Sample name"),
    POP      = get("Population code"),
    SUPERPOP = get("Superpopulation code")
)]
superpops <- c("AFR","AMR","EAS","EUR","SAS")
log_msg(paste("  Metadata loaded:", nrow(meta), "samples"))

# ---------------------------------------------------------------------------
# Process each large block
# ---------------------------------------------------------------------------
all_block_summaries <- list()

for (i in seq_len(nrow(blocks_large))) {
    blk <- blocks_large[i]
    block_num <- blk$BLOCK_ID
    bp1 <- blk$BP1
    bp2 <- blk$BP2
    size_kb <- round(blk$BLOCK_SIZE_BP / 1000, 1)

    log_msg(paste0("  Block ", block_num, ": ", bp1, "-", bp2, " (", size_kb, " kb)"))

    # Get block SNP IDs from PLINK blocks.det SNPS column (pipe-delimited)
    blk_snp_ids <- unlist(strsplit(blk$SNPS, "|", fixed=TRUE))

    # Subset hap matrix to block SNPs by SNP ID
    filt <- hap[hap$SNP %in% blk_snp_ids, ]
    if (nrow(filt) == 0) {
        log_msg(paste0("    No SNPs found in block ", block_num, ", skipping"))
        next
    }
    log_msg(paste0("    SNPs in block: ", nrow(filt)))

    # Set SNP as rownames, drop non-genotype columns
    rownames(filt) <- filt$SNP
    filt_geno <- subset(filt, select = -c(CHROM, SNP, POS, REF, ALT))

    # Transpose: rows = sample.copy, cols = SNPs
    filt_t <- t(filt_geno)

    # Collapse each row to a haplotype pattern string (e.g. "0-1-0-0")
    pattern <- apply(filt_t, 1, function(row) paste(row, collapse="-"))

    # Annotate with sample and copy number
    summary_df <- data.frame(
        pattern = pattern,
        stringsAsFactors = FALSE
    )
    summary_df$id <- rownames(summary_df)

    # Split id into sample and copy (e.g. HG00096.1 -> HG00096, 1)
    summary_df$sample <- sub("\\.([12])$", "", summary_df$id)
    summary_df$copy   <- sub("^.*\\.([12])$", "\\1", summary_df$id)

    # Add block metadata
    summary_df$HAPID        <- paste0("HapGrp", sprintf("%02d", block_num))
    summary_df$CHR          <- blk$CHR
    summary_df$BP1          <- bp1
    summary_df$BP2          <- bp2
    summary_df$KB           <- blk$KB
    summary_df$NSNPS        <- nrow(filt)
    summary_df$genotypeSNPs <- paste(filt$SNP, collapse="|")

    # Join population metadata
    summary_df <- merge(summary_df,
                        meta[, c("ID","POP","SUPERPOP")],
                        by.x="sample", by.y="ID", all.x=TRUE)

    # Compute hom/het/total counts per pattern per superpopulation
    # A sample is het for a pattern if it appears in exactly one copy
    counts_by_sp <- function(df, group_col) {
        result <- df %>%
            dplyr::group_by(HAPID, CHR, BP1, BP2, KB, NSNPS, genotypeSNPs,
                            pattern, sample, !!rlang::sym(group_col)) %>%
            dplyr::summarize(diplotype = ifelse(dplyr::n() == 2, "hom", "het"),
                             .groups="drop") %>%
            dplyr::group_by(HAPID, CHR, BP1, BP2, KB, NSNPS, genotypeSNPs,
                            pattern, !!rlang::sym(group_col)) %>%
            dplyr::summarize(hom   = sum(diplotype=="hom"),
                             het   = sum(diplotype=="het"),
                             total = dplyr::n_distinct(sample),
                             .groups="drop")
        result
    }

    # Compute total counts across all superpops
    total_counts <- counts_by_sp(
        summary_df %>% dplyr::mutate(SUPERPOP="Total"), "SUPERPOP"
    ) %>%
        dplyr::rename(Hom=hom, Het=het, Total=total) %>%
        dplyr::select(-SUPERPOP)

    # Compute per-superpop counts and pivot to wide format
    sp_counts <- counts_by_sp(summary_df, "SUPERPOP") %>%
        tidyr::pivot_wider(
            names_from  = SUPERPOP,
            values_from = c(hom, het, total),
            names_glue  = "{SUPERPOP}_{.value}",
            values_fill = 0
        ) %>%
        dplyr::select(-HAPID, -CHR, -BP1, -BP2, -KB, -NSNPS, -genotypeSNPs)

    totals_superpop <- dplyr::left_join(total_counts, sp_counts, by="pattern") %>%
        dplyr::arrange(dplyr::desc(Total))

    # Build output table per block: one row per haplotype pattern
    blk_df <- totals_superpop
    blk_df$block_id       <- block_num
    blk_df$snp_ids        <- paste(filt$SNP, collapse=";")
    blk_df$snp_positions  <- paste(filt$POS, collapse=";")
    blk_df$snp_refs       <- paste(filt$REF, collapse=";")
    blk_df$snp_alts       <- paste(filt$ALT, collapse=";")

    # ---------------------------------------------------------------------------
    # Table S1 additional columns
    # ---------------------------------------------------------------------------
    # haplotypeID: block_num.rank (e.g. 1.1, 1.2) ordered by Total descending
    blk_df <- blk_df %>%
        dplyr::mutate(
            haplotypeID = paste0(block_num, ".", dplyr::row_number()),
            HapBlock    = block_num,
            chrom       = blk$CHR,
            blockStart  = bp1,
            blockEnd    = bp2,
            KB          = blk$KB,
            n_SNPs      = nrow(filt),
            # n_SNPsAlt: number of alt alleles (1s) in this haplotype's pattern
            n_SNPsAlt   = sapply(pattern, function(p)
                              sum(as.integer(strsplit(p, "-")[[1]]))),
            # ASOtargetIsAlt: TRUE if any ASO target SNP carries alt allele in this haplotype
            ASOtargetIsAlt = sapply(pattern, function(p) {
                alleles <- as.integer(strsplit(p, "-")[[1]])
                idx <- which(snp_pos_key(filt$SNP) %in% baseline_keys)
                length(idx) > 0 && any(alleles[idx] == 1, na.rm=TRUE)
            }),
            RefAltAllelePattern = pattern,
            SNPs = paste(filt$SNP, collapse="|")
        )
    outfile <- file.path(opt$outdir,
                         sprintf("block_%02d_%d_%d.tsv", block_num, bp1, bp2))
    write.table(blk_df, outfile, sep="\t", row.names=FALSE, quote=FALSE)
    log_msg(paste0("    Saved: ", basename(outfile)))
    log_msg(paste0("    Distinct haplotypes: ", length(unique(blk_df$pattern))))
    log_msg(paste0("    Haplotypes with total >= ", opt$viz_min_het, ": ",
                   sum(blk_df$Total >= opt$viz_min_het, na.rm=TRUE)))

    hap_combined <- NULL
    p_ld         <- NULL

    # -----------------------------------------------------------------------
    # Haplotype pattern plot
    # colored by allele, filtered to total >= viz_min_het)
    # -----------------------------------------------------------------------
    to_plot <- blk_df %>%
        dplyr::filter(!is.na(Total) & Total >= opt$viz_min_het) %>%
        dplyr::mutate(
            AllelePatternID = dplyr::row_number(),
            HapPattern = paste0("Hap", block_num, ".", AllelePatternID)
        ) %>%
        dplyr::mutate(
            refalt  = strsplit(pattern, "-"),
            SNP     = strsplit(snp_ids, ";"),
            pos_str = strsplit(snp_positions, ";"),
            ref_str = strsplit(snp_refs, ";"),
            alt_str = strsplit(snp_alts, ";")
        ) %>%
        tidyr::unnest(c(refalt, SNP, pos_str, ref_str, alt_str)) %>%
        dplyr::mutate(
            pos     = as.integer(pos_str),
            allele  = ifelse(refalt == "0", "reference",
                      ifelse(nchar(ref_str) != nchar(alt_str), "indel", alt_str)),
            aso_snp = ifelse(snp_pos_key(SNP) %in% baseline_keys, SNP, "not ASO SNV")
        ) %>%
        dplyr::mutate(
            allele  = factor(allele, levels=c("A","T","G","C","indel","reference")),
            HapPattern = forcats::fct_reorder(HapPattern, Total)
        )
    
    hap_start = to_plot %>% pull(pos) %>% min()
    hap_end = to_plot %>% pull(pos) %>% max()
    
    if (nrow(to_plot) > 0) {
      allele_colors <- c("A"="green3","T"="firebrick1","G"="gold",
                         "C"="dodgerblue2","indel"="black","reference"="grey81")
  
      # Build per-haplotype hom/het counts aligned to plot order
      hap_order <- levels(to_plot$HapPattern)  # ordered by fct_reorder (ascending total)
      table_data <- to_plot %>%
          dplyr::distinct(HapPattern, .keep_all=TRUE) %>%
          dplyr::select(HapPattern, Hom, Het) %>%
          dplyr::rename(Haplotype=HapPattern) %>%
          dplyr::mutate(Haplotype=factor(Haplotype, levels=hap_order)) %>%
          dplyr::arrange(Haplotype)
  
      # One diamond per baseline ASO SNP actually present in this block, so the
      # scale follows the configured pair rather than naming it.
      aso_levels <- setdiff(unique(to_plot$aso_snp), "not ASO SNV")
      aso_shapes <- c(stats::setNames(rep(23, length(aso_levels)), aso_levels),
                      "not ASO SNV" = 21)

      p <- ggplot(to_plot, aes(x=pos, y=HapPattern, group=HapPattern)) +
          geom_line(color="grey81") +
          geom_point(aes(fill=allele, shape=aso_snp), size=3,
                     position=position_jitter(width=5, height=0)) +
          scale_fill_manual(values=allele_colors, name="Allele") +
          scale_shape_manual(
              values=aso_shapes,
              name="ASO SNV") +
          guides(fill=guide_legend(override.aes=list(shape=21))) +
          xlab("Position on chromosome 2") +
          ylab("Haplotype Pattern") +
          ggtitle(paste0("Haplotype Block ", block_num,
                         " chr2:", hap_start, "-", hap_end, 
                         " (", size_kb, " kb)") ) +
          theme_minimal(base_size=12) +
          theme(
              plot.title = element_text(size = 12, face = "bold"),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              axis.title.x=element_text(margin=margin(t=10)),
              axis.title.y=element_blank(),
          )

      # Build side table as ggplot to guarantee y-axis alignment
      table_data_plot <- table_data %>%
          dplyr::mutate(
            Haplotype = factor(Haplotype, levels=hap_order),
            Hom = as.character(Hom),
            Het = as.character(Het)
          ) %>%
          tidyr::pivot_longer(cols=c(Hom, Het),
                              names_to="col", values_to="val") %>%
          dplyr::mutate(col = factor(col, levels=c("Haplotype","Hom","Het")))
  
      side_table_plot <- ggplot(table_data_plot,
                                aes(x=col, y=Haplotype, label=val)) +
          geom_text(size=3, hjust=0.5) +
          # Column headers
          annotate("text",
                   x=c("Hom","Het"),
                   y=length(hap_order) + 0.7,
                   label=c("Hom","Het"),
                   fontface="bold", size=3, hjust=0.5) +
          scale_x_discrete(limits = c("Hom", "Het"), expand=expansion(add=c(1.2,0.5))) +
          scale_y_discrete(limits=hap_order) +
          coord_cartesian(clip="off") +
          theme_void() +
          theme(
              plot.margin=margin(0, 0, 0, 0)
          )
  
      combined <- patchwork::wrap_plots(
          side_table_plot, p,
          ncol=2, widths=c(0.9, 4)
      )
      hap_combined <- combined  # store for later panel combination
  
      plot_file <- file.path(opt$outdir,
                             sprintf("block_%02d_haplotype_plot.pdf", block_num))
      plot_svg <- file.path(opt$outdir,
                             sprintf("block_%02d_haplotype_plot.svg", block_num))
      ggplot2::ggsave(plot_file, combined,
                      width  = max(8.5, nrow(filt) * 0.2),
                      height = max(4, length(unique(to_plot$HapPattern)) * 0.4))
      ggplot2::ggsave(plot_svg, combined,
                      width  = max(8.5, nrow(filt) * 0.2),
                      height = max(4, length(unique(to_plot$HapPattern)) * 0.4))
      log_msg(paste0("    Haplotype plot: ", basename(plot_file)))
    }

    # -----------------------------------------------------------------------
    # LD heatmap (D')
    # -----------------------------------------------------------------------
    # Define a function to keep every 5th break
    every_nth <- function(x, n=2) {
      if (length(x) > 50) {
      return(x[seq(1, length(x), by = n)])
      } else {
      return(x)
      }
    }
    

    if (!is.null(ld) && nrow(filt) > 1) {
        # Detect D' column name — PLINK uses DP or D_prime depending on version
        dp_col <- intersect(c("DP","D_prime","Dprime"), colnames(ld))
        if (length(dp_col) > 0) {
            dp_col <- dp_col[1]
            # restrict to happlot variants (informative SNPs)
            blk_positions <- filt$POS
            label_size = ifelse(length(blk_positions)>50, 6.5, 10)
            ld_blk <- ld[ld$BP_A %in% blk_positions &
                         ld$BP_B %in% blk_positions, ]
             
            if (nrow(ld_blk) > 0) {
                ld_plot_df <- data.frame(
                    PosA   = factor(ld_blk$BP_A),
                    PosB   = factor(ld_blk$BP_B),
                    Dprime = ld_blk[[dp_col]]
                )
                p_ld <- ggplot(ld_plot_df, aes(x=PosA, y=PosB, fill=Dprime)) +
                    geom_tile() +
#                    scale_y_discrete(breaks = every_nth) # reduce labels to every nth position
                    scale_fill_viridis_c(limits=c(0,1), option="H",
                                         direction=1, name="D'") +
                    theme_minimal(base_size=12) +
                    theme(axis.text.x=element_blank(),
                          axis.title.x=element_blank(),
                          axis.ticks.x=element_blank(),
                          axis.text.y=element_text(size=label_size),
                          panel.grid.major=element_blank(),
                          legend.position=c(0.9,0.3),
                  panel.grid.minor=element_blank()) +
                    labs(x="Position", y="Position on Chr2")  # no title on LD plot

                ld_file_out <- file.path(opt$outdir,
                                          sprintf("block_%02d_ld_heatmap.pdf", block_num))
                ld_svg_out <- file.path(opt$outdir,
                                          sprintf("block_%02d_ld_heatmap.svg", block_num))
                ggsave(ld_file_out, p_ld,
                       width=6,
                       height=6)
                ggsave(ld_svg_out, p_ld,
                       width=6,
                       height=6)
                log_msg(paste0("    LD heatmap: ", basename(ld_file_out)))
            }
        } else {
            log_msg("    WARNING: D' column not found in LD file, skipping LD heatmap")
        }
    }

    # -----------------------------------------------------------------------
    # Combined PNG: haplotype plot (A) + LD heatmap (B)
    # -----------------------------------------------------------------------
    if (!is.null(hap_combined) && !is.null(p_ld)) {
        combined_ab <- patchwork::wrap_plots(
            hap_combined, p_ld,
            ncol=1, heights=c(2, 1)
        ) +
        patchwork::plot_annotation(tag_levels="A")

        hap_w  <- max(8.5, nrow(filt) * 0.2)
        combined_png <- file.path(opt$outdir,
                                  sprintf("block_%02d_combined.png", block_num))
        ggplot2::ggsave(combined_png, combined_ab,
                        width  = hap_w,
                        height = hap_w * 0.75,
                        dpi    = 300)
        log_msg(paste0("    Combined PNG: ", basename(combined_png)))
    }

    all_block_summaries[[length(all_block_summaries)+1]] <- blk_df
}

# ---------------------------------------------------------------------------
# Write combined Table S1 equivalent
# ---------------------------------------------------------------------------
if (length(all_block_summaries) > 0) {
    combined <- do.call(dplyr::bind_rows, all_block_summaries)
    write.table(combined,
                file.path(opt$outdir, "all_blocks_haplotype_summary.tsv"),
                sep="\t", row.names=FALSE, quote=FALSE)
    log_msg(paste("Combined haplotype summary saved:",
                  file.path(opt$outdir, "all_blocks_haplotype_summary.tsv")))
}

log_msg("parse_haplotypes.R complete")
