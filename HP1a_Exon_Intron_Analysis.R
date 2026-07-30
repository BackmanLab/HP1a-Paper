
##################################################################
# Author: Lucas Carter
# Email: lucascarter2025@u.northwestern.edu
# PI: Vadim Backman
# Description:
# Final HP1a metagene analysis pipeline.
# Three conditions plotted as colored lines: Auxin, Hypoxia,
# Auxin+Hypoxia (all log2 ratios vs Control).
#
# Three plots:
#   Plot 1: All genes — no faceting
#   Plot 2: Faceted by expression group
#           (Low+Unexpressed merged, Medium, High)
#   Plot 3: Faceted by E/I group (Low/Medium/High E/I)
#
# Each plot has three stacked panels:
#   Panel 1: TSS-TES gene body (scale-regions: -2kb | body | +2kb)
#   Panel 2: Exon metagene (5' to 3' within exon, averaged across exons)
#   Panel 3: Intron metagene (5' to 3' within intron)
#
# No ChIP stratification in this pipeline.
##################################################################

##--------------------------------------------------------------## 0. Libraries

suppressPackageStartupMessages({
  library(rtracklayer)
  library(plyranges)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(AnnotationDbi)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(BiocParallel)
  library(patchwork)
})

reduce <- GenomicRanges::reduce

##--------------------------------------------------------------## 1. Parameters

GTF_PATH  <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/gtf/Homo_sapiens.GRCh38.115.gtf"
EXPR_FILE <- "/projects/b1042/BackmanLab/Lucas/eric_splicing/expr/HCT116_RNAseq_quantification.tsv"
BW_DIR    <- "/projects/b1042/BackmanLab/Tiffany/HP1a/lucas_analysis/coverage/log2ratio"

BW_FILES <- c(
  Auxin         = file.path(BW_DIR, "Auxin_vs_Control.log2ratio.bigWig"),
  Hypoxia       = file.path(BW_DIR, "Hypoxia_vs_Control.log2ratio.bigWig"),
  Auxin_Hypoxia = file.path(BW_DIR, "Auxin_Hypoxia_vs_Control.log2ratio.bigWig")
)

COND_LEVELS <- c("Auxin", "Hypoxia", "Auxin_Hypoxia")
COLORS <- c(
  Auxin         = "#0072B2",
  Hypoxia       = "#E69F00",
  Auxin_Hypoxia = "#CC79A7"
)
LABELS <- c(
  Auxin         = "Auxin",
  Hypoxia       = "Hypoxia",
  Auxin_Hypoxia = "Auxin + Hypoxia"
)

TPM_THRESHOLD <- 1
STD_CHROMS    <- paste0("chr", c(1:22, "X", "Y"))
N_WORKERS     <- 7

## Scale-regions body bins
N_BINS_UP     <- 10
N_BINS_BODY   <- 50
N_BINS_DN     <- 10
FLANK         <- 2000
MIN_BODY      <- N_BINS_BODY + 1L

## Exon and intron bins
N_BINS_EXON   <- 50
N_BINS_INTRON <- 100
MIN_EXON      <- N_BINS_EXON   + 1L
MIN_INTRON    <- N_BINS_INTRON + 1L

## Stratification levels
EXPR_LEVELS <- c("Low/Unexpressed", "Medium", "High")
EI_LEVELS   <- c("Low E/I", "Medium E/I", "High E/I")

##--------------------------------------------------------------## 2. Validate paths

message("Checking paths...")
for (f in c(GTF_PATH, EXPR_FILE, unname(BW_FILES))) {
  if (!file.exists(f)) stop("File not found: ", f)
}
message("  All files found.")

##--------------------------------------------------------------## 3. Build TxDb

message("Building TxDb...")
txdb     <- makeTxDbFromGFF(GTF_PATH, format = "gtf")
genes_gr <- genes(txdb)
genes_gr$gene_id <- sub("\\..*", "", names(genes_gr))
seqlevelsStyle(genes_gr) <- "UCSC"
genes_gr <- keepSeqlevels(genes_gr,
                          intersect(seqlevels(genes_gr), STD_CHROMS),
                          pruning.mode = "coarse")

## Exons by gene
exons_by_gene <- exonsBy(txdb, by = "gene")
names(exons_by_gene) <- sub("\\..*", "", names(exons_by_gene))
seqlevelsStyle(exons_by_gene) <- "UCSC"
exons_by_gene <- keepSeqlevels(exons_by_gene,
                               intersect(seqlevels(exons_by_gene), STD_CHROMS),
                               pruning.mode = "coarse")

## Introns by gene
tx2gene        <- AnnotationDbi::select(txdb,
                                        keys    = keys(txdb, "TXNAME"),
                                        columns = c("TXNAME", "GENEID"),
                                        keytype = "TXNAME")
introns_by_tx  <- intronsByTranscript(txdb, use.names = TRUE)
introns_unlist <- unlist(introns_by_tx, use.names = TRUE)
intron_gids    <- tx2gene$GENEID[match(names(introns_unlist), tx2gene$TXNAME)]
introns_by_gene <- splitAsList(introns_unlist,
                               sub("\\..*", "", intron_gids))
seqlevelsStyle(introns_by_gene) <- "UCSC"
introns_by_gene <- keepSeqlevels(introns_by_gene,
                                 intersect(seqlevels(introns_by_gene),
                                           STD_CHROMS),
                                 pruning.mode = "coarse")
rm(introns_by_tx, introns_unlist, intron_gids); gc()
message("  TxDb built: ", length(genes_gr), " genes")

##--------------------------------------------------------------## 4. Expression annotation

message("Loading TPM data...")
expr <- read_tsv(EXPR_FILE, show_col_types = FALSE) %>%
  dplyr::mutate(
    gene_id = sub("\\..*", "",
                  sapply(strsplit(target_id, "\\|"), `[`, 2))
  ) %>%
  dplyr::select(gene_id, tpm) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(TPM = sum(tpm, na.rm = TRUE), .groups = "drop")

genes_meta <- dplyr::tibble(gene_id = genes_gr$gene_id) %>%
  dplyr::left_join(expr, by = "gene_id") %>%
  dplyr::mutate(TPM = tidyr::replace_na(TPM, 0))

tpm_expressed <- genes_meta$TPM[genes_meta$TPM >= TPM_THRESHOLD]
genes_meta <- genes_meta %>%
  dplyr::mutate(
    ## Merge Unexpressed and Low into one category
    expr_group = dplyr::case_when(
      TPM < quantile(tpm_expressed, 2/3, na.rm = TRUE) ~ "Low/Unexpressed",
      TPM < quantile(tpm_expressed, 2/3, na.rm = TRUE) ~ "Medium",
      TRUE                                              ~ "High"
    ),
    ## Correct case_when — unexpressed + low together, then medium, then high
    expr_group = dplyr::case_when(
      TPM < quantile(tpm_expressed, 1/3, na.rm = TRUE) ~ "Low/Unexpressed",
      TPM < quantile(tpm_expressed, 2/3, na.rm = TRUE) ~ "Medium",
      TRUE                                              ~ "High"
    ),
    expr_group = factor(expr_group, levels = EXPR_LEVELS)
  )

cat("Expression group counts:\n")
print(table(genes_meta$expr_group, useNA = "always"))

##--------------------------------------------------------------## 5. E/I ratio

message("Computing E/I ratio...")
exon_bp <- sum(width(reduce(exons_by_gene)))

gene_lengths <- dplyr::tibble(
  gene_id     = sub("\\..*", "", names(genes_gr)),
  gene_length = as.numeric(width(genes_gr))
)

exon_bp_df <- dplyr::tibble(
  gene_id = sub("\\..*", "", names(exon_bp)),
  exon_bp = as.numeric(exon_bp)
) %>%
  dplyr::left_join(gene_lengths, by = "gene_id") %>%
  dplyr::mutate(
    intron_bp = gene_length - exon_bp,
    ei_ratio  = dplyr::if_else(intron_bp > 0,
                               exon_bp / intron_bp,
                               NA_real_),
    ei_group  = dplyr::case_when(
      is.na(ei_ratio)                                    ~ "High E/I",
      ei_ratio < quantile(ei_ratio, 0.33, na.rm = TRUE) ~ "Low E/I",
      ei_ratio < quantile(ei_ratio, 0.67, na.rm = TRUE) ~ "Medium E/I",
      TRUE                                               ~ "High E/I"
    ),
    ei_group = factor(ei_group, levels = EI_LEVELS)
  )

message("E/I quantile breaks:")
message(paste(
  c("0%", "33%", "67%", "100%"),
  round(quantile(exon_bp_df$ei_ratio, c(0, 0.33, 0.67, 1),
                 na.rm = TRUE), 3),
  sep = ": ", collapse = " | "
))
cat("E/I group counts:\n")
print(table(exon_bp_df$ei_group, useNA = "always"))

##--------------------------------------------------------------## 6. Attach metadata to genes_gr

genes_meta <- genes_meta %>%
  dplyr::left_join(
    exon_bp_df %>% dplyr::select(gene_id, ei_ratio, ei_group),
    by = "gene_id"
  )

genes_gr$TPM        <- genes_meta$TPM
genes_gr$expr_group <- genes_meta$expr_group
genes_gr$ei_group   <- genes_meta$ei_group

##--------------------------------------------------------------## 7. safe_tile() helper

safe_tile <- function(gr, n_bins, label = "features") {
  tiled       <- tile(gr, n = n_bins)
  tile_counts <- lengths(tiled)
  bad         <- tile_counts != n_bins
  if (any(bad)) {
    message("  Dropping ", sum(bad), " ", label)
    gr    <- gr[!bad]
    tiled <- tiled[!bad]
  }
  list(gr = gr, tiled = tiled)
}

##--------------------------------------------------------------## 8. Bin construction helpers
## All bins carry expr_group and ei_group as character columns.

attach_bins_metadata <- function(flat, gr, n_bins, nf, feature_label) {
  bin_idx     <- rep.int(seq_len(n_bins), nf)
  feature_idx <- rep.int(seq_len(nf), rep.int(n_bins, nf))
  is_minus    <- as.logical(strand(gr) == "-")
  minus_mask  <- rep.int(is_minus, rep.int(n_bins, nf))
  bin_idx[minus_mask] <- (n_bins + 1L) - bin_idx[minus_mask]
  
  mcols(flat) <- DataFrame(
    feature_idx = feature_idx,
    bin_idx     = bin_idx,
    feature     = feature_label,
    gene_id     = rep.int(gr$gene_id,
                          rep.int(n_bins, nf)),
    expr_group  = rep.int(as.character(gr$expr_group),
                          rep.int(n_bins, nf)),
    ei_group    = rep.int(as.character(gr$ei_group),
                          rep.int(n_bins, nf))
  )
  flat
}

## TSS-TES scale-regions body bins
make_body_bins <- function(genes_gr) {
  feat_gr <- genes_gr[width(genes_gr) > MIN_BODY]
  feat_gr <- keepSeqlevels(feat_gr,
                           intersect(seqlevels(feat_gr), STD_CHROMS),
                           pruning.mode = "coarse")
  res     <- safe_tile(feat_gr, N_BINS_BODY, "body bins")
  feat_gr <- res$gr; tiled <- res$tiled
  flat    <- unlist(tiled, use.names = FALSE)
  attach_bins_metadata(flat, feat_gr, N_BINS_BODY, length(feat_gr), "Body")
}

## Upstream flank bins
make_upstream_bins <- function(genes_gr) {
  tss_gr  <- GenomicRanges::resize(genes_gr, width = 1, fix = "start")
  feat_gr <- GenomicRanges::flank(tss_gr, width = FLANK,
                                  start = TRUE, both = FALSE)
  feat_gr <- feat_gr[width(feat_gr) > N_BINS_UP]
  feat_gr <- keepSeqlevels(feat_gr,
                           intersect(seqlevels(feat_gr), STD_CHROMS),
                           pruning.mode = "coarse")
  res     <- safe_tile(feat_gr, N_BINS_UP, "upstream bins")
  feat_gr <- res$gr; tiled <- res$tiled
  flat    <- unlist(tiled, use.names = FALSE)
  attach_bins_metadata(flat, feat_gr, N_BINS_UP, length(feat_gr), "Upstream")
}

## Downstream flank bins
make_downstream_bins <- function(genes_gr) {
  tes_gr  <- GenomicRanges::resize(genes_gr, width = 1, fix = "end")
  feat_gr <- GenomicRanges::flank(tes_gr, width = FLANK,
                                  start = FALSE, both = FALSE)
  feat_gr <- feat_gr[width(feat_gr) > N_BINS_DN]
  feat_gr <- keepSeqlevels(feat_gr,
                           intersect(seqlevels(feat_gr), STD_CHROMS),
                           pruning.mode = "coarse")
  res     <- safe_tile(feat_gr, N_BINS_DN, "downstream bins")
  feat_gr <- res$gr; tiled <- res$tiled
  flat    <- unlist(tiled, use.names = FALSE)
  attach_bins_metadata(flat, feat_gr, N_BINS_DN, length(feat_gr), "Downstream")
}

## Exon bins — each reduced exon tiled independently
make_exon_bins <- function(exons_by_gene, genes_gr) {
  exons_red  <- GenomicRanges::reduce(exons_by_gene)
  flat_exons <- unlist(exons_red, use.names = TRUE)
  flat_exons$gene_id <- sub("\\..*", "", names(flat_exons))
  flat_exons <- flat_exons[width(flat_exons) > MIN_EXON]
  seqlevelsStyle(flat_exons) <- "UCSC"
  flat_exons <- keepSeqlevels(flat_exons,
                              intersect(seqlevels(flat_exons), STD_CHROMS),
                              pruning.mode = "coarse")
  
  expr_map <- setNames(as.character(genes_gr$expr_group), genes_gr$gene_id)
  ei_map   <- setNames(as.character(genes_gr$ei_group),   genes_gr$gene_id)
  flat_exons$expr_group <- expr_map[flat_exons$gene_id]
  flat_exons$ei_group   <- ei_map[flat_exons$gene_id]
  
  res        <- safe_tile(flat_exons, N_BINS_EXON, "exon bins")
  flat_exons <- res$gr; tiled <- res$tiled
  flat       <- unlist(tiled, use.names = FALSE)
  attach_bins_metadata(flat, flat_exons, N_BINS_EXON,
                       length(flat_exons), "Exon")
}

## Intron bins
make_intron_bins <- function(introns_by_gene, genes_gr) {
  introns_red  <- GenomicRanges::reduce(introns_by_gene)
  flat_introns <- unlist(introns_red, use.names = TRUE)
  flat_introns$gene_id <- sub("\\..*", "", names(flat_introns))
  flat_introns <- flat_introns[width(flat_introns) > MIN_INTRON]
  seqlevelsStyle(flat_introns) <- "UCSC"
  flat_introns <- keepSeqlevels(flat_introns,
                                intersect(seqlevels(flat_introns), STD_CHROMS),
                                pruning.mode = "coarse")
  
  expr_map <- setNames(as.character(genes_gr$expr_group), genes_gr$gene_id)
  ei_map   <- setNames(as.character(genes_gr$ei_group),   genes_gr$gene_id)
  flat_introns$expr_group <- expr_map[flat_introns$gene_id]
  flat_introns$ei_group   <- ei_map[flat_introns$gene_id]
  
  res          <- safe_tile(flat_introns, N_BINS_INTRON, "intron bins")
  flat_introns <- res$gr; tiled <- res$tiled
  flat         <- unlist(tiled, use.names = FALSE)
  attach_bins_metadata(flat, flat_introns, N_BINS_INTRON,
                       length(flat_introns), "Intron")
}

##--------------------------------------------------------------## 9. Signal extraction

extract_one_chrom <- function(chrom, bw_gr, bins) {
  bc <- bw_gr[as.character(seqnames(bw_gr)) == chrom]
  bn <- bins[as.character(seqnames(bins))   == chrom]
  if (!length(bc) || !length(bn)) return(NULL)
  tryCatch(
    plyranges::join_overlap_inner(bn, bc) %>%
      dplyr::as_tibble() %>%
      dplyr::rename(ow = width) %>%
      dplyr::group_by(gene_id, feature, feature_idx, bin_idx,
                      expr_group, ei_group) %>%
      dplyr::summarise(
        mean_signal = sum(score * ow) / sum(ow),
        .groups     = "drop"
      ),
    error = function(e) NULL
  )
}

extract_signal <- function(bw_gr, bins) {
  chroms <- intersect(as.character(unique(seqnames(bins))),
                      as.character(unique(seqnames(bw_gr))))
  BiocParallel::register(BiocParallel::MulticoreParam(workers = N_WORKERS))
  dplyr::bind_rows(
    BiocParallel::bplapply(chroms, extract_one_chrom,
                           bw_gr = bw_gr, bins = bins)
  )
}

##--------------------------------------------------------------## 10. Profile summarisation
## Body/upstream/downstream: x_pos from bin position + offset
## Exon/intron: x_pos = bin_idx (no offset)
## Exon summarisation: two-step average
##   Step 1: average across exon instances per gene per bin_idx
##   Step 2: average across genes

add_x_pos_body <- function(df) {
  df %>% dplyr::mutate(
    x_pos = dplyr::case_when(
      feature == "Upstream"   ~ as.numeric(bin_idx),
      feature == "Body"       ~ as.numeric(bin_idx) + N_BINS_UP,
      feature == "Downstream" ~ as.numeric(bin_idx) + N_BINS_UP + N_BINS_BODY,
      TRUE                    ~ as.numeric(bin_idx)
    )
  )
}

summarise_body <- function(df, sample_name) {
  df %>%
    dplyr::mutate(condition = sample_name) %>%
    dplyr::group_by(condition, feature, expr_group, ei_group, x_pos) %>%
    dplyr::summarise(
      mean    = mean(mean_signal, na.rm = TRUE),
      se      = sd(mean_signal,   na.rm = TRUE) /
        sqrt(sum(!is.na(mean_signal))),
      n_genes = sum(!is.na(mean_signal)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(upper = mean + 1.96*se, lower = mean - 1.96*se)
}

summarise_exon <- function(df, sample_name) {
  ## Two-step: collapse exon instances within gene first
  df %>%
    dplyr::mutate(condition = sample_name,
                  x_pos     = as.numeric(bin_idx)) %>%
    dplyr::group_by(condition, feature, expr_group,
                    ei_group, gene_id, x_pos) %>%
    dplyr::summarise(
      gene_mean = mean(mean_signal, na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    dplyr::group_by(condition, feature, expr_group, ei_group, x_pos) %>%
    dplyr::summarise(
      mean    = mean(gene_mean, na.rm = TRUE),
      se      = sd(gene_mean,   na.rm = TRUE) /
        sqrt(sum(!is.na(gene_mean))),
      n_genes = sum(!is.na(gene_mean)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(upper = mean + 1.96*se, lower = mean - 1.96*se)
}

summarise_intron <- function(df, sample_name) {
  df %>%
    dplyr::mutate(condition = sample_name,
                  x_pos     = as.numeric(bin_idx)) %>%
    dplyr::group_by(condition, feature, expr_group, ei_group, x_pos) %>%
    dplyr::summarise(
      mean    = mean(mean_signal, na.rm = TRUE),
      se      = sd(mean_signal,   na.rm = TRUE) /
        sqrt(sum(!is.na(mean_signal))),
      n_genes = sum(!is.na(mean_signal)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(upper = mean + 1.96*se, lower = mean - 1.96*se)
}

##--------------------------------------------------------------## 11. Plot theme and helpers

N_TOT    <- N_BINS_UP + N_BINS_BODY + N_BINS_DN
x_breaks_body <- c(1, N_BINS_UP + 0.5,
                   N_BINS_UP + N_BINS_BODY + 0.5, N_TOT)
x_labels_body <- c(paste0("-", FLANK/1000, "kb"), "TSS",
                   "TES", paste0("+", FLANK/1000, "kb"))
boundary_lines <- data.frame(
  xintercept = c(N_BINS_UP + 0.5, N_BINS_UP + N_BINS_BODY + 0.5)
)
x_breaks_sub <- c(1, N_BINS_EXON/2,   N_BINS_EXON)
x_labels_sub <- c("5'", "50%", "3'")
x_breaks_int <- c(1, N_BINS_INTRON/2, N_BINS_INTRON)
x_labels_int <- c("5'", "50%", "3'")

hp1a_theme <- theme_classic(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 11),
    legend.position  = "top",
    legend.title     = element_blank(),
    legend.key.size  = unit(0.8, "lines"),
    axis.line        = element_line(color = "grey40"),
    axis.ticks       = element_line(color = "grey40"),
    axis.text        = element_text(color = "grey20"),
    axis.title       = element_text(color = "grey20"),
    panel.spacing    = unit(1.1, "lines"),
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 11)
  )

make_n_labels <- function(gdf, facet_col = NULL) {
  if (is.null(facet_col)) {
    dplyr::tibble(
      label = paste0("n = ", nrow(gdf))
    )
  } else {
    gdf %>%
      dplyr::group_by(.data[[facet_col]]) %>%
      dplyr::summarise(n_genes = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(label = paste0("n = ", n_genes)) %>%
      dplyr::rename(facet_var = 1)
  }
}

## Build a single body/exon/intron panel ggplot
make_panel <- function(df, x_breaks, x_labels, x_lab,
                       y_lab, title_str, subtitle_str,
                       facet_col = NULL, label_df = NULL,
                       show_boundary = FALSE) {
  df <- df %>%
    dplyr::mutate(condition = factor(condition, levels = COND_LEVELS))
  
  if (!is.null(facet_col))
    df <- df %>% dplyr::mutate(facet_var = .data[[facet_col]])
  
  p <- ggplot(df, aes(x = x_pos, y = mean,
                      color = condition, fill = condition)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.10, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dotted",
               color = "grey60", linewidth = 0.4) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_color_manual(values = COLORS, labels = LABELS) +
    scale_fill_manual(values  = COLORS, labels = LABELS) +
    labs(x = x_lab, y = y_lab,
         title = title_str, subtitle = subtitle_str) +
    hp1a_theme
  
  if (show_boundary)
    p <- p + geom_vline(data = boundary_lines,
                        aes(xintercept = xintercept),
                        linetype = "dashed", color = "grey50",
                        linewidth = 0.4, inherit.aes = FALSE)
  
  if (!is.null(label_df)) {
    if (!is.null(facet_col)) {
      p <- p + geom_text(data = label_df,
                         aes(x = Inf, y = Inf, label = label),
                         inherit.aes = FALSE,
                         hjust = 1.1, vjust = 1.5,
                         size = 3.0, color = "grey40", fontface = "italic")
    } else {
      p <- p + annotate("text", x = Inf, y = Inf, label = label_df$label,
                        hjust = 1.1, vjust = 1.5,
                        size = 3.0, color = "grey40", fontface = "italic")
    }
  }
  
  if (!is.null(facet_col))
    p <- p + facet_wrap(~ facet_var, nrow = 1, scales = "free_y")
  
  p
}

## Collapse strat columns for a given facet variable (or none)
collapse <- function(df, facet_col = NULL, feature_filter) {
  df <- df %>% dplyr::filter(feature == feature_filter)
  group_vars <- c("condition", "feature", "x_pos")
  if (!is.null(facet_col)) group_vars <- c(group_vars, facet_col)
  df %>%
    dplyr::group_by(across(all_of(group_vars))) %>%
    dplyr::summarise(
      mean = mean(mean, na.rm = TRUE),
      se   = mean(se,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(upper = mean + 1.96*se, lower = mean - 1.96*se)
}

## Assemble three-panel figure
assemble_figure <- function(body_df, exon_df, intron_df,
                            gdf, facet_col = NULL,
                            title_str, subtitle_str) {
  lev <- if (!is.null(facet_col)) {
    if (facet_col == "expr_group") EXPR_LEVELS else EI_LEVELS
  } else NULL
  
  if (!is.null(facet_col)) {
    body_df   <- body_df   %>%
      dplyr::mutate(facet_var = factor(.data[[facet_col]], levels = lev))
    exon_df   <- exon_df   %>%
      dplyr::mutate(facet_var = factor(.data[[facet_col]], levels = lev))
    intron_df <- intron_df %>%
      dplyr::mutate(facet_var = factor(.data[[facet_col]], levels = lev))
    label_df  <- make_n_labels(gdf, facet_col) %>%
      dplyr::mutate(facet_var = factor(facet_var, levels = lev))
  } else {
    label_df <- make_n_labels(gdf)
  }
  
  p_body <- make_panel(
    body_df,
    x_breaks = x_breaks_body, x_labels = x_labels_body,
    x_lab = NULL, y_lab = "log2(COND/CTRL)",
    title_str    = "Gene body (TSS → TES)",
    subtitle_str = paste0("-", FLANK/1000, "kb  |  scaled gene body  |  +",
                          FLANK/1000, "kb"),
    facet_col      = if (!is.null(facet_col)) "facet_var" else NULL,
    label_df       = label_df,
    show_boundary  = TRUE
  )
  
  p_exon <- make_panel(
    exon_df,
    x_breaks = x_breaks_sub, x_labels = x_labels_sub,
    x_lab = "Position within exon (5' to 3')",
    y_lab = "log2(COND/CTRL)",
    title_str    = "Exon signal",
    subtitle_str = "Averaged across all exons per gene",
    facet_col    = if (!is.null(facet_col)) "facet_var" else NULL,
    label_df     = NULL
  )
  
  p_intron <- make_panel(
    intron_df,
    x_breaks = x_breaks_int, x_labels = x_labels_int,
    x_lab = "Position within intron (5' to 3')",
    y_lab = "log2(COND/CTRL)",
    title_str    = "Intron signal",
    subtitle_str = "Averaged across all introns per gene",
    facet_col    = if (!is.null(facet_col)) "facet_var" else NULL,
    label_df     = NULL
  )
  
  p_body / p_exon / p_intron +
    plot_layout(heights = c(2, 1.5, 1.5), guides = "collect") +
    plot_annotation(
      title    = title_str,
      subtitle = subtitle_str,
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(color = "grey40", size = 11)
      )
    ) & theme(legend.position = "top")
}

##--------------------------------------------------------------## 12. Build bins and extract signal

message("Building bins...")
up_bins     <- make_upstream_bins(genes_gr)
body_bins   <- make_body_bins(genes_gr)
dn_bins     <- make_downstream_bins(genes_gr)
exon_bins   <- make_exon_bins(exons_by_gene, genes_gr)
intron_bins <- make_intron_bins(introns_by_gene, genes_gr)

meta_bins <- c(up_bins, body_bins, dn_bins)
rm(up_bins, body_bins, dn_bins); gc()

message("Extracting signal for all conditions...")
profile_list <- lapply(names(BW_FILES), function(cond_name) {
  cat("  Processing:", cond_name, "\n")
  bw_gr <- plyranges::read_bigwig(BW_FILES[cond_name])
  seqlevelsStyle(bw_gr) <- "UCSC"
  bw_gr <- keepSeqlevels(bw_gr,
                         intersect(seqlevels(bw_gr), STD_CHROMS),
                         pruning.mode = "coarse")
  
  body_raw   <- extract_signal(bw_gr, meta_bins)   %>% add_x_pos_body()
  exon_raw   <- extract_signal(bw_gr, exon_bins)
  intron_raw <- extract_signal(bw_gr, intron_bins)
  rm(bw_gr); gc()
  
  list(
    body   = summarise_body(body_raw,   cond_name),
    exon   = summarise_exon(exon_raw,   cond_name),
    intron = summarise_intron(intron_raw, cond_name)
  )
})

rm(meta_bins, exon_bins, intron_bins); gc()

body_df   <- dplyr::bind_rows(lapply(profile_list, `[[`, "body")) %>%
  dplyr::mutate(condition = factor(condition, levels = COND_LEVELS))
exon_df   <- dplyr::bind_rows(lapply(profile_list, `[[`, "exon")) %>%
  dplyr::mutate(condition = factor(condition, levels = COND_LEVELS))
intron_df <- dplyr::bind_rows(lapply(profile_list, `[[`, "intron")) %>%
  dplyr::mutate(condition = factor(condition, levels = COND_LEVELS))

rm(profile_list); gc()

## genes metadata tibble for labels
gdf <- dplyr::as_tibble(mcols(genes_gr)) %>%
  dplyr::mutate(
    expr_group = factor(expr_group, levels = EXPR_LEVELS),
    ei_group   = factor(ei_group,   levels = EI_LEVELS)
  )

SUBTITLE_BASE <- "log2(Condition/Control)  |  HP1a  |  HCT116"

##--------------------------------------------------------------## 13. Plot 1 — All genes, no faceting

message("Plotting Plot 1 — all genes...")

p1 <- assemble_figure(
  body_df   = collapse(body_df,   NULL, "Body") %>%
    dplyr::bind_rows(collapse(body_df, NULL, "Upstream")) %>%
    dplyr::bind_rows(collapse(body_df, NULL, "Downstream")),
  exon_df   = collapse(exon_df,   NULL, "Exon"),
  intron_df = collapse(intron_df, NULL, "Intron"),
  gdf       = gdf,
  facet_col = NULL,
  title_str    = "Plot 1 — All genes",
  subtitle_str = SUBTITLE_BASE
)
print(p1)

##--------------------------------------------------------------## 14. Plot 2 — Faceted by expression group

message("Plotting Plot 2 — expression group faceting...")

p2 <- assemble_figure(
  body_df   = dplyr::bind_rows(
    collapse(body_df,   "expr_group", "Body"),
    collapse(body_df,   "expr_group", "Upstream"),
    collapse(body_df,   "expr_group", "Downstream")
  ),
  exon_df   = collapse(exon_df,   "expr_group", "Exon"),
  intron_df = collapse(intron_df, "expr_group", "Intron"),
  gdf       = gdf,
  facet_col = "expr_group",
  title_str    = "Plot 2 — Faceted by expression level",
  subtitle_str = SUBTITLE_BASE
)
print(p2)

##--------------------------------------------------------------## 15. Plot 3 — Faceted by E/I group

message("Plotting Plot 3 — E/I group faceting...")

p3 <- assemble_figure(
  body_df   = dplyr::bind_rows(
    collapse(body_df,   "ei_group", "Body"),
    collapse(body_df,   "ei_group", "Upstream"),
    collapse(body_df,   "ei_group", "Downstream")
  ),
  exon_df   = collapse(exon_df,   "ei_group", "Exon"),
  intron_df = collapse(intron_df, "ei_group", "Intron"),
  gdf       = gdf,
  facet_col = "ei_group",
  title_str    = "Plot 3 — Faceted by E/I group",
  subtitle_str = SUBTITLE_BASE
)
print(p3)
