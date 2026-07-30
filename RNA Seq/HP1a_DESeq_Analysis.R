library(tximport)
library(DESeq2)
library(pheatmap)

# --- File paths ---
basedir <- "C:/Users/tkx1376/OneDrive - Northwestern University/RNA Sequencing/HP1a/seq_raw_output"
samples_file <- file.path(basedir, "samples.csv")

# --- Load sample metadata ---
samples <- read.csv(samples_file, header=TRUE, stringsAsFactors=FALSE)

# Set condition factor with Control as baseline
samples$condition <- factor(samples$condition, 
                            levels=c("Control","Auxin","Hypoxia","Auxin_Hypoxia"))

# Use clean sample names
samples$sample_name <- samples$samples
rownames(samples) <- samples$sample_name

# Check for duplicates
if(any(duplicated(samples$sample_name))){
  stop("Duplicate sample names found!")
}

# --- Build file list ---
files <- file.path(basedir, paste0(samples$sample_name, "_R1_001Aligned.genes.results"))
names(files) <- samples$sample_name

# Check files exist
missing_files <- files[!file.exists(files)]
if(length(missing_files) > 0){
  stop("Missing files:\n", paste(missing_files, collapse="\n"))
}

# Check each file has rows and expected_count
valid_files <- sapply(files, function(f){
  df <- read.table(f, header=TRUE)
  all(c("expected_count") %in% colnames(df)) && nrow(df) > 0
})
if(!all(valid_files)){
  stop("Some files are invalid:\n", paste(names(files)[!valid_files], collapse="\n"))
}

# Keep only valid files
files <- files[valid_files]
samples <- samples[names(files), ]

cat("Using", length(files), "valid RSEM files:\n")
print(names(files))

# --- Import RSEM ---
txi.rsem <- tximport(files, type="rsem", txIn=FALSE, txOut=FALSE)

# --- Filter out zero-length or unexpressed genes ---
zero_genes <- (apply(txi.rsem$abundance, 1, max) == 0) | 
  (apply(txi.rsem$length, 1, min) == 0)

cat("Filtering out", sum(zero_genes), "genes with zero expression or zero length\n")

txi_filtered <- list(
  counts = txi.rsem$counts[!zero_genes, ],
  abundance = txi.rsem$abundance[!zero_genes, ],
  length = txi.rsem$length[!zero_genes, ],
  countsFromAbundance = txi.rsem$countsFromAbundance
)

# --- Create DESeq2 dataset ---
ddsTxi <- DESeqDataSetFromTximport(txi_filtered,
                                   colData = samples,
                                   design = ~condition)

# --- Optional: filter very low counts ---
keep <- rowSums(counts(ddsTxi)) >= 10
ddsTxi <- ddsTxi[keep, ]

cat("Number of genes kept for DESeq2:", nrow(ddsTxi), "\n")

# --- Run DESeq2 ---
dds <- DESeq(ddsTxi)

# --- Extract contrasts ---
res_auxin   <- results(dds, contrast=c("condition","Auxin","Control"))
res_hypoxia <- results(dds, contrast=c("condition","Hypoxia","Control"))
res_auxhyp  <- results(dds, contrast=c("condition","Auxin_Hypoxia","Control"))

# Variance-stabilizing transformation
vsd <- vst(dds, blind=FALSE)

#PCA plot
plotPCA(vsd, intgroup="condition")

# Plot Pearson's correlation

# Extract transformed expression matrix
expr_matrix <- assay(vsd)

files <- c(
  Control_1       = file.path(basedir, "Control_1_S25_L005_R1_001Aligned.genes.results"),
  Control_2       = file.path(basedir, "Control_2_S26_L005_R1_001Aligned.genes.results"),
  Control_3       = file.path(basedir, "Control_3_S27_L005_R1_001Aligned.genes.results"),
  Auxin_1         = file.path(basedir, "Auxin_1_S28_L005_R1_001Aligned.genes.results"),
  Auxin_2         = file.path(basedir, "Auxin_2_S29_L005_R1_001Aligned.genes.results"),
  Auxin_3         = file.path(basedir, "Auxin_3_S30_L005_R1_001Aligned.genes.results"),
  Hypoxia_1       = file.path(basedir, "Hypoxia_1_S31_L005_R1_001Aligned.genes.results"),
  Hypoxia_2       = file.path(basedir, "Hypoxia_2_S32_L005_R1_001Aligned.genes.results"),
  Hypoxia_3       = file.path(basedir, "Hypoxia_3_S33_L005_R1_001Aligned.genes.results"),
  Auxin_Hypoxia_1 = file.path(basedir, "Auxin_Hypoxia_1_S34_L005_R1_001Aligned.genes.results"),
  Auxin_Hypoxia_2 = file.path(basedir, "Auxin_Hypoxia_2_S35_L005_R1_001Aligned.genes.results"),
  Auxin_Hypoxia_3 = file.path(basedir, "Auxin_Hypoxia_3_S36_L005_R1_001Aligned.genes.results")
)
CONDITION_MAP <- c(
  Control_1_S25_L005 = "Control", Control_2_S26_L005 = "Control", Control_3_S27_L005 = "Control",
  Auxin_1_S28_L005   = "Auxin",   Auxin_2_S29_L005   = "Auxin",  Auxin_3_S30_L005   = "Auxin",
  Hypoxia_1_S31_L005 = "Hypoxia", Hypoxia_2_S32_L005 = "Hypoxia", Hypoxia_3_S33_L005 = "Hypoxia",
  Auxin_Hypoxia_1_S34_L005 = "Auxin_Hypoxia",
  Auxin_Hypoxia_2_S35_L005 = "Auxin_Hypoxia",
  Auxin_Hypoxia_3_S36_L005 = "Auxin_Hypoxia"
)

# -----------------------------
# Compute sample correlations
# -----------------------------
sample_cor <- cor(expr_matrix)

# View correlation matrix
print(sample_cor)

## Annotation bar showing condition per sample
annotation_col <- data.frame(
  Condition = CONDITION_MAP[colnames(sample_cor)],
  row.names = colnames(sample_cor)
)

library(RColorBrewer)

CONDITION_COLORS <- c(
  Control       = "#999999",
  Auxin         = "#0072B2",
  Hypoxia       = "#E69F00",
  Auxin_Hypoxia = "#CC79A7"
)

annotation_colors <- list(Condition = CONDITION_COLORS)

## Color scale — range slightly below min correlation to 1
cor_min  <- floor(min(sample_cor) * 100) / 100
breaks   <- seq(cor_min, 1, length.out = 101)
colors   <- colorRampPalette(
  rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu"))
)(100)

clean_names <- gsub("_S[0-9]+_L[0-9]+", "", colnames(sample_cor))

pheatmap(
  sample_cor,
  color            = colors,
  breaks           = breaks,
  annotation_col   = annotation_col,
  annotation_row   = annotation_col,
  annotation_colors = annotation_colors,
  cluster_rows     = TRUE,
  cluster_cols     = TRUE,
  labels_col = clean_names,
  labels_row = clean_names,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "complete",
  show_rownames    = TRUE,
  show_colnames    = TRUE,
  fontsize         = 8,
  fontsize_row     = 10,
  fontsize_col     = 10,
  display_numbers  = TRUE,
  number_format    = "%.3f",
  number_color     = "grey20",
  fontsize_number  = 8,
  main             =  "Sample Correlation Plot",
  border_color     = NA,
  cellwidth        = 25,
  cellheight       = 25
)

# -----------------------------
# Plot correlation heatmap
# -----------------------------
pheatmap(
  sample_cor,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Sample Correlation Heatmap"
)

# --- Save results ---
write.csv(as.data.frame(res_auxin),
          file=file.path(basedir,"DESeq2_Auxin_vs_Control.csv"))
write.csv(as.data.frame(res_hypoxia),
          file=file.path(basedir,"DESeq2_Hypoxia_vs_Control.csv"))
write.csv(as.data.frame(res_auxhyp),
          file=file.path(basedir,"DESeq2_AuxinHypoxia_vs_Control.csv"))
