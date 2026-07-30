# =============================================

# AUXIN-ONLY GO ENRICHMENT + BUBBLE PLOTS

# =============================================

# ======================

# 0️⃣ Load required packages

# ======================

if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "biomaRt", "AnnotationDbi", "GO.db"), update=FALSE, ask=FALSE)
install.packages(c("dplyr", "ggplot2"))

library(biomaRt)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)

# ======================

# 1️⃣ Define directories & file

# ======================

base_dir <- "C:/Users/tkx1376/OneDrive - Northwestern University/RNA Sequencing/HP1a/rsem"
output_dir <- file.path(base_dir, "GO_BubblePlots")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

res_file <- file.path(base_dir, "DESeq2_Auxin_vs_Control.csv")

# ======================

# 2️⃣ Map Ensembl → Entrez

# ======================

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

map_ensembl_to_entrez <- function(ensembl_ids){
  mapping <- getBM(
    attributes = c("ensembl_gene_id","entrezgene_id"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl
  )
  mapping <- mapping[!is.na(mapping$entrezgene_id), ]
  return(mapping)
}

add_entrez <- function(df){
  mapping <- map_ensembl_to_entrez(rownames(df))
  df$EntrezID <- mapping$entrezgene_id[match(rownames(df), mapping$ensembl_gene_id)]
  df <- df[!is.na(df$EntrezID), ]
  return(df)
}

# ======================

# 3️⃣ Load DESeq2 + split up/down

# ======================

df <- read.csv(res_file, row.names = 1)
df <- add_entrez(df)

up_genes <- df %>% filter(log2FoldChange > 0) %>% pull(EntrezID)
down_genes <- df %>% filter(log2FoldChange < 0) %>% pull(EntrezID)

# ======================

# 4️⃣ Run GO enrichment

# ======================

run_enrich <- function(entrez_ids){
  if(length(entrez_ids) == 0) return(NULL)
  enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
}

ego_up <- run_enrich(up_genes)
ego_down <- run_enrich(down_genes)

# ======================

# 5️⃣ Prep for plotting

# ======================

prep_for_plot <- function(go_obj, regulation){
  if(is.null(go_obj)) return(NULL)
  go_df <- as.data.frame(go_obj)
  
  go_df %>%
    mutate(
      GeneRatio = sapply(strsplit(as.character(GeneRatio), "/"),
                         function(x) as.numeric(x[1])/as.numeric(x[2])),
      adj_pval = p.adjust,
      Regulation = regulation
    ) %>%
    dplyr::select(Description, GeneRatio, adj_pval, Regulation)
}

plot_df <- bind_rows(
  prep_for_plot(ego_up, "Upregulated"),
  prep_for_plot(ego_down, "Downregulated")
)

# ======================

# 6️⃣ Select top GO terms

# ======================

top_terms <- function(df, regulation){
  df %>%
    filter(Regulation == regulation) %>%
    arrange(adj_pval) %>%
    pull(Description) %>%
    unique() %>%
    head(15)
}

top_up_terms <- top_terms(plot_df, "Upregulated")
top_down_terms <- top_terms(plot_df, "Downregulated")

plot_up <- plot_df %>%
  filter(Regulation == "Upregulated", Description %in% top_up_terms)

plot_down <- plot_df %>%
  filter(Regulation == "Downregulated", Description %in% top_down_terms)

# Factor ordering

plot_up$Description <- factor(plot_up$Description, levels = rev(top_up_terms))
plot_down$Description <- factor(plot_down$Description, levels = rev(top_down_terms))

# ======================

# 7️⃣ Plot: Upregulated

# ======================

p_up <- ggplot(plot_up, aes(x = "Auxin", y = Description)) +
  geom_point(aes(size = GeneRatio, color = adj_pval), alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red", name = "Adjusted p-value") +
  scale_size(range = c(3, 10), name = "Gene Ratio") +
  guides(
    size = guide_legend(order = 1),      # 👈 TOP
    color = guide_colorbar(order = 2)    # 👈 BOTTOM
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.6)
  ) +
  labs(title = "GO Enrichment Analysis of Upregulated Genes in Auxin", x = NULL, y = NULL)

# ======================

# 8️⃣ Plot: Downregulated

# ======================

p_down <- ggplot(plot_down, aes(x = "Auxin", y = Description)) +
  geom_point(aes(size = GeneRatio, color = adj_pval), alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red", name = "Adjusted p-value") +
  scale_size(range = c(3, 10), name = "Gene Ratio") +
  guides(
    size = guide_legend(order = 1),      # 👈 TOP
    color = guide_colorbar(order = 2)    # 👈 BOTTOM
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.6) 
  )+
  labs(title = "GO Enrichment Analysis of Downregulated Genes in Auxin", x = NULL, y = NULL)

# ======================

# 💾 Save plots

# ======================

ggsave(file.path(output_dir, "GO Analysis_Auxin_Upregulated.png"), p_up, width = 10, height = 7)
ggsave(file.path(output_dir, "GO_Analysis_Auxin_Downregulated.png"), p_down, width = 10, height = 7)

# ======================

# 👀 Display

# ======================

print(p_up)
print(p_down)

# PLOT HEATMAP
# =============================================
# HEATMAP: CONTROL vs AUXIN (65 DEGs)
# =============================================

library(dplyr)
library(pheatmap)
library(readr)
library(purrr)
library(tibble)
library(grid)

# =============================================================
# 0️⃣ Load DESeq2 Results
# =============================================================

res_auxin <- read.csv("C:/Users/tkx1376/OneDrive - Northwestern University/RNA Sequencing/HP1a/rsem/DESeq2_Auxin_vs_Control.csv")
colnames(res_auxin)[1] <- "gene_id"

# =============================================================
# 1️⃣ Get significant genes (padj < 0.05)
# =============================================================

sig_genes <- res_auxin %>%
  filter(!is.na(padj), padj < 0.05) %>%
  pull(gene_id)

# =============================================================
# 2️⃣ TPM files
# =============================================================

files_ctrl <- c(
  "C:/Users/tkx1376/OneDrive - Northwestern University/RNA Sequencing/HP1a/seq_raw_output/Control_1_S25_L005_R1_001Aligned.genes.results",
  "C:/Users/tkx1376/OneDrive - Northwestern University/RNA Sequencing/HP1a/seq_raw_output/Control_2_S26_L005_R1_001Aligned.genes.results",
  "C:/Users/tkx1376/OneDrive - Northwestern University/RNA Sequencing/HP1a/seq_raw_output/Control_3_S27_L005_R1_001Aligned.genes.results"
)

files_auxin <- c(
  "C:/Users/tkx1376/OneDrive - Northwestern University/RNA Sequencing/HP1a/seq_raw_output/Auxin_1_S28_L005_R1_001Aligned.genes.results",
  "C:/Users/tkx1376/OneDrive - Northwestern University/RNA Sequencing/HP1a/seq_raw_output/Auxin_2_S29_L005_R1_001Aligned.genes.results",
  "C:/Users/tkx1376/OneDrive - Northwestern University/RNA Sequencing/HP1a/seq_raw_output/Auxin_3_S30_L005_R1_001Aligned.genes.results"
)

# =============================================================
# 3️⃣ Read replicates
# =============================================================

read_replicates <- function(files, label){
  mat_list <- lapply(seq_along(files), function(i){
    df <- read_tsv(files[i], show_col_types = FALSE)
    tibble(
      gene_id = df$gene_id,
      !!paste0(label, " ", i) := log2(df$TPM + 1)
    )
  })
  reduce(mat_list, full_join, by = "gene_id")
}

ctrl_mat  <- read_replicates(files_ctrl, "Control")
auxin_mat <- read_replicates(files_auxin, "Auxin")

# =============================================================
# 4️⃣ Combine + filter
# =============================================================

expr_mat <- reduce(list(ctrl_mat, auxin_mat), full_join, by = "gene_id")

expr_mat <- expr_mat %>%
  filter(gene_id %in% sig_genes)

gene_ids <- expr_mat$gene_id
expr_mat <- expr_mat %>% select(-gene_id)
expr_mat <- as.matrix(expr_mat)
rownames(expr_mat) <- gene_ids

# =============================================================
# 5️⃣ Select top 65 most variable genes
# =============================================================

gene_variance <- apply(expr_mat, 1, var)

top_n <- min(65, nrow(expr_mat))
top_idx <- order(gene_variance, decreasing = TRUE)[1:top_n]

expr_mat <- expr_mat[top_idx, ]

# =============================================================
# 🧾 Output directory
# =============================================================

output_dir <- "C:/Users/tkx1376/OneDrive - Northwestern University/RNA Sequencing/HP1a/Final Heatmaps"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# =============================================================
# 🧾 Export Ensembl IDs
# =============================================================

write.csv(
  data.frame(Ensembl_ID = rownames(expr_mat)),
  file = file.path(output_dir, "65_DEG_Ensembl_IDs.csv"),
  row.names = FALSE
)

# =============================================================
# 🔥 Manual gene mapping (FIXED VERSION)
# =============================================================

gene_map <- c(
  "ENSG00000066827"="ZFAT",
  "ENSG00000067715"="SYT1",
  "ENSG00000089177"="KIF16B",
  "ENSG00000107077"="KDM4C",
  "ENSG00000109911"="ELP4",
  "ENSG00000112964"="GHR",
  "ENSG00000113391"="ARH2A",
  "ENSG00000114315"="HES1",
  "ENSG00000120088"="CRHR1",
  "ENSG00000120519"="SLC10A7",
  "ENSG00000120738"="EGR1",
  "ENSG00000122877"="EGR2",
  "ENSG00000123106"="CCDC91",
  "ENSG00000125740"="FOSB",
  "ENSG00000132549"="VPS13B",
  "ENSG00000134532"="SOX5",
  "ENSG00000134775"="FHOD3",
  "ENSG00000138411"="HECW2",
  "ENSG00000138696"="BMPR1B",
  "ENSG00000140961"="OSGIN1",
  "ENSG00000144619"="CNTN4",
  "ENSG00000144642"="RBMS3",
  "ENSG00000145416"="MARCHF1",
  "ENSG00000145743"="FBXL17",
  "ENSG00000145996"="CDKAL1",
  "ENSG00000146350"="TBC1D32",
  "ENSG00000147202"="DIAPH2",
  "ENSG00000149489"="ROM1",
  "ENSG00000150471"="ADGRL3",
  "ENSG00000150672"="DLG2",
  "ENSG00000151276"="MAGI1",
  "ENSG00000151338"="MIPOL1",
  "ENSG00000151422"="FER",
  "ENSG00000152061"="RABGAP1L",
  "ENSG00000152208"="GRID2",
  "ENSG00000154655"="L3MBTL4",
  "ENSG00000156140"="ADAMTS3",
  "ENSG00000158528"="PPP1R9A",
  "ENSG00000160201"="U2AF1",
  "ENSG00000164008"="C1orf50",
  "ENSG00000165895"="ARHGAP42",
  "ENSG00000169855"="ROBO1",
  "ENSG00000169946"="ZFPM2",
  "ENSG00000170345"="FOS",
  "ENSG00000170579"="DLGAP1",
  "ENSG00000171723"="GPHN",
  "ENSG00000172296"="SPTLC3",
  "ENSG00000172780"="RAB43",
  "ENSG00000174891"="RSRC1",
  "ENSG00000178105"="DDX10",
  "ENSG00000182771"="GRID1",
  "ENSG00000184304"="PRKD1",
  "ENSG00000184903"="IMMP2L",
  "ENSG00000188352"="FOCAD",
  "ENSG00000196275"="GTF2IRD2",
  "ENSG00000204072"="ARMCX7P",
  "ENSG00000210151"="MT-TS1",
  "ENSG00000212643"="ZRSR2P1",
  "ENSG00000215421"="ZNF407",
  "ENSG00000219023"="RPS15AP19",
  "ENSG00000248643"="RBM14-RBM4",
  "ENSG00000259426"="KIF23-AS1",
  "ENSG00000277125"="PMS2P14",
  "ENSG00000277367"="NEK2P1",
  "ENSG00000284788"="ENSG00000284788"
)

# =============================================================
# 🔥 APPLY FIXED MAPPING (NO WARNINGS VERSION)
# =============================================================

rownames(expr_mat) <- gsub("\\..*", "", rownames(expr_mat))

mapped <- gene_map[rownames(expr_mat)]
mapped <- as.character(mapped)

missing_idx <- which(is.na(mapped))

if (length(missing_idx) > 0) {
  mapped[missing_idx] <- rownames(expr_mat)[missing_idx]
}

rownames(expr_mat) <- mapped

# =============================================================
# 6️⃣ Scale + cluster
# =============================================================

scaled_mat <- t(scale(t(expr_mat)))

hc <- hclust(dist(scaled_mat))
scaled_mat <- scaled_mat[hc$order, ]

# =============================================================
# 7️⃣ Force Control first
# =============================================================

col_order <- c(
  grep("^Control", colnames(scaled_mat), value = TRUE),
  grep("^Auxin", colnames(scaled_mat), value = TRUE)
)

scaled_mat <- scaled_mat[, col_order]

# =============================================================
# 8️⃣ Save heatmap
# =============================================================

output_file <- file.path(output_dir, "Control_vs_Auxin.png")

heat_colors <- colorRampPalette(c("blue","white","red"))(100)

pheatmap(
  scaled_mat,
  color = heat_colors,
  border_color = NA,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  labels_row = rownames(scaled_mat),
  fontsize_row = 8,
  fontsize_col = 12,
  main = "",
  filename = output_file,
  width = 10,
  height = 10
)
