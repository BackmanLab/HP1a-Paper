# =============================================

# COMBINED GO: HYPOXIA vs AUXIN+HYPOXIA

# =============================================

# ======================

# 0️⃣ Packages

# ======================

library(biomaRt)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)

# ======================

# 1️⃣ Paths

# ======================

base_dir <- "C:/Users/tkx1376/OneDrive - Northwestern University/RNA Sequencing/HP1a/rsem"
output_dir <- file.path(base_dir, "GO_BubblePlots")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

files <- list(
  "Hypoxia" = file.path(base_dir, "DESeq2_Hypoxia_vs_Control.csv"),
  "Auxin+Hypoxia" = file.path(base_dir, "DESeq2_AuxinHypoxia_vs_Control.csv")
)

# ======================

# 2️⃣ Mapping

# ======================

ensembl <- useEnsembl("genes", "hsapiens_gene_ensembl")

map_ids <- function(ids){
  getBM(
    attributes = c("ensembl_gene_id","entrezgene_id"),
    filters = "ensembl_gene_id",
    values = ids,
    mart = ensembl
  ) %>% filter(!is.na(entrezgene_id))
}

add_entrez <- function(df){
  m <- map_ids(rownames(df))
  df$EntrezID <- m$entrezgene_id[match(rownames(df), m$ensembl_gene_id)]
  df %>% filter(!is.na(EntrezID))
}

# ======================

# 3️⃣ GO enrichment

# ======================

run_go <- function(g){
  if(length(g)==0) return(NULL)
  enrichGO(g, OrgDb=org.Hs.eg.db, keyType="ENTREZID",
           ont="BP", pAdjustMethod="BH", qvalueCutoff=0.05, readable=TRUE)
}

prep <- function(obj, contrast, reg){
  if(is.null(obj)) return(NULL)
  as.data.frame(obj) %>%
    mutate(
      GeneRatio = sapply(strsplit(GeneRatio,"/"),
                         function(x) as.numeric(x[1])/as.numeric(x[2])),
      adj_pval = p.adjust,
      Contrast = contrast,
      Regulation = reg
    ) %>%
    select(Description, GeneRatio, adj_pval, Contrast, Regulation)
}

# ======================

# 4️⃣ Run BOTH conditions

# ======================

plot_df <- bind_rows(lapply(names(files), function(name){
  
  df <- read.csv(files[[name]], row.names = 1) %>% add_entrez()
  
  up <- df %>% filter(log2FoldChange > 0) %>% pull(EntrezID)
  down <- df %>% filter(log2FoldChange < 0) %>% pull(EntrezID)
  
  ego_up <- run_go(up)
  ego_down <- run_go(down)
  
  bind_rows(
    prep(ego_up, name, "Upregulated"),
    prep(ego_down, name, "Downregulated")
  )
}))

# ======================

# 5️⃣ Select shared top terms

# ======================

top_terms <- function(df, reg){
  df %>%
    filter(Regulation == reg) %>%
    arrange(adj_pval) %>%
    pull(Description) %>%
    unique() %>%
    head(15)
}

top_up <- top_terms(plot_df, "Upregulated")
top_down <- top_terms(plot_df, "Downregulated")

plot_up <- plot_df %>% filter(Regulation=="Upregulated", Description %in% top_up)
plot_down <- plot_df %>% filter(Regulation=="Downregulated", Description %in% top_down)

# factor order

plot_up$Contrast <- factor(plot_up$Contrast, levels = c("Hypoxia","Auxin+Hypoxia"))
plot_down$Contrast <- factor(plot_down$Contrast, levels = c("Hypoxia","Auxin+Hypoxia"))

plot_up$Description <- factor(plot_up$Description, levels = rev(top_up))
plot_down$Description <- factor(plot_down$Description, levels = rev(top_down))

# ======================

# 6️⃣ PLOTS

# ======================

p_up <- ggplot(plot_up, aes(x=Contrast, y=Description)) +
  geom_point(aes(size=GeneRatio, color=adj_pval), alpha=0.8) +
  scale_color_gradient(low="blue", high="red", name="Adjusted p-value") +
  scale_size(range=c(3,10), name="Gene Ratio") +
  guides(size=guide_legend(order=1), color=guide_colorbar(order=2)) +
  theme_minimal(base_size=18) +
  theme(plot.title=element_text(hjust=0.6)) +
  labs(title="GO Enrichment Analysis of Upregulated Genes in Hypoxia and Combo", x=NULL, y=NULL)

p_down <- ggplot(plot_down, aes(x=Contrast, y=Description)) +
  geom_point(aes(size=GeneRatio, color=adj_pval), alpha=0.8) +
  scale_color_gradient(low="blue", high="red", name="Adjusted p-value") +
  scale_size(range=c(3,10), name="Gene Ratio") +
  guides(size=guide_legend(order=1), color=guide_colorbar(order=2)) +
  theme_minimal(base_size=18) +
  theme(plot.title=element_text(hjust=0.6)) +
  labs(title="GO Enrichment Analysis of Downregulated Genes in Hypoxia and Combo", x=NULL, y=NULL)

# ======================

# 💾 Save

# ======================

ggsave(file.path(output_dir,"GO_Analysis_Hypoxia_Combo_Upregulated.png"), p_up, width=12, height=8)
ggsave(file.path(output_dir,"GO_Analysis_Hypoxia_Combo_Downregulated.png"), p_down, width=12, height=8)

print(p_up)
print(p_down)

library(dplyr)
library(pheatmap)
library(readr)
library(tibble)
library(grid)

# =============================================================
# 0️⃣ Load DESeq2 Results
# =============================================================

res_hypoxia <- read.csv("C:/Users/tkx1376/OneDrive - Northwestern University/RNA Sequencing/HP1a/rsem/DESeq2_Hypoxia_vs_Control.csv")
res_ah <- read.csv("C:/Users/tkx1376/OneDrive - Northwestern University/RNA Sequencing/HP1a/rsem/DESeq2_AuxinHypoxia_vs_Control.csv")

colnames(res_hypoxia)[1] <- "gene_id"
colnames(res_ah)[1] <- "gene_id"

sig_genes <- unique(c(
  res_hypoxia %>% filter(!is.na(padj), padj < 0.05) %>% pull(gene_id),
  res_ah %>% filter(!is.na(padj), padj < 0.05) %>% pull(gene_id)
))

# =============================================================
# 1️⃣ TPM FILES
# =============================================================

files_ctrl <- c(
  "C:/Users/tkx1376/OneDrive - Northwestern University/RNA Sequencing/HP1a/seq_raw_output/Control_1_S25_L005_R1_001Aligned.genes.results",
  "C:/Users/tkx1376/OneDrive - Northwestern University/RNA Sequencing/HP1a/seq_raw_output/Control_2_S26_L005_R1_001Aligned.genes.results",
  "C:/Users/tkx1376/OneDrive - Northwestern University/RNA Sequencing/HP1a/seq_raw_output/Control_3_S27_L005_R1_001Aligned.genes.results"
)


files_hypoxia <- c(
  "C:/Users/tkx1376/OneDrive - Northwestern University/RNA Sequencing/HP1a/seq_raw_output/Hypoxia_1_S31_L005_R1_001Aligned.genes.results",
  "C:/Users/tkx1376/OneDrive - Northwestern University/RNA Sequencing/HP1a/seq_raw_output/Hypoxia_2_S32_L005_R1_001Aligned.genes.results",
  "C:/Users/tkx1376/OneDrive - Northwestern University/RNA Sequencing/HP1a/seq_raw_output/Hypoxia_3_S33_L005_R1_001Aligned.genes.results"
)

files_ah <- c(
  "C:/Users/tkx1376/OneDrive - Northwestern University/RNA Sequencing/HP1a/seq_raw_output/Auxin_Hypoxia_1_S34_L005_R1_001Aligned.genes.results",
  "C:/Users/tkx1376/OneDrive - Northwestern University/RNA Sequencing/HP1a/seq_raw_output/Auxin_Hypoxia_2_S35_L005_R1_001Aligned.genes.results",
  "C:/Users/tkx1376/OneDrive - Northwestern University/RNA Sequencing/HP1a/seq_raw_output/Auxin_Hypoxia_3_S36_L005_R1_001Aligned.genes.results"
)
# =============================================================
# 2️⃣ SAFE READER (NO MERGE, NO DUPLICATION POSSIBLE)
# =============================================================

read_group <- function(files, label){
  
  mats <- lapply(seq_along(files), function(i){
    
    df <- read_tsv(files[i], show_col_types = FALSE)
    df$gene_id <- gsub("\\..*", "", df$gene_id)
    
    df <- df[, c("gene_id", "TPM")]
    
    # CLEAN LABELS (NO DASHES, HUMAN READABLE)
    pretty_label <- if(label == "AuxinHypoxia") {
      "Auxin + Hypoxia"
    } else {
      label
    }
    
    colnames(df)[2] <- paste0(pretty_label, " ", i)
    
    df
  })
  
  mat_list <- lapply(mats, function(x){
    m <- as.data.frame(x)
    rownames(m) <- m$gene_id
    m <- m[, -1, drop = FALSE]
    m
  })
  
  do.call(cbind, mat_list)
}
# =============================================================
# 3️⃣ BUILD MATRICES
# =============================================================

ctrl_mat <- read_group(files_ctrl, "Control")
hyp_mat  <- read_group(files_hyp, "Hypoxia")
ah_mat   <- read_group(files_ah, "AuxinHypoxia")

# =============================================================
# 🚨 HARD SAFETY CHECK (catches duplication instantly)
# =============================================================

ctrl_mat <- ctrl_mat[, !duplicated(colnames(ctrl_mat))]
hyp_mat  <- hyp_mat[, !duplicated(colnames(hyp_mat))]
ah_mat   <- ah_mat[, !duplicated(colnames(ah_mat))]

# =============================================================
# 4️⃣ ALIGN GENES (CRITICAL BIOLOGICAL STEP)
# =============================================================

common_genes <- Reduce(intersect, list(
  rownames(ctrl_mat),
  rownames(hyp_mat),
  rownames(ah_mat)
))

ctrl_mat <- ctrl_mat[common_genes, , drop = FALSE]
hyp_mat  <- hyp_mat[common_genes, , drop = FALSE]
ah_mat   <- ah_mat[common_genes, , drop = FALSE]

# FORCE UNIQUE COLUMN NAMES BEFORE BINDING
colnames(ctrl_mat) <- make.unique(colnames(ctrl_mat))
colnames(hyp_mat)  <- make.unique(colnames(hyp_mat))
colnames(ah_mat)   <- make.unique(colnames(ah_mat))

# SAFE BIND (NO AUTO .1 GENERATION)
expr_mat <- cbind(
  ctrl_mat,
  hyp_mat,
  ah_mat
)

# FINAL SAFETY CLEAN (removes any leftover accidental duplicates)
expr_mat <- expr_mat[, !duplicated(colnames(expr_mat))]

# =============================================================
# 5️⃣ FILTER SIGNIFICANT GENES
# =============================================================

expr_mat <- expr_mat[rownames(expr_mat) %in% sig_genes, ]

# =============================================================
# 6️⃣ SELECT TOP VARIABLE GENES
# =============================================================

gene_variance <- apply(expr_mat, 1, var)

top_n <- min(65, nrow(expr_mat))
expr_mat <- expr_mat[order(gene_variance, decreasing = TRUE)[1:top_n], ]

# =============================================================
# 7️⃣ EXPORT ENSEMBL IDS
# =============================================================

output_dir <- "C:/Users/tkx1376/OneDrive - Northwestern University/RNA Sequencing/HP1a/Final Heatmaps"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

write.csv(
  data.frame(Ensembl_ID = rownames(expr_mat)),
  file.path(output_dir, "Combo_EnsemblIDs.csv"),
  row.names = FALSE
)

gene_map <- c(
  "ENSG00000026025"="VIM",
  "ENSG00000065978"="YBX1",
  "ENSG00000070756"="PABPC1",
  "ENSG00000080824"="HSP90AA1",
  "ENSG00000087086"="FTL",
  "ENSG00000096384"="HSP90AB1",
  "ENSG00000100097"="LGALS1",
  "ENSG00000100867"="DHRS2",
  "ENSG00000102317"="RBM3",
  "ENSG00000104529"="EEF1D",
  "ENSG00000105193"="RPS16",
  "ENSG00000106211"="HSPB1",
  "ENSG00000108821"="COL1A1",
  "ENSG00000111640"="GAPDH",
  "ENSG00000111669"="TPI1",
  "ENSG00000113140"="SPARC",
  "ENSG00000118523"="CCN2",
  "ENSG00000120885"="CLU",
  "ENSG00000123416"="TUBA1B",
  "ENSG00000133112"="TPT1",
  "ENSG00000134333"="LDHA",
  "ENSG00000137331"="IER3",
  "ENSG00000137818"="RPLP1",
  "ENSG00000140988"="RPS2",
  "ENSG00000142534"="RPS11",
  "ENSG00000142541"="RPL13A",
  "ENSG00000142871"="CCN1",
  "ENSG00000145592"="RPL37",
  "ENSG00000147403"="RPL10",
  "ENSG00000149257"="SERPINH1",
  "ENSG00000149925"="ALDOA",
  "ENSG00000156508"="EEF1A1",
  "ENSG00000163191"="S100A11",
  "ENSG00000163430"="FSTL1",
  "ENSG00000167996"="FTH1",
  "ENSG00000174444"="RPL4",
  "ENSG00000182718"="ANXA2",
  "ENSG00000182774"="RPS17",
  "ENSG00000184009"="ACTG1",
  "ENSG00000189060"="H1-0",
  "ENSG00000196205"="EEF1A1P5",
  "ENSG00000197903"="H2BC12",
  "ENSG00000198695"="MT-ND6",
  "ENSG00000198727"="MT-CYB",
  "ENSG00000198763"="MT-ND2",
  "ENSG00000198786"="MT-ND5",
  "ENSG00000198804"="MT-CO1",
  "ENSG00000198840"="MT-ND3",
  "ENSG00000198886"="MT-ND4",
  "ENSG00000205542"="TMSB4X",
  "ENSG00000210082"="MT-RNR2",
  "ENSG00000210135"="MT-TN",
  "ENSG00000210140"="MT-TC",
  "ENSG00000210144"="MT-TY",
  "ENSG00000210151"="MT-TS1",
  "ENSG00000212907"="MT-ND4L",
  "ENSG00000213626"="LBH",
  "ENSG00000229117"="RPL41",
  "ENSG00000231500"="RPS18",
  "ENSG00000237550"="ENSG00000237550",
  "ENSG00000240342"="RPS2P5",
  "ENSG00000279864"="NCOR1P4",
  "ENSG00000280614"="ENSG00000280614",
  "ENSG00000280800"="ENSG00000280800",
  "ENSG00000281181"="ENSG00000281181"
)

rownames(expr_mat) <- gsub("\\..*", "", rownames(expr_mat))

mapped <- gene_map[rownames(expr_mat)]
mapped <- as.character(mapped)

mapped[is.na(mapped)] <- rownames(expr_mat)[is.na(mapped)]

# 🔥 CRITICAL FIX
rownames(expr_mat) <- make.unique(mapped)

# =============================================================
# 8️⃣ SCALE (ROW Z-SCORE)
# =============================================================

scaled_mat <- t(scale(t(expr_mat)))

scaled_mat[scaled_mat > 2] <- 2
scaled_mat[scaled_mat < -2] <- -2

# =============================================================
# 9️⃣ BIOLOGICAL COLUMN ORDER
# =============================================================

scaled_mat <- scaled_mat[, c(
  grep("Control", colnames(scaled_mat), value = TRUE),
  grep("Hypoxia", colnames(scaled_mat), value = TRUE),
  grep("Auxin \\+ Hypoxia", colnames(scaled_mat), value = TRUE)
)]

# =============================================================
# 🔟 FINAL HEATMAP (BIOLOGICAL STRUCTURE)
# =============================================================
# clean whitespace
colnames(scaled_mat) <- trimws(colnames(scaled_mat))

# 🔥 REMOVE the unwanted ".1" duplicates directly
scaled_mat <- scaled_mat[, !grepl("\\.1$", colnames(scaled_mat))]

col_order <- c(
  grep("^Control", colnames(scaled_mat), value = TRUE),
  grep("^Hypoxia", colnames(scaled_mat), value = TRUE),
  grep("^Auxin \\+ Hypoxia", colnames(scaled_mat), value = TRUE)
)

scaled_mat <- scaled_mat[, col_order]
heat_colors <- colorRampPalette(c("blue","white","red"))(100)

pheatmap(
  scaled_mat,
  
  cluster_rows = TRUE,        # biological grouping
  cluster_cols = FALSE,
  
  treeheight_row = 0,
  treeheight_col = 0,
  border_color = NA,
  
  gaps_row = NULL,
  gaps_col = NULL,
  
  show_colnames = TRUE,
  show_rownames = TRUE,
  
  color = heat_colors,
  fontsize_row = 7,
  fontsize_col = 10,
  
  legend_breaks = seq(-2, 2, by = 0.5),
  legend_labels = seq(-2, 2, by = 0.5),
  
  filename = file.path(output_dir, "FINAL_3Condition_Heatmap.png"),
  width = 11,
  height = 10
)

