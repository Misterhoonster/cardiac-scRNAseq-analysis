library(Seurat)
library(DESeq2)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(data.table)
library(annotables)
library(EnhancedVolcano)
library(tibble)
library(scuttle)

## load data

# export to DESeq2 format
TE_list <- readRDS("~/cardiac_seq_analysis/data/TE/TE_list.rds")

te_smoking <- readRDS("~/cardiac_seq_analysis/data/TE/smoking_tes_clustered.rds")

# convert to SCE
te_smoking <- as.SingleCellExperiment(te_smoking)
colData(te_smoking)

# aggregate by cluster and treatment
sum_by <- 
  c("orig.ident")

summed <- 
  aggregateAcrossCells(te_smoking, id=colData(te_smoking)[,sum_by],
                       use.altexps = NULL)

# extract counts matrix
smoking_counts <- assay(summed, "counts")

# te_smoking <- subset(te_smoking, features = rownames(te_smoking)[rownames(te_smoking) %in% TE_list$final_TEs])
# te_smoking_counts <- te_smoking[["SCT"]]@counts
# 
# # reduce size of counts so it can be converted to matrix
# te_smoking <- CreateSeuratObject(te_smoking_counts, min.cells = 5)
# te_smoking_counts <- te_smoking[["RNA"]]@counts
# 
# # convert sparse counts matrix to df
# counts.df <- te_smoking_counts %>% as.matrix %>% t %>% as.data.frame
# counts.df <- tibble::rownames_to_column(counts.df, "cellnames")
# clusterassignments <- data.frame(te_smoking[["orig.ident"]])
# clusterassignments <- tibble::rownames_to_column(clusterassignments, "cellnames")
# counts.df <- merge(clusterassignments, counts.df, by = "cellnames")
# 
# rownames(counts.df) <- counts.df$cellnames
# counts.df$cellnames <- NULL
# 
# counts.df <- counts.df %>% group_by(orig.ident) %>%
#   summarise(across(everything(), ~ sum(., is.na(.), 0)))
# 
# counts.df <- as.data.frame(counts.df)
# rownames(counts.df) <- counts.df$orig.ident
# counts.df$orig.ident <- NULL
# 
# counts_t <- as.data.frame(t(counts.df))
# 
# smoking_counts <- counts_t

saveRDS(smoking_counts, file = "~/cardiac_seq_analysis/data/TE/smoking_deseq_counts.rds")

## create metadata df

# Age column
Age <- c("Young", "Young", "Old", "Middle", "Middle", "Old", "Middle", "Old",
         "Middle", "Middle")

# Sex column
Sex <- c("F", "F", "F", "M", "F", "M", "F", "F", "F", "M")

# Condition column
Condition <- c("Control", "Control", "Control", "Control", "Control", "Control",
               "Smoking", "Smoking", "Smoking", "Smoking")

smoking_metadata <- data.frame(Condition, Sex, Age)

# add sample names
rownames(smoking_metadata) <- c("4638-2n-1", "4638-heart", "5828-heart", "5087-heart-5prime",
                                "936-heart-5prime", "5919-2n-all", "1156-heart-5prime", 
                                "5874-heart-5prime", "5111-heart-5prime", "604-heart-5prime")

# save metadata
saveRDS(smoking_metadata, file = "~/cardiac_seq_analysis/data/TE/smoking_deseq_metadata.rds")

# load counts and metadata
smoking_counts <- readRDS('~/cardiac_seq_analysis/data/TE/smoking_deseq_counts_sct.rds')
smoking_metadata <- readRDS('~/cardiac_seq_analysis/data/TE/smoking_deseq_metadata.rds')

head(smoking_counts)
str(smoking_counts)

head(smoking_metadata)

# re-order columns of raw counts df
reorder_idx <- match(rownames(smoking_metadata), colnames(smoking_counts))
reorder_counts <- smoking_counts[, reorder_idx]

# create DESeq2 object
dds_smoking <- DESeqDataSetFromMatrix(countData = reorder_counts,
                                      colData = smoking_metadata,
                                      design = ~ Age + Sex + Condition)

# setwd to images folder
setwd("~/cardiac_seq_analysis/images/TE/deseq2/")

# log transform counts
vsd_all <- vst(dds_smoking, blind = F)

# plot pca
pc1 <- plotPCA(vsd_all, intgroup = "Age") + geom_label_repel(aes(label = name))
ggsave("smoking_pca_age.png",
       plot = pc1,
       width = 5, height = 5)

pc2 <- plotPCA(vsd_all, intgroup = "Sex") + geom_label_repel(aes(label = name))
ggsave("smoking_pca_sex.png",
       plot = pc2,
       width = 5, height = 5)

pc3 <- plotPCA(vsd_all, intgroup = "Condition") + geom_label_repel(aes(label = name))
ggsave("smoking_pca_condition.png",
       plot = pc3,
       width = 5, height = 5)

# run deseq analysis
dds_smoking <- DESeq(dds_smoking)

# plot dispersion estimates
plotDispEsts(dds_smoking)
ggsave("smoking_disp_ests.png",
       width = 10, height = 10)

# calculate smoking results
smoking_res <- results(dds_smoking,
                       contrast = c("Condition", "Smoking", "Control"),
                       alpha = 0.05)

# lfc shrinkage
smoking_res <- lfcShrink(dds_smoking,
                         coef = 5,
                         res = smoking_res)

# MA plot
plotMA(smoking_res)

# summarize results
summary(smoking_res)

# add annotations
smoking_res_all <- data.frame(smoking_res) %>%
  rownames_to_column(var = "symbol")

# smoking_res_all <- left_join(x = smoking_res_all,
#                              y = grch38[, c("symbol", "description")],
#                              by = "symbol")

# select significant genes with padj < 0.05
smoking_res_sig <- subset(smoking_res_all, padj < 0.05) %>%
  data.frame() %>% arrange(padj, log2FoldChange)
write.csv(smoking_res_sig, file = "~/cardiac_seq_analysis/data/TE/smoking_te_markers_deseq2.csv")


smoking_sig_te <- smoking_res_sig %>% filter(symbol %in% TE_list$final_TEs) %>% arrange(padj, log2FoldChange)
write.csv(smoking_sig_te, file = "~/cardiac_seq_analysis/data/TE/smoking_te_markers_deseq2_with_genes.csv")

# volcano plot
smoking_res_all <- smoking_res_all %>%
  mutate(threshold = padj < 0.05)

volcano <- EnhancedVolcano(smoking_res_all,
                lab = smoking_res_all$symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05)
ggsave("smoking_volcano.png",
       plot = volcano,
       width = 12, height = 10)

## heatmap
# Determine the size factors to use for normalization
dds_smoking <- estimateSizeFactors(dds_smoking)

# Extract the normalized counts
smoking_normalized_counts <- counts(dds_smoking, normalized=TRUE)

# Subset normalized counts to significant genes
sig_norm_counts <- smoking_normalized_counts[smoking_res_sig$symbol, ]

# Choose heatmap color palette
heat_colors <- brewer.pal(n = 6, name = "YlOrRd")

# Plot heatmap
heatmap <- pheatmap(sig_norm_counts,
         color = heat_colors, 
         cluster_rows = T,
         cluster_cols = F,
         show_rownames = F,
         annotation = select(smoking_metadata, Condition), 
         scale = "row")
ggsave("smoking_heatmap.png",
       plot = heatmap,
       width = 12, height = 10)

