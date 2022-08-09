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

load("~/cardiac_seq_analysis/data/smoking/choudhury_cells.Robj")

## export to DESeq2 format
cardiac_trimmed_counts <- choudhuryCells[["RNA"]]@counts

# convert sparse counts matrix to df
counts.df <- cardiac_trimmed_counts %>% as.matrix %>% t %>% as.data.frame
counts.df <- tibble::rownames_to_column(counts.df, "cellnames")
clusterassignments <- data.frame(choudhuryCells[["orig.ident"]])
clusterassignments <- tibble::rownames_to_column(clusterassignments, "cellnames")
counts.df <- merge(clusterassignments, counts.df, by = "cellnames")

rownames(counts.df) <- counts.df$cellnames
counts.df$cellnames <- NULL

counts.df <- counts.df %>% group_by(orig.ident) %>%
  summarise(across(everything(), ~ sum(., is.na(.), 0)))

counts.df <- as.data.frame(counts.df)
rownames(counts.df) <- counts.df$orig.ident
counts.df$orig.ident <- NULL

counts_t <- transpose(counts.df)
rownames(counts_t) <- colnames(counts.df)
colnames(counts_t) <- rownames(counts.df)

smoking_counts <- counts_t

save(smoking_counts, file = "~/cardiac_seq_analysis/data/smoking/choudhury_counts_by_ident.Robj")

## create metadata df

# Age column 
Age <- c("Middle", "Middle", "Old", "Old", "Middle")

# Condition column
Condition <- c("Smoking", "Smoking", "Control", "Smoking", "Control")

smoking_metadata <- data.frame(Condition, Age)

# add sample names
rownames(smoking_metadata) <- c("1156", "5111", "5828-2n", "5874", "936")

# save metadata
save(smoking_metadata, file = "~/cardiac_seq_analysis/data/smoking/choudhury_metadata.Robj")

# load counts and metadata
load('~/cardiac_seq_analysis/data/smoking/choudhury_counts_by_ident.Robj')
load('~/cardiac_seq_analysis/data/smoking/choudhury_metadata.Robj')

head(smoking_counts)
str(smoking_counts)

head(smoking_metadata)

# re-order columns of raw counts df
reorder_idx <- match(rownames(smoking_metadata), colnames(smoking_counts))
reorder_counts <- smoking_counts[, reorder_idx]

# create DESeq2 object
dds_smoking <- DESeqDataSetFromMatrix(countData = reorder_counts,
                                      colData = smoking_metadata,
                                      design = ~ Age + Condition)

# log transform counts
vsd_all <- vst(dds_smoking, blind = F)

# plot pca
png("~/cardiac_seq_analysis/images/smoking/ch_smoking_pca_age.png",
    width = 7*100, height = 5*100,
    res=100)
plotPCA(vsd_all, intgroup = "Age")
dev.off()

png("~/cardiac_seq_analysis/images/smoking/ch_smoking_pca_condition.png",
    width = 7*100, height = 5*100,
    res=100)
plotPCA(vsd_all, intgroup = "Condition")
dev.off()

# run deseq analysis
dds_smoking <- DESeq(dds_smoking)

# plot dispersion estimates
png("~/cardiac_seq_analysis/images/smoking/ch_smoking_disp_ests.png",
    width = 10*100, height = 10*100,
    res=100)
plotDispEsts(dds_smoking)
dev.off()

# calculate smoking results
smoking_res <- results(dds_smoking,
                       contrast = c("Condition", "Smoking", "Control"),
                       alpha = 0.05,
                       lfcThreshold = 0.1)

# lfc shrinkage
smoking_res <- lfcShrink(dds_smoking,
                         coef = 3,
                         res = smoking_res)

# MA plot
plotMA(smoking_res)

# summarize results
summary(smoking_res)

# add annotations
smoking_res_all <- data.frame(smoking_res) %>%
  rownames_to_column(var = "symbol")

smoking_res_all <- left_join(x = smoking_res_all,
                             y = grch38[, c("symbol", "description")],
                             by = "symbol")

# select significant genes with padj < 0.05
smoking_res_sig <- subset(smoking_res_all, padj < 0.05) %>%
  data.frame()

# volcano plot
smoking_res_all <- smoking_res_all %>%
  mutate(threshold = padj < 0.05)

png("~/cardiac_seq_analysis/images/smoking/ch_smoking_volcano.png",
    width = 12*100, height = 10*100,
    res=100)
EnhancedVolcano(smoking_res_all,
                lab = smoking_res_all$symbol,
                x = 'log2FoldChange',
                y = 'pvalue')
dev.off()

# ggplot(smoking_res_all) +
#   geom_point(aes(x = log2FoldChange, y = -log10(padj),
#                  color = threshold)) +
#   xlab("log2 fold change") +
#   ylab("-log10 adjusted p-value") +
#   ylim(0, 5) +
#   xlim(-5, 5) +
#   theme(legend.position = "none",
#         plot.title = element_text(size = rel(1.5), hjust = 0.5),
#         axis.title = element_text(size = rel(1.25)))

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
png("~/cardiac_seq_analysis/images/smoking/ch_smoking_heatmap.png",
    width = 12*100, height = 10*100,
    res=100)
pheatmap(sig_norm_counts,
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = select(smoking_metadata, Condition), 
         scale = "row")
dev.off()

# extract the top genes with padj values
ch_top_de_genes <- smoking_res_sig %>%
  arrange(padj) %>%
  select(symbol, padj, log2FoldChange)
write.table(ch_top_de_genes,
            "~/cardiac_seq_analysis/data/smoking/ch_smoking_DESeq_markers.csv",
            sep=",", quote=F, row.names=F)
save(ch_top_de_genes, file = '~/cardiac_seq_analysis/data/smoking/ch_top_de_genes.Robj')

# load de genes from integrated analysis
load("~/cardiac_seq_analysis/data/smoking/top_de_genes.Robj")

overlap_idx <- which(top_de_genes$symbol %in% ch_top_de_genes$symbol)
overlap_genes <- top_de_genes[overlap_idx,]

# save overlapping genes
saveRDS(overlap_genes, file = "~/cardiac_seq_analysis/data/smoking/overlap_top_de_genes.rds")

# save table w/o descriptions
write.table(overlap_genes[-c(2)], file = "~/cardiac_seq_analysis/data/smoking/overlap_top_de_genes.csv")
