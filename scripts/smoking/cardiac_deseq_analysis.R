library(Seurat)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(annotables)
library(EnhancedVolcano)
library(tibble)

# load counts and metadata
load('~/cardiac_seq_analysis/data/smoking/deseq_rawcounts.Robj')
load('~/cardiac_seq_analysis/data/smoking/deseq_metadata.Robj')

head(smoking_counts)
str(smoking_counts)

head(smoking_metadata)

# re-order columns of raw counts df
reorder_idx <- match(rownames(smoking_metadata), colnames(smoking_counts))
reorder_counts <- smoking_counts[, reorder_idx]

# create DESeq2 object
dds_smoking <- DESeqDataSetFromMatrix(countData = reorder_counts,
                                      colData = smoking_metadata,
                                      design = ~ Age + Source + Condition)

# log transform counts
vsd_all <- vst(dds_smoking, blind = T)

# plot pca
png("~/cardiac_seq_analysis/images/smoking/smoking_pca_age.png",
    width = 10*100, height = 3*100,
    res=100)
plotPCA(vsd_all, intgroup = "Age")
dev.off()

png("~/cardiac_seq_analysis/images/smoking/smoking_pca_source.png",
    width = 10*100, height = 3*100,
    res=100)
plotPCA(vsd_all, intgroup = "Source")
dev.off()

png("~/cardiac_seq_analysis/images/smoking/smoking_pca_condition.png",
    width = 10*100, height = 3*100,
    res=100)
plotPCA(vsd_all, intgroup = "Condition")
dev.off()

# run deseq analysis
dds_smoking <- DESeq(dds_smoking)

# plot dispersion estimates
png("~/cardiac_seq_analysis/images/smoking/smoking_disp_ests.png",
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
                         coef = 4,
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

png("~/cardiac_seq_analysis/images/smoking/smoking_volcano.png",
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

# extract the top 6 genes with padj values
top_de_genes <- smoking_res_sig %>%
  arrange(padj) %>%
  select(symbol, description, padj, log2FoldChange)
write.table(top_de_genes,
            "~/cardiac_seq_analysis/data/smoking/smoking_DESeq_markers.csv",
            sep=",", quote=F, row.names=F)
save(top_de_genes, file = '~/cardiac_seq_analysis/data/smoking/top_de_genes.Robj')

# volcano plot w/ overlapping genes
overlap_genes <- readRDS("~/cardiac_seq_analysis/data/smoking/overlap_top_de_genes.rds")

png("~/cardiac_seq_analysis/images/smoking/overlap_smoking_volcano.png",
    width = 12*100, height = 10*100,
    res=100)
EnhancedVolcano(smoking_res_all,
                lab = smoking_res_all$symbol,
                selectLab = overlap_genes$symbol,
                x = 'log2FoldChange',
                y = 'pvalue')
dev.off()
