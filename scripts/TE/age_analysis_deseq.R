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

te_age <- readRDS("~/cardiac_seq_analysis/data/TE/age_te.rds")

# convert to SCE
te_age <- as.SingleCellExperiment(te_age)
colData(te_age)

# aggregate by cluster and treatment
sum_by <- 
  c("orig.ident")

summed <- 
  aggregateAcrossCells(te_age, id=colData(te_age)[,sum_by],
                       use.altexps = NULL)

# extract counts matrix
age_counts <- assay(summed, "counts")

# # te_age <- subset(te_age, features = rownames(te_age)[rownames(te_age) %in% TE_list$final_TEs])
# 
# # reduce size of counts so it can be converted to matrix
# # te_age <- CreateSeuratObject(te_age_counts, min.cells = 6)
# te_age_counts <- te_age[["RNA"]]@counts
# 
# # convert sparse counts matrix to df
# counts.df <- te_age_counts %>% as.matrix %>% t %>% as.data.frame
# counts.df <- tibble::rownames_to_column(counts.df, "cellnames")
# clusterassignments <- data.frame(te_age[["orig.ident"]])
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
# age_counts <- counts_t

saveRDS(age_counts, file = "~/cardiac_seq_analysis/data/TE/age_deseq_counts.rds")

## create metadata df

# Age column
Age <- c("Young", "Young", "Old", "Middle", "Middle", "Old")

# Sex column
Sex <- c("F", "F", "F", "M", "F", "M")

age_metadata <- data.frame(Age, Sex)

# add sample names
rownames(age_metadata) <- c("4638-2n-1", "4638-heart", "5828-heart", "5087-heart-5prime",
                                "936-heart-5prime", "5919-2n-all")

# save metadata
saveRDS(age_metadata, file = "~/cardiac_seq_analysis/data/TE/age_deseq_metadata.rds")

# load counts and metadata
age_counts <- readRDS('~/cardiac_seq_analysis/data/TE/age_deseq_counts.rds')
age_metadata <- readRDS('~/cardiac_seq_analysis/data/TE/age_deseq_metadata.rds')

head(age_counts)
str(age_counts)

head(age_metadata)

# set wd to images folder
setwd("~/cardiac_seq_analysis/images/TE/deseq2/")

# re-order columns of raw counts df
reorder_idx <- match(rownames(age_metadata), colnames(age_counts))
reorder_counts <- age_counts[, reorder_idx]

# create DESeq2 object
dds_age <- DESeqDataSetFromMatrix(countData = reorder_counts,
                                      colData = age_metadata,
                                      design = ~ Sex + Age)

# log transform counts
vsd_all <- vst(dds_age, blind = F)

# plot pca
pc1 <- plotPCA(vsd_all, intgroup = "Sex") + geom_label_repel(aes(label = name))
ggsave("age_sex_pca.png",
       plot = pc1,
       width = 5, height = 5)

pc2 <- plotPCA(vsd_all, intgroup = "Age") + geom_label_repel(aes(label = name))
ggsave("age_pca.png",
       plot = pc2,
       width = 5, height = 5)

# run deseq analysis
dds_age <- DESeq(dds_age)

# plot dispersion estimates
disp <- plotDispEsts(dds_age)
ggsave("age_disp_ests.png",
       width = 10, height = 10)

# calculate smoking results
age_ym <- results(dds_age,
                   contrast = c("Age", "Young", "Middle"),
                       alpha = 0.05,
                       lfcThreshold = 0.1)

age_yo <- results(dds_age,
                   contrast = c("Age", "Young", "Old"),
                   alpha = 0.05,
                   lfcThreshold = 0.1)

age_mo <- results(dds_age,
                   contrast = c("Age", "Middle", "Old"),
                   alpha = 0.05,
                   lfcThreshold = 0.1)

# lfc shrinkage
age_res_ym <- lfcShrink(dds_age,
                         contrast = c("Age", "Young", "Middle"),
                         res = age_ym,
                         type = "normal")

age_res_yo <- lfcShrink(dds_age,
                         contrast = c("Age", "Young", "Old"),
                         res = age_yo,
                         type = "normal")

age_res_mo <- lfcShrink(dds_age,
                         contrast = c("Age", "Middle", "Old"),
                         res = age_mo,
                         type = "normal")

# MA plot
plotMA(age_res_ym)
plotMA(age_res_yo)
plotMA(age_res_mo)

# summarize results
summary(age_res_ym)
summary(age_res_yo)
summary(age_res_mo)

# add annotations
age_res_ym <- data.frame(age_res_ym) %>%
  rownames_to_column(var = "symbol")
age_res_yo <- data.frame(age_res_yo) %>%
  rownames_to_column(var = "symbol")
age_res_mo <- data.frame(age_res_mo) %>%
  rownames_to_column(var = "symbol")

age_sig_ym <- age_res_ym %>% filter(symbol %in% TE_list$final_TEs, padj < 0.05) %>% arrange(padj)
write.csv(age_sig_ym, file = "~/cardiac_seq_analysis/data/TE/age_te_markers_deseq2_ym.csv")

age_sig_yo <- age_res_yo %>% filter(symbol %in% TE_list$final_TEs, padj < 0.05) %>% arrange(padj)
write.csv(age_sig_yo, file = "~/cardiac_seq_analysis/data/TE/age_te_markers_deseq2_yo.csv")

age_sig_mo <- age_res_mo %>% filter(symbol %in% TE_list$final_TEs, padj < 0.05) %>% arrange(padj)
write.csv(age_sig_mo, file = "~/cardiac_seq_analysis/data/TE/age_te_markers_deseq2_mo.csv")

# volcano plot
age_res_all <- age_res_all %>%
  mutate(threshold = padj < 0.05)

EnhancedVolcano(age_res_all,
                lab = age_res_all$symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05)
ggsave("age_volcano.png",
       width = 12, height = 10)

## heatmap
# Determine the size factors to use for normalization
dds_age <- estimateSizeFactors(dds_age)

# Extract the normalized counts
age_normalized_counts <- counts(dds_age, normalized=TRUE)

# Choose heatmap color palette
heat_colors <- brewer.pal(n = 6, name = "YlOrRd")

# Plot heatmap

heat_ym <- pheatmap(age_normalized_counts[age_sig_ym$symbol, ],
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = select(age_metadata, Age), 
         scale = "row")
ggsave("age_heatmap_ym.png",
       plot = heat_ym,
    width = 12, height = 10)

# too few markers to heatmap
pheatmap(age_normalized_counts[age_sig_yo$symbol,,drop=F],
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = select(age_metadata, Age), 
         scale = "row")
ggsave("age_heatmap_yo.png",
       width = 12, height = 10)

heat_mo <- pheatmap(age_normalized_counts[age_sig_mo$symbol, ],
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = select(age_metadata, Age), 
         scale = "row")
ggsave("age_heatmap_mo.png",
       plot = heat_mo,
       width = 12, height = 10)
