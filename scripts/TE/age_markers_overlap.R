# set wd to data folder
setwd("~/cardiac_seq_analysis/data/TE/")

# load in results from DESeq + Seurat analyses
y <- read.csv("young_age_te_markers.csv", row.names = "X") %>% rownames_to_column("symbol")
m <- read.csv("middle_age_te_markers.csv", row.names = "X") %>% rownames_to_column("symbol")
o <- read.csv("old_age_te_markers.csv", row.names = "X") %>% rownames_to_column("symbol")

ym <- read.csv("age_te_markers_deseq2_ym.csv", row.names = "X")
yo <- read.csv("age_te_markers_deseq2_yo.csv", row.names = "X")
mo <- read.csv("age_te_markers_deseq2_mo.csv", row.names = "X")

# preserve markers found in both analyses
y_overlap_names <- intersect(union(ym$symbol, yo$symbol), y$symbol)
m_overlap_names <- intersect(union(ym$symbol, mo$symbol), m$symbol)
o_overlap_names <- intersect(union(yo$symbol, mo$symbol), o$symbol)

y_overlap <- y %>% filter(symbol %in% y_overlap_names) %>% arrange(p_val_adj, desc(avg_log2FC))
m_overlap <- m %>% filter(symbol %in% m_overlap_names) %>% arrange(p_val_adj, desc(avg_log2FC))
o_overlap <- o %>% filter(symbol %in% o_overlap_names) %>% arrange(p_val_adj, desc(avg_log2FC))

# save overlapping markers
write.csv(y_overlap, "results/young_overlap_markers.csv")
write.csv(m_overlap, "results/middle_overlap_markers.csv")
write.csv(o_overlap, "results/old_overlap_markers.csv")
