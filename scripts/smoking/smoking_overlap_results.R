library(Seurat)
library(ggplot2)
library(patchwork)

# check for overlaps between all smoking DEG lists

#### let's just do GO enrichment don't have to deal with labeling every single gene

# load all DEG data
ch_seurat <- read.table("~/cardiac_seq_analysis/data/smoking/results/ch_seurat.csv", sep=",", header = T)
ch_deseq <- read.table("~/cardiac_seq_analysis/data/smoking/results/ch_deseq.csv", sep=",", header = T)
integrated_seurat <- read.table("~/cardiac_seq_analysis/data/smoking/results/integrated_seurat.csv", sep=",", header = T)
integrated_deseq <- read.table("~/cardiac_seq_analysis/data/smoking/results/integrated_deseq.csv", sep=",", header = T)

# display data
head(ch_seurat)
head(ch_deseq)
head(integrated_seurat)
head(integrated_deseq)

# find overlapping genes
integrated_overlap <- integrated_seurat[which(integrated_seurat$symbol %in% integrated_deseq$symbol),]

ch_overlap <- ch_seurat[which(ch_seurat$symbol %in% ch_deseq$symbol),]

all_overlap <- integrated_overlap[which(integrated_overlap$symbol %in% ch_overlap$symbol),]

# save all results
write.table(ch_seurat, '~/cardiac_seq_analysis/data/smoking/results/ch_seurat.csv',
                         sep=",", quote = F, row.names = F)
write.table(ch_deseq, '~/cardiac_seq_analysis/data/smoking/results/ch_deseq.csv',
                         sep=",", quote = F, row.names = F)
write.table(integrated_seurat, '~/cardiac_seq_analysis/data/smoking/results/integrated_seurat.csv',
                         sep=",", quote = F, row.names = F)
write.table(integrated_deseq, '~/cardiac_seq_analysis/data/smoking/results/integrated_deseq.csv',
                         sep=",", quote = F, row.names = F)

write.table(integrated_overlap, '~/cardiac_seq_analysis/data/smoking/results/integrated_overlap.csv',
            sep=",", quote = F, row.names = F)
write.table(ch_overlap, '~/cardiac_seq_analysis/data/smoking/results/ch_overlap.csv',
            sep=",", quote = F, row.names = F)
write.table(all_overlap, '~/cardiac_seq_analysis/data/smoking/results/all_overlap.csv',
            sep=",", quote = F, row.names = F)


### wikipathways analysis w/ clusterprofiler
## download human genome database
library(org.Hs.eg.db)
library(clusterProfiler)
library(cowplot)
library(enrichplot)
hs <- org.Hs.eg.db

setwd("~/cardiac_seq_analysis/images/smoking/results/")

add_pathways <- function(df, name) {
  print(paste("Adding pathways for", name))
  
  # split list into up and down regulated
  marker_up <- df[df$avg_log2FC > 0,]$symbol
  marker_down <- df[df$avg_log2FC < 0,]$symbol
  
  # change symbol names to official values
  og <- c("MT-", "^CYB$", "CO1", "CO2", "CO3", "H1FX", "LINC01268", "DUSP27",
          "C5orf17", "LINC01481", "LINC00599", "C15orf41", "GRASP")
  new <- c("", "CYTB", "COX1", "COX2", "COX3", "H1-10", "MROCKI", "DUSP29",
           "LINC02899", "PRANCR", "MIR124-1HG", "CDIN1", "TAMALIN")
  
  for (i in 1:length(og)) {
    marker_up <- gsub(og[i], new[i], marker_up)
    marker_down <- gsub(og[i], new[i], marker_down)
  }
  
  # GO enrichment analysis
  ego_up <- enrichGO(gene          = marker_up,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "CC",
                  keyType       = "SYMBOL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
  ego_down <- enrichGO(gene          = marker_down,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "CC",
                     keyType       = "SYMBOL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)
  
  # save up and down GOs
  write.table(head(ego_up, 50), quote = F, sep = ",", row.names = F,
          file = paste0("~/cardiac_seq_analysis/data/smoking/results/",
                        name, "_GO_up.csv"))
  write.table(head(ego_down, 50), quote = F, sep = ",", row.names = F,
          file = paste0("~/cardiac_seq_analysis/data/smoking/results/",
                        name, "_GO_down.csv"))
  
  # create graph of sig pathways
  ego_up <- pairwise_termsim(ego_up)
  graph_up <- emapplot(ego_up, cex_label_category = 0.8)
  ggsave(paste0(name, "_graph_up.png"), plot = graph_up, width = 10, height = 10)
  
  ego_down <- pairwise_termsim(ego_down)
  graph_down <- emapplot(ego_down, cex_label_category = 0.8)
  ggsave(paste0(name, "_graph_down.png"), plot = graph_down, width = 10, height = 10)
  
  # dot plots
  dotplot(ego_up)
  ggsave(paste0(name, "_dot_down.png"), width = 10, height = 10)
  
  dotplot(ego_down)
  ggsave(paste0(name, "_dot_up.png"), width = 10, height = 10)
}

combined_df <- list(list(ch_seurat, "ch_seurat"),
                    list(integrated_seurat, "integrated_seurat"))

for (df in combined_df) {
  add_pathways(df[[1]], df[[2]])
}


## look for distribution of RGCC and FABP4

# load choudhury cells Seurat object
load("~/cardiac_seq_analysis/data/smoking/choudhury_cells.Robj")
setwd("~/cardiac_seq_analysis/images/smoking/")

# Plot cell type clusters
DimPlot(choudhuryCells, reduction = "umap", label = T)

# Violin plot of sig genes

plots <- VlnPlot(choudhuryCells, features = c("RGCC", "FABP4"), group.by = "Condition")

for(i in 1:length(plots)) {
  plots[[i]] <- plots[[i]] + geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_minimal()
}
ggsave("rgcc_fabp4_vln_smoking_cum.png", width = 10, height = 5)

p1 <- VlnPlot(subset(choudhuryCells, Condition == "Smoking"), features = c("RGCC", "FABP4")) + 
ggsave("rgcc_fabp4_vln_smoking.png", width = 10, height = 5)

p2 <- VlnPlot(subset(choudhuryCells, Condition == "Control"), features = c("RGCC", "FABP4"))
ggsave("rgcc_fabp4_vln_control.png", width = 10, height = 5)

# Dot plot of sig genes
DotPlot(choudhuryCells, features = c("RGCC", "FABP4"))
ggsave("rgcc_fabp4_dot.png", width = 7, height = 5)

# Feature plot of sig genes
features <- FeaturePlot(choudhuryCells, features = c("RGCC", "FABP4"), label = T,
            split.by = "Condition")
ggsave("rgcc_fabp4_feature.png", plot = features, width = 10, height = 10)

## check sig pathways
# load pathways
ch_pathways_down <- read.table("~/cardiac_seq_analysis/data/smoking/results/ch_seurat_GO_down.csv",
                               sep = ",", header = T)
ch_pathways_up <- read.table("~/cardiac_seq_analysis/data/smoking/results/ch_seurat_GO_up.csv",
                             sep = ",", header = T)

integrated_pathways_down <- read.table("~/cardiac_seq_analysis/data/smoking/results/integrated_seurat_GO_down.csv",
                               sep = ",", header = T)
integrated_pathways_up <- read.table("~/cardiac_seq_analysis/data/smoking/results/integrated_seurat_GO_up.csv",
                             sep = ",", header = T)

# find overlapping pathways
overlap_pathways_down <- ch_pathways_down[ch_pathways_down$ID %in% integrated_pathways_down$ID,]
overlap_pathways_up <- ch_pathways_up[ch_pathways_up$ID %in% integrated_pathways_up$ID,]

overlap_pathways_down$geneID <- gsub("ND2", "MT-ND2", overlap_pathways_down$geneID)

# reformat gene lists for each pathway
overlap_pathways_up <- overlap_pathways_up %>% transmute(ID, Description, gene_list = strsplit(geneID, "/"))
overlap_pathways_down <- overlap_pathways_down %>% transmute(ID, Description, gene_list = strsplit(geneID, "/"))

overlap_pathways_up$Description <- gsub(" ", ".", overlap_pathways_up$Description)
overlap_pathways_down$Description <- gsub(" ", ".", overlap_pathways_down$Description)
overlap_pathways_down$Description <- gsub("-", ".", overlap_pathways_down$Description)

# add composite scores for each gene set
for (i in 1:nrow(overlap_pathways_up)) {
  choudhuryCells <- AddModuleScore(choudhuryCells, features = overlap_pathways_up$gene_list[i],
                                   name = overlap_pathways_up$Description[i])
}

for (i in 1:nrow(overlap_pathways_down)) {
  choudhuryCells <- AddModuleScore(choudhuryCells, features = overlap_pathways_down$gene_list[i],
                                   name = overlap_pathways_down$Description[i],
                                   search = T)
}

# various plots of data
DotPlot(choudhuryCells, features = paste0(overlap_pathways_down$Description, "1"))
DotPlot(choudhuryCells, features = paste0(overlap_pathways_up$Description, "1"), split.by = "Condition")

FeaturePlot(choudhuryCells, features = paste0(overlap_pathways_down$Description, "1")[1:6], label = T,
            split.by = "Condition", min.cutoff = "q5")
VlnPlot(choudhuryCells, features = paste0(overlap_pathways_down$Description, "1")[1:6],
                    split.by = "Condition")
