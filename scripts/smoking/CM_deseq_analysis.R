library(DESeq2)
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(data.table)
library(annotables)
library(EnhancedVolcano)
library(tibble)

# load choudhury cells Seurat object
load("~/cardiac_seq_analysis/data/smoking/choudhury_cells.Robj")

# DESeq pipeline
run_deseq <- function(celltype) {
  prepare_rawcounts(celltype)
  dds <- create_dds(celltype)
  plot_pca(dds, celltype)
  find_degs(dds, celltype)
}

# prepare rawcounts matrix for DESeq
prepare_rawcounts <- function(celltype) {
  print(paste("Preparing rawcounts for", celltype))
  
  ## export to DESeq2 format
  counts <- NULL
  if (celltype == "T cells") {
    counts <- subset(choudhuryCells, idents = "NK/T Cells")
  }
  else {
    counts <- subset(choudhuryCells, idents = celltype)
  }
  counts <- counts[["RNA"]]@counts
  
  # convert sparse counts matrix to df
  counts.df <- counts %>% as.matrix %>% t %>% as.data.frame
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
  
  counts_t <- as.data.frame(t(counts.df))
  
  saveRDS(counts_t, file = paste0("~/cardiac_seq_analysis/data/smoking/celltypes/",
                                  celltype, "/rawcounts.rds"))
}

create_dds <- function(celltype) {
  print(paste("creating DDS for", celltype))
  
  # load rawcounts and metadata
  rawcounts <- readRDS(paste0("~/cardiac_seq_analysis/data/smoking/celltypes/",
                              celltype, "/rawcounts.rds"))
  
  load('~/cardiac_seq_analysis/data/smoking/choudhury_metadata.Robj')
  
  # include only overlapping metadata
  overlap_idx <- which(rownames(smoking_metadata) %in% colnames(rawcounts))
  smoking_metadata <- smoking_metadata[rownames(smoking_metadata)[overlap_idx],]
  
  # re-order columns of raw counts df
  reorder_idx <- match(rownames(smoking_metadata), colnames(rawcounts))
  
  reorder_counts <- rawcounts[, reorder_idx]
  
  # create DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = reorder_counts,
                                   colData = smoking_metadata,
                                   design = ~ Age + Condition)
  return(dds)
}

plot_pca <- function(dds, celltype) {
  print(paste("plotting PCA for", celltype))
  
  tryCatch({
    # log transform counts
    vsd_all <- DESeq2::vst(dds, blind = T)
    
    # plot pca
    pca <- plotPCA(vsd_all, intgroup = "Age")
    ggsave(paste0("~/cardiac_seq_analysis/images/smoking/celltypes/",
                  celltype, "/pca_age.png"),
           plot = pca,
           width = 7, height = 5)
    
    pca <- plotPCA(vsd_all, intgroup = "Condition")
    ggsave(paste0("~/cardiac_seq_analysis/images/smoking/celltypes/",
                  celltype, "/pca_condition.png"),
           plot = pca,
           width = 7, height = 5)
  },
    error = function(e) {
      print(e)
      print(paste("PCA could not be performed on", celltype))
    }
  )
  
}

find_degs <- function(dds, celltype, alpha = 0.05, lfc = 0.1) {
  print(paste("finding DEGs for", celltype))
  
  tryCatch({
    # run deseq analysis
    dds <- DESeq(dds)
    
    # plot dispersion estimates
    # DESeq2::plotDispEsts(dds)
    # ggsave(paste0("~/cardiac_seq_analysis/images/smoking/celltypes/", celltype,
    #               "/disp_ests.png"),
    #        width = 10, height = 10)
    
    # calculate smoking results
    smoking_res <- results(dds,
                           contrast = c("Condition", "Smoking", "Control"),
                           alpha = alpha,
                           lfcThreshold = lfc)
    
    # lfc shrinkage
    smoking_res <- lfcShrink(dds,
                             coef = 3,
                             res = smoking_res)
    
    # turn symbols into column
    smoking_res_all <- data.frame(smoking_res) %>%
      rownames_to_column(var = "symbol")
    
    # select significant genes with padj < 0.05
    smoking_res_sig <- subset(smoking_res_all, padj < 0.05) %>%
      data.frame()
    
    # extract the top genes with padj values
    ch_top_de_genes <- smoking_res_sig %>%
      arrange(padj) %>%
      dplyr::select(symbol, padj, log2FoldChange)
    write.table(ch_top_de_genes,
                paste0("~/cardiac_seq_analysis/data/smoking/celltypes/",
                       celltype, "/deseq_markers.csv"),
                sep=",", quote=F, row.names=F)
    saveRDS(ch_top_de_genes, file = paste0("~/cardiac_seq_analysis/data/smoking/celltypes/",
                                           celltype, "/deseq_markers.rds"))
    
    # plot volcano
    plot_volcano(smoking_res_all, celltype)
    
    # plot heatmap
    plot_heatmap(dds, smoking_res_sig, celltype)
  },
  error = function(e) {
    print(e)
    print(paste("Error while running DESeq for", celltype))
  })
}

plot_volcano <- function (smoking_res_all, celltype) {
  print(paste("plotting volcano for", celltype))
  
  # add column showing if sig or not
  smoking_res_all <- smoking_res_all %>%
    mutate(threshold = padj < 0.05)
  
  volcano <- EnhancedVolcano(smoking_res_all,
                  lab = smoking_res_all$symbol,
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  pCutoff = 0.05)
  ggsave(paste0("~/cardiac_seq_analysis/images/smoking/celltypes/",
                celltype, "/volcano.png"),
         plot = volcano,
         width = 12, height = 10)
}

plot_heatmap <- function (dds, smoking_res_sig, celltype) {
  print(paste("plotting heatmap for", celltype))
  
  # load smoking metadata
  load('~/cardiac_seq_analysis/data/smoking/choudhury_metadata.Robj')
  
  tryCatch({
    # Determine the size factors to use for normalization
    dds <- estimateSizeFactors(dds)
    
    # Extract the normalized counts
    smoking_normalized_counts <- counts(dds, normalized=TRUE)
    
    # Subset normalized counts to significant genes
    sig_norm_counts <- smoking_normalized_counts[smoking_res_sig$symbol, ]
    
    # Choose heatmap color palette
    heat_colors <- brewer.pal(n = 6, name = "YlOrRd")
    
    # Plot heatmap
    heatmap <- pheatmap(sig_norm_counts,
                        color = heat_colors, 
                        cluster_rows = T, 
                        show_rownames = F,
                        annotation = dplyr::select(smoking_metadata, Condition), 
                        scale = "row")
    ggsave(paste0("~/cardiac_seq_analysis/images/smoking/celltypes/",
                  celltype, "/heatmap.png"),
           plot = heatmap,
           width = 12, height = 10)
  },
  error = function(e) {
    print(e)
    print("Heatmap could not be plotted")
  })
  
}

# run DESeq pipeline for all cell types
celltypes = c("Cardiomyocytes", "Fibroblasts", "Endothelium", "Myeloid",
              "Neuron", "Pericytes", "Smooth Muscle", "Mast Cells",
              "NK/T Cells", "Lymphatic", "B Cells", "Adipocytes")
celltypes = c("T cells", "Lymphatic", "B Cells", "Adipocytes")

for (celltype in celltypes) {
  run_deseq(celltype)
}

