library(Seurat)
library(dplyr)
library(ggplot2)
library(sctransform)
library(patchwork)
library(RColorBrewer)
library(tibble)
library(data.table)

# combined raw data
load("~/cardiac_seq_analysis/data/smoking/cardiac_combined.Robj")

## add metadata

# condition meta
cardiac_combined$Condition <- ifelse(cardiac_combined$orig.ident %in% c('1156',
                                                                   '5111',
                                                                   '5874'), "Smoking", "Control")

# age group meta
cardiac_combined$Age <- ifelse(cardiac_combined$orig.ident %in% c('1156',
                                                                           '5111',
                                                                           '936',
                                                                           'TWCM-11-42',
                                                                           'TWCM-13-1',
                                                                           'TWCM-13-101',
                                                                           'TWCM-13-104',
                                                                           'TWCM-13-132'),"Middle","Old")

# run PCA and UMAP
cardiac_combined <- RunPCA(cardiac_combined, npcs=80, verbose = T)
ElbowPlot(object = cardiac_combined, ndims=80) 

cardiac_combined <- RunUMAP(cardiac_combined, reduction = "pca", dims = 1:80)

# save
save(cardiac_combined, file = "~/cardiac_seq_analysis/smoking/cardiac_combined_PCA.Robj")

# clustering
cardiac_combined <- FindNeighbors(object = cardiac_combined, reduction='pca', dims = 1:80)
cardiac_combined <- FindClusters(object = cardiac_combined,
                                 force.recalc = T,
                                 save.SNN = T,
                               resolution = c(0.3,0.4,0.5,0.6,0.7,0.8))

# DE (finding and saving markers)
markers <- FindAllMarkers(object = cardiac_combined,
                          only.pos = TRUE,
                          min.pct = 0.10,
                          thresh.use = 0.10)
write.table(markers, "~/cardiac_seq_analysis/data/smoking/cardiac_combined_markers.csv", sep=",", quote=F, row.names=F)

# save
save(cardiac_combined, file = "~/cardiac_seq_analysis/data/smoking/cardiac_combined_DE.Robj")

## cell cluster identification
# CM: 4, 6, 8, 12
VlnPlot(cardiac_combined, features = c("RYR2"))
FeaturePlot(cardiac_combined, features = c("RYR2"), reduction = "umap",
            label = TRUE)
# SM: 15
VlnPlot(cardiac_combined, features = c("MYH11"))
FeaturePlot(cardiac_combined, features = c("MYH11"), reduction = "umap",
            label = TRUE)

# FB: 2, 5, 11, 13, 17
VlnPlot(cardiac_combined, features = c("DCN"))
FeaturePlot(cardiac_combined, features = c("DCN"), reduction = "umap",
            label = TRUE)
# Endo: 1, 7, 9, 10, 14, 16
VlnPlot(cardiac_combined, features = c("VWF"))
FeaturePlot(cardiac_combined, features = c("VWF"), reduction = "umap",
            label = TRUE)
# LEC: 23
VlnPlot(cardiac_combined, features = c("CCL21"))
FeaturePlot(cardiac_combined, features = c("CCL21"), reduction = "umap",
            label = TRUE)

# Neuron: 21, 28
VlnPlot(cardiac_combined, features = c("NRXN1"))
FeaturePlot(cardiac_combined, features = c("NRXN1"), reduction = "umap",
            label = TRUE)

# Epi: 24
VlnPlot(cardiac_combined, features = c("HAS1"))
FeaturePlot(cardiac_combined, features = c("HAS1"), reduction = "umap",
            label = TRUE)

# Mast: 25
VlnPlot(cardiac_combined, features = c("KIT"))
FeaturePlot(cardiac_combined, features = c("KIT"), reduction = "umap",
            label = TRUE)

# B cells: NA
VlnPlot(cardiac_combined, features = c("MZB1"))
FeaturePlot(cardiac_combined, features = c("MZB1"), reduction = "umap",
            label = TRUE)

# Adipo: 26
VlnPlot(cardiac_combined, features = c("PLIN1"))
FeaturePlot(cardiac_combined, features = c("PLIN1"), reduction = "umap",
            label = TRUE)

# Myeloid: 3, 27
VlnPlot(cardiac_combined, features = c("C1QC"))
FeaturePlot(cardiac_combined, features = c("C1QC"), reduction = "umap",
            label = TRUE)
# PC: 0, 19, 20, 22
VlnPlot(cardiac_combined, features = c("KCNJ8"))
FeaturePlot(cardiac_combined, features = c("KCNJ8"), reduction = "umap",
            label = TRUE)
# EC: 16
VlnPlot(cardiac_combined, features = c("LEPR"))
FeaturePlot(cardiac_combined, features = c("LEPR"), reduction = "umap",
            label = TRUE)

# NK/T cells: 18
VlnPlot(cardiac_combined, features = c("CD3E"))
FeaturePlot(cardiac_combined, features = c("CD3E"), reduction = "umap",
            label = TRUE)

# add cell type metadata
current.cluster.ids<-levels(cardiac_combined)               

new.cluster.ids<-c('Pericytes',
                   'Endothelium',
                   'Fibroblasts',
                   'Myeloid',
                   'Cardiomyocytes',
                   'Fibroblasts',
                   'Cardiomyocytes',
                   'Endothelium',
                   'Cardiomyocytes',
                   'Endothelium',
                   'Endothelium',
                   'Fibroblasts',
                   'Cardiomyocytes',
                   'Fibroblasts',
                   'Endothelium',
                   'Smooth_Muscle',
                   'Endothelium',
                   'Fibroblasts',
                   'T/NK_Cells',
                   'Pericytes',
                   'Pericytes',
                   'Neurons',
                   'Pericytes',
                   'Lymphatic',
                   'Epicardium',
                   'Mast_Cells',
                   'Adipocytes',
                   'Myeloid',
                   'Neurons')

Idents(cardiac_combined) <- plyr::mapvalues(x = Idents(cardiac_combined), from = current.cluster.ids, to = new.cluster.ids)

# save
save(cardiac_combined, file = "~/cardiac_seq_analysis/data/smoking/cardiac_combined_clustered.Robj")

# calculate smoking v non-smoking cell type distribution
control <- subset(cardiac_combined, Condition == "Control")
smoking <- subset(cardiac_combined, Condition == "Smoking")

control_cell_counts <- sort(summary(Idents(control)))
control_cell_percent <- control_cell_counts/sum(control_cell_counts) * 100

smoking_cell_counts <- sort(summary(Idents(smoking)))
smoking_cell_percent <- smoking_cell_counts/sum(smoking_cell_counts) * 100

control_percents_df <- data.frame(cell_type = names(control_cell_counts),
                                type = "Control",
                                percent = control_cell_percent,
                                row.names = NULL)
smoking_percents_df <- data.frame(cell_type = names(control_cell_counts),
                                  type = "Smoking",
                             percent = smoking_cell_percent,
                             row.names = NULL)
cell_percents_df <- rbind(control_percents_df, smoking_percents_df)

## Visualizations

# UMAP
png("~/cardiac_seq_analysis/images/smoking/UMAP_smoking.png",
    width = 10*100, height = 8*100,
    res=100)
DimPlot(object = smoking, reduction = 'umap', label = T, pt.size=0.1)
dev.off()

png("~/cardiac_seq_analysis/images/smoking/UMAP_control.png",
    width = 10*100, height = 8*100,
    res=100)
DimPlot(object = control, reduction = 'umap', label = T, pt.size=0.1)
dev.off()

# Violin plots
marker_genes <- c("KCNJ8", "VWF", "DCN", "C1QC", "RYR2", "MYH11", "CD3E", "NRXN1",
  "CCL21", "HAS1", "KIT", "PLIN1")

png("~/cardiac_seq_analysis/images/smoking/smoking_violin.png",
    width = 15*100, height = 8*100,
    res=100)
VlnPlot(cardiac_combined, features = marker_genes, pt.size = 0)
dev.off()

# Dot plots
png("~/cardiac_seq_analysis/images/smoking/smoking_dotplot.png",
    width = 12*100, height = 8*100,
    res=100)
DotPlot(cardiac_combined, features = marker_genes)
dev.off()

# Bar plot of cell composition by smoking status
png("~/cardiac_seq_analysis/images/smoking/smoking_status_barplot.png",
    width = 15*100, height = 8*100,
    res=100)
ggplot(cell_percents_df, aes(x=reorder(cell_type, -percent), y=percent, fill=type)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=paste0(round(percent, 1), "%")), position=position_dodge(width = 0.9), vjust=-0.25) +
  labs(title = "Heart cell composition by smoking status", x = "\nCell type", y = "% total", fill = "Condition") +
  scale_fill_brewer(palette="Paired") + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# find DEGs between smoking and control
cell_idents <- Idents(cardiac_combined)
saveRDS(cell_idents, file = "~/cardiac_seq_analysis/data/smoking/cell_idents.rds")

Idents(cardiac_combined) <- cardiac_combined[["Condition"]]

smoking_markers <- FindMarkers(cardiac_combined, ident.1 = "Smoking",
                               ident.2 = "Control", min.pct = 0.25,
                               test.use = "wilcox")

# change rownames to column and save to table
smoking_markers <- tibble::rownames_to_column(smoking_markers, "symbol")
write.table(smoking_markers, "~/cardiac_seq_analysis/data/smoking/smoking_DE_markers.csv", sep=",", quote=F, row.names=F)

### wikipathways analysis w/ clusterprofiler
## download human genome database
library(org.Hs.eg.db)
hs <- org.Hs.eg.db

smoking_markers <- read.table("~/cardiac_seq_analysis/data/smoking/smoking_DE_markers.csv", sep = ",", header = T)

## all DEGs
# convert symbols to entrez id
marker_names <- smoking_markers$gene_symbol
marker_names <- gsub("MT-", "", marker_names)
marker_names[59] = "LINC02899"
marker_names[191] <- "H1-10"
marker_names[199] <- "PRANCR"
marker_names[497] = "CYTB"
marker_names[516] = "COX2"
marker_names[536] = "COX3"
marker_names[560] = "COX1"
marker_names[577] = "MIR124-1HG"
marker_names[588] = "CDIN1"
marker_names[636] <- "TAMALIN"
marker_names <- select(hs, 
                       keys = marker_names,
                       columns = c("ENTREZID"),
                       keytype = "SYMBOL")
marker_names <- na.omit(marker_names)
marker_names <- marker_names$ENTREZID

# find pathways and plot
wp = enrichWP(marker_names, organism = "Homo sapiens") 
dotplot(wp, showCategory=20)
ggsave("~/cardiac_seq_analysis/images/smoking/smoking_all_pathways.png", width = 10,
       height = 8)

## enriched in smoking
# convert symbols to entrez id
marker_names <- smoking_markers[smoking_markers$avg_log2FC > 0,]$gene_symbol
marker_names[29] <- "H1-10"
marker_names[30] <- "PRANCR"
marker_names[111] <- "TAMALIN"
marker_names <- select(hs,
                       keys = marker_names,
                       columns = c("ENTREZID"),
                       keytype = "SYMBOL")
marker_names <- marker_names$ENTREZID

# find enriched pathways and plot
wp = enrichWP(marker_names, organism = "Homo sapiens") 
dotplot(wp, showCategory=20)
ggsave("~/cardiac_seq_analysis/images/smoking/smoking_all_enriched_pathways.png", width = 10,
       height = 10)

## downregulated in smoking
# convert symbols to entrez ids
marker_names <- smoking_markers[smoking_markers$avg_log2FC < 0,]$gene_symbol
marker_names <- gsub("MT-", "", marker_names)
marker_names[55] = "LINC02899"
marker_names[412] = "CYTB"
marker_names[427] = "COX2"
marker_names[443] = "COX3"
marker_names[461] = "COX1"
marker_names[477] = "MIR124-1HG"
marker_names[488] = "CDIN1"
marker_names <- select(hs, 
                       keys = marker_names,
                       columns = c("ENTREZID"),
                       keytype = "SYMBOL")
marker_names <- na.omit(marker_names)
marker_names <- marker_names$ENTREZID

# find symbols and plot
wp = enrichWP(marker_names, organism = "Homo sapiens") 
dotplot(wp, showCategory=20)
ggsave("~/cardiac_seq_analysis/images/smoking/smoking_all_downregulated_pathways.png", width = 10,
       height = 8)

## export to DESeq2 format

# reduce size of counts matrix
genes.use <- names(which(rowSums(cardiac_combined[["RNA"]]@counts) >= 15))
cardiac_trimmed_counts <- cardiac_combined[["RNA"]]@counts[genes.use,]

# convert sparse counts matrix to df
counts.df <- cardiac_trimmed_counts %>% as.matrix %>% t %>% as.data.frame
counts.df <- tibble::rownames_to_column(counts.df, "cellnames")
clusterassignments <- data.frame(cardiac_combined[["orig.ident"]])
clusterassignments <- tibble::rownames_to_column(clusterassignments, "cellnames")
counts.df <- merge(clusterassignments, counts.df, by = "cellnames")

rownames(counts.df) <- counts.df$cellnames
counts.df$cellnames <- NULL

counts.df <- counts.df %>% group_by(orig.ident) %>%
  summarise(across(everything(), ~ sum(., is.na(.), 0)))

rownames(counts.df) <- counts.df$orig.ident
counts.df$orig.ident <- NULL
smoking_counts <- counts.df

save(smoking_counts, file = "~/cardiac_seq_analysis/data/smoking/deseq_rawcounts.Robj")


## create metadata df

# Age column
Age <- c("Middle", "Middle", "Old", "Middle", "Middle", "Middle", "Middle",
         "Middle", "Middle", "Old", "Old", "Old")

# Condition column
Condition <- c("Smoking", "Smoking", "Smoking", "Control", "Control", "Control",
               "Control", "Control", "Control", "Control", "Control", "Control")

# Data source
Source <- c("Ch", "Ch", "Ch", "Ch", "Na", "Na", "Na", "Na", "Na", "Ch", "Na", "Na")

smoking_metadata <- data.frame(Condition, Age, Source)

# add sample names
rownames(smoking_metadata) <- c("1156", "5111", "5874", "936", "TWCM-11-42",
                                "TWCM-13-1", "TWCM-13-101", "TWCM-13-104",
                                "TWCM-13-132", "5828", "TWCM-10-68", "TWCM-11-104")

save(smoking_metadata, file = "~/cardiac_seq_analysis/data/smoking/deseq_rawcounts.Robj")

