library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(tibble)
library(sctransform)
library(purrr)
library(clusterProfiler)

## CREDIT: https://github.com/alkoenig/Atlas_of_Human_Heart_Failure_Lavine/blob/main/DCM_Nuclei_Seurat.R

# set wd
setwd("~/cardiac_seq_analysis/samples/")

# gather data
dataFolders<-as.matrix(read.csv("ploidy_folders.txt", header=FALSE))
dataFolders<-as.list(dataFolders[,1])

sampleNames <-as.matrix(read.csv("ploidy_names.txt", header=FALSE))
sampleNames <-as.list(sampleNames[,1])

raw.data<-as.list(NULL)
for(i in 1:length(dataFolders)){
  raw.data[i]<-Read10X(data.dir=dataFolders[[i]]) 
}

objs<-as.list(NULL)
for(i in 1:length(raw.data)){
  objs[i] <- CreateSeuratObject(counts=raw.data[[i]], project = sampleNames[[i]])
  objs[[i]]$orig.ident<-sampleNames[[i]]                              
  objs[[i]]<-RenameCells(objs[[i]],add.cell.id = sampleNames[[i]])
}

rm(raw.data)

# merge objects
ploidy_cells <- reduce(objs, merge)
rm(objs)

# save Raw object
saveRDS(ploidy_cells, file="~/cardiac_seq_analysis/data/ploidy/raw_ploidy_cells.rds")

load("~/cardiac_seq_analysis/data/ploidy/raw_ploidy_cells.rds")

# QC
ploidy_cells <- PercentageFeatureSet(ploidy_cells, pattern = "^MT-", col.name = "percent.mito")
VlnPlot(ploidy_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

ploidy_cells <- subset(ploidy_cells, percent.mito<10)
ploidy_cells <- subset(ploidy_cells, nCount_RNA > 1000 & nCount_RNA < 10000)

## saving
saveRDS(ploidy_cells, file="~/cardiac_seq_analysis/data/ploidy/ploidy_cells.rds")

# normalize
ploidy_cells <- SCTransform(ploidy_cells, vars.to.regress = "percent.mito", verbose = FALSE, conserve.memory=TRUE)

## PCA
ploidy_cells <- RunPCA(ploidy_cells, npcs=80, verbose = TRUE)
ElbowPlot(object = ploidy_cells, ndims=80) 

## TSNE
ploidy_cells <- RunTSNE(object = ploidy_cells, dims = 1:80)

## UMAP
ploidy_cells <- RunUMAP(object = ploidy_cells, dims = 1:80, verbose=FALSE)

## clustering
ploidy_cells <- FindNeighbors(object = ploidy_cells, reduction='pca', dims = 1:80)
ploidy_cells <- FindClusters(object = ploidy_cells,
                               resolution = c(0.3,0.4,0.5,0.6,0.7,0.8))

Idents(ploidy_cells)<-ploidy_cells$SCT_snn_res.0.8
DimPlot(object = ploidy_cells, reduction = 'umap', pt.size=0.5, label = T, split.by = 'orig.ident')
ggsave("~/cardiac_seq_analysis/images/ploidy/UMAP_Ploidy_Cells_0.8.png", width=20, height=10)

# Define marker genes
marker_genes <- c("DCN", "VWF", "RYR2", "KCNJ8", "MYH11", "C1QC", "NRXN1",
                  "CCL21", "PLIN1")

#Add Names Meta
current.cluster.ids<-levels(ploidy_cells)                  

new.cluster.ids<-c('Fibroblasts',
                   'Endothelium',
                   'Cardiomyocytes',
                   'Endothelium',
                   'Pericytes',
                   'Fibroblasts',
                   'Endothelium',
                   'Cardiomyocytes',
                   'Cardiomyocytes',
                   'Cardiomyocytes',
                   'Fibroblasts',
                   'Smooth_Muscle',
                   'Endothelium',
                   'Fibroblasts',
                   'Myeloid',
                   'Cardiomyocytes',
                   'Neuron',
                   'Pericytes',
                   'Lymphatic',
                   'Adipocytes')

Idents(ploidy_cells) <- plyr::mapvalues(x = Idents(ploidy_cells), from = current.cluster.ids, to = new.cluster.ids)

# save after clustering
saveRDS(ploidy_cells, file="~/cardiac_seq_analysis/data/ploidy/ploidy_cells.rds")

# Dot plots
png("~/cardiac_seq_analysis/images/ploidy/ploidy_dotplot.png",
    width = 12*100, height = 8*100,
    res=100)
DotPlot(ploidy_cells, features = marker_genes) + RotatedAxis()
dev.off()

# Bar plot of cell composition by smoking status
png("~/cardiac_seq_analysis/images/ploidy/ploidy_barplot.png",
    width = 15*100, height = 8*100,
    res=100)
ggplot(cell_percents_df, aes(x=reorder(cell_type, -percent), y=percent, fill=type)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=paste0(round(percent, 1), "%")), position=position_dodge(width = 0.9), vjust=-0.25) +
  labs(title = "Heart cell composition by ploidy", x = "\nCell type", y = "% total", fill = "Condition") +
  scale_fill_brewer(palette="Paired") + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# find DEGs between smoking and control (all cell types)

# cell_idents <- Idents(ploidy_cells)
# saveRDS(cell_idents, file = "~/cardiac_seq_analysis/data/ploidy_cell_idents.rds")
# 
# Idents(ploidy_cells) <- cardiac_combined[["Condition"]]

ploidy_markers <- FindMarkers(ploidy_cells,
                              ident.1 = "5828-2n",
                              ident.2 = "5828-6n",
                               group.by = "orig.ident",
                               min.pct = 0.25,
                               test.use = "wilcox")

ploidy_markers <- tibble::rownames_to_column(ploidy_markers, "gene_symbol")
write.table(ploidy_markers, "~/cardiac_seq_analysis/data/ploidy/ploidy_DE_markers.csv", sep=",", quote=F, row.names=F)

# subset CMs
ploidy_cardio = subset(ploidy_cells, idents = "Cardiomyocytes")

DimPlot(ploidy_cardio, label = T, split.by = 'orig.ident')
DimPlot(ploidy_cardio, label = T, group.by = 'orig.ident')

ploidy_cardio_markers <- FindMarkers(ploidy_cardio,
                              ident.1 = "5828-2n",
                              ident.2 = "5828-6n",
                              group.by = "orig.ident",
                              min.pct = 0.25,
                              test.use = "wilcox")
ploidy_cardio_markers <- tibble::rownames_to_column(ploidy_cardio_markers, "symbol")
write.table(ploidy_cardio_markers, "~/cardiac_seq_analysis/data/ploidy/ploidy_cardio_DE_markers.csv", sep=",", quote=F, row.names=F)

## saving
saveRDS(ploidy_cells, file="~/cardiac_seq_analysis/data/ploidy/ploidy_cells.rds")

### wikipathways analysis w/ clusterprofiler
## download human genome database
library(org.Hs.eg.db)
hs <- org.Hs.eg.db

## enriched in 2n
# convert symbols to entrez id
marker_names <- rownames(ploidy_cardio_markers[ploidy_cardio_markers$avg_log2FC > 0,])
marker_names <- gsub("MT-", "", marker_names)
marker_names[2] <- "CYTB"
marker_names[8] <- "COX3"
marker_names[9] <- "COX2"
marker_names[11] <- "COX1"
marker_names <- select(hs, 
       keys = marker_names,
       columns = c("ENTREZID"),
       keytype = "SYMBOL")
marker_names <- marker_names$ENTREZID

# find enriched pathways and plot
wp = enrichWP(marker_names, organism = "Homo sapiens") 
dotplot(wp, showCategory=20)
ggsave("~/cardiac_seq_analysis/images/ploidy/2n_enriched_pathways.png", width = 10,
       height = 8)

## enriched in 6n
# convert symbols to entrez ids
marker_names <- rownames(ploidy_cardio_markers[ploidy_cardio_markers$avg_log2FC < 0,])
marker_names <- gsub("MT-", "", marker_names)
marker_names <- select(hs, 
                       keys = marker_names,
                       columns = c("ENTREZID"),
                       keytype = "SYMBOL")
marker_names <- na.omit(marker_names)
marker_names <- marker_names$ENTREZID

# find symbols and plot
wp = enrichWP(marker_names, organism = "Homo sapiens") 
dotplot(wp, showCategory=20)
ggsave("~/cardiac_seq_analysis/images/ploidy/6n_enriched_pathways.png", width = 10,
       height = 8)

#------------------------------------

### Ploidy analysis w/ all 2n + 6n samples

# gather data
dataFolders<-list("4638-2n-1", "4638-2n-2", "5828-2n",
                  "1940-6n", "5828-6n", "5919-6n")

sampleNames <-list("4638-2n-1", "4638-2n-2", "5828-2n",
                   "1940-6n", "5828-6n", "5919-6n")

raw.data<-as.list(NULL)
for(i in 1:length(dataFolders)){
  raw.data[i]<-Read10X(data.dir=dataFolders[[i]]) 
}

objs<-as.list(NULL)
for(i in 1:length(raw.data)){
  objs[i] <- CreateSeuratObject(counts=raw.data[[i]], project = sampleNames[[i]])
  objs[[i]]$orig.ident<-sampleNames[[i]]                              
  objs[[i]]<-RenameCells(objs[[i]],add.cell.id = sampleNames[[i]])
}

rm(raw.data)

# merge objects
ploidy_cells <- reduce(objs, merge)
rm(objs)

# QC
ploidy_cells <- PercentageFeatureSet(ploidy_cells, pattern = "^MT-", col.name = "percent.mito")
VlnPlot(subset_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


s_4638_1 <- WhichCells(ploidy_cells,
                       expression = nCount_RNA > 2000 & nCount_RNA < 8000 & percent.mito < 5 & orig.ident == "4638-2n-1")
s_4638_2 <- WhichCells(ploidy_cells,
                       expression = nCount_RNA < 10000 & percent.mito < 5 & orig.ident == "4638-2n-2")
s_5828_2 <- WhichCells(ploidy_cells,
                       expression = nCount_RNA < 6000 & percent.mito < 10 & orig.ident == "5828-2n")
s_5828_6 <- WhichCells(ploidy_cells,
                       expression = nCount_RNA < 10000 & percent.mito < 5 & orig.ident == "5828-6n")
s_1940 <- WhichCells(ploidy_cells,
                     expression = nCount_RNA < 15000 & percent.mito < 5 & orig.ident == "1940-6n")
s_5919 <- WhichCells(ploidy_cells,
                     expression = nCount_RNA < 20000 & percent.mito < 5 & orig.ident == "5919-6n")

subset_cells <- subset(ploidy_cells,
                       cells = intersect(WhichCells(ploidy_cells),
                                       c(s_4638_1, s_4638_2, s_5828_2,
                                         s_5828_6, s_1940, s_5919)))

## saving
ploidy_cells <- subset_cells
saveRDS(ploidy_cells, file="~/cardiac_seq_analysis/data/ploidy/ploidy_cells_all.rds")

# normalize
ploidy_cells <- SCTransform(ploidy_cells, vars.to.regress = "percent.mito", verbose = FALSE, conserve.memory=TRUE)

## PCA
ploidy_cells <- RunPCA(ploidy_cells, npcs=80, verbose = TRUE)
ElbowPlot(object = ploidy_cells, ndims=80) 

## TSNE
ploidy_cells <- RunTSNE(object = ploidy_cells, dims = 1:80)

## UMAP
ploidy_cells <- RunUMAP(object = ploidy_cells, dims = 1:80, verbose=FALSE)

## clustering
ploidy_cells <- FindNeighbors(object = ploidy_cells, reduction='pca', dims = 1:80)
ploidy_cells <- FindClusters(object = ploidy_cells,
                             resolution = c(0.3,0.4,0.5,0.6,0.7,0.8))

Idents(ploidy_cells)<-ploidy_cells$SCT_snn_res.0.5

# add ploidy metadata
ploidy_cells$Ploidy <- ifelse(ploidy_cells$orig.ident %in% c('4638-2n-1', 
                                                               '4638-2n-2', 
                                                               '5828-2n'), "2n", "6n")

# Define marker genes
marker_genes <- c("DCN", "VWF", "RYR2", "MYH11", "KCNJ8", "C1QC", "NRXN1", "CCL21", "PLIN1")

#Add Names Meta
current.cluster.ids<-levels(ploidy_cells)                  

new.cluster.ids<-c('Fibroblasts',
                   'Endothelium',
                   'Cardiomyocytes',
                   'Cardiomyocytes',
                   'Endothelium',
                   'Smooth_Muscle',
                   'Endothelium',
                   'Cardiomyocytes',
                   'Fibroblasts',
                   'Cardiomyocytes',
                   'Cardiomyocytes',
                   'Pericytes',
                   'Endothelium',
                   'Myeloid',
                   'Neuron',
                   'Cardiomyocytes',
                   'Pericytes',
                   'Lymphatic',
                   'Adipocytes',
                   'Endothelium')

Idents(ploidy_cells) <- plyr::mapvalues(x = Idents(ploidy_cells), from = current.cluster.ids, to = new.cluster.ids)

# plot 2n and 6n clusters
DimPlot(object = ploidy_cells, reduction = 'umap', pt.size=0.5, label = T,
        split.by = "Ploidy")
ggsave("~/cardiac_seq_analysis/images/ploidy/UMAP_Ploidy_Cells_All_0.5.png", width=20, height=10)

# save after clustering
saveRDS(ploidy_cells, file="~/cardiac_seq_analysis/data/ploidy/ploidy_cells_all.rds")

# calculate ploidy cell type distribution
diploid <- subset(ploidy_cells, Ploidy == "2n")
hexaploid <- subset(ploidy_cells, Ploidy == "6n")

diploid_cell_counts <- sort(summary(Idents(diploid)))
diploid_cell_percent <- diploid_cell_counts/sum(diploid_cell_counts) * 100

hexaploid_cell_counts <- sort(summary(Idents(hexaploid)))
hexaploid_cell_percent <- hexaploid_cell_counts/sum(hexaploid_cell_counts) * 100

diploid_percents_df <- data.frame(cell_type = names(diploid_cell_counts),
                                  ploidy = "2n",
                                  percent = diploid_cell_percent,
                                  row.names = NULL)
hexaploid_percents_df <- data.frame(cell_type = names(hexaploid_cell_counts),
                                  ploidy = "6n",
                                  percent = hexaploid_cell_percent,
                                  row.names = NULL)
cell_percents_df <- rbind(diploid_percents_df, hexaploid_percents_df)

# plot as paired bar graph
ggplot(cell_percents_df, aes(x=reorder(cell_type, -percent), y=percent, fill=ploidy)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=paste0(round(percent, 1), "%")), position=position_dodge(width = 0.9), vjust=-0.25) +
  labs(title = "Heart cell composition by ploidy", x = "\nCell type", y = "% total", fill = "Condition") +
  scale_fill_brewer(palette="Paired") + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("~/cardiac_seq_analysis/images/ploidy/heart_cell_composition_by_ploidy.png", width = 15,
       height = 10)

# Dot plots
png("~/cardiac_seq_analysis/images/ploidy/ploidy_dotplot_all.png",
    width = 12*100, height = 8*100,
    res=100)
DotPlot(ploidy_cells, features = marker_genes) + RotatedAxis()
dev.off()

# find DEGs between smoking and control (all cell types)

# cell_idents <- Idents(ploidy_cells)
# saveRDS(cell_idents, file = "~/cardiac_seq_analysis/data/ploidy_cell_idents.rds")
# 
# Idents(ploidy_cells) <- cardiac_combined[["Condition"]]

ploidy_markers <- FindMarkers(ploidy_cells,
                              ident.1 = "2n",
                              ident.2 = "6n",
                              group.by = "Ploidy",
                              min.pct = 0.25,
                              test.use = "wilcox")

ploidy_markers <- tibble::rownames_to_column(ploidy_markers, "symbol")
write.table(ploidy_markers, "~/cardiac_seq_analysis/data/ploidy/ploidy_DE_markers.csv", sep=",", quote=F, row.names=F)

# subset CMs
ploidy_cardio = subset(ploidy_cells, idents = "Cardiomyocytes")

DimPlot(ploidy_cardio, label = T, split.by = 'Ploidy', group.by = 'orig.ident')

ploidy_markers <- FindMarkers(ploidy_cardio,
                              ident.1 = "2n",
                              ident.2 = "6n",
                              group.by = "Ploidy",
                              min.pct = 0.25,
                              test.use = "wilcox")
ploidy_cardio_markers <- tibble::rownames_to_column(ploidy_markers, "symbol")
write.table(ploidy_cardio_markers, "~/cardiac_seq_analysis/data/ploidy/ploidy_all_DE_markers.csv", sep=",", quote=F, row.names=F)

### wikipathways analysis w/ clusterprofiler
## download human genome database
library(org.Hs.eg.db)
hs <- org.Hs.eg.db

## enriched in 2n
# convert symbols to entrez id
marker_names <- ploidy_cardio_markers[ploidy_cardio_markers$avg_log2FC > 0,]$symbol
marker_names <- gsub("MT-", "", marker_names)
marker_names[16] <- "CYTB"
marker_names[11] <- "COX3"
marker_names[8] <- "COX2"
marker_names[7] <- "COX1"
marker_names <- select(hs,
                       keys = marker_names,
                       columns = c("ENTREZID"),
                       keytype = "SYMBOL")
marker_names <- marker_names$ENTREZID

# find enriched pathways and plot
wp = enrichWP(marker_names, organism = "Homo sapiens") 
dotplot(wp, showCategory=20)
ggsave("~/cardiac_seq_analysis/images/ploidy/2n_all_enriched_pathways.png", width = 10,
       height = 10)

## enriched in 6n
# convert symbols to entrez ids
marker_names <- ploidy_cardio_markers[ploidy_cardio_markers$avg_log2FC < 0,]$symbol
marker_names <- select(hs, 
                       keys = marker_names,
                       columns = c("ENTREZID"),
                       keytype = "SYMBOL")
marker_names <- na.omit(marker_names)
marker_names <- marker_names$ENTREZID

# find symbols and plot
wp = enrichWP(marker_names, organism = "Homo sapiens") 
dotplot(wp, showCategory=20)
ggsave("~/cardiac_seq_analysis/images/ploidy/6n_all_enriched_pathways.png", width = 10,
       height = 8)

# check for overlap between single sample and all sample analyses
single_sample_markers <- read.table("~/cardiac_seq_analysis/data/ploidy/ploidy_cardio_DE_markers.csv", sep = ",",
                                    header = T)
all_sample_markers <- read.table("~/cardiac_seq_analysis/data/ploidy/ploidy_all_DE_markers.csv", sep = ",",
                                    header = T)
head(single_sample_markers)
head(all_sample_markers)

# p-vals and log2fc are from all_samples
overlap_markers <- all_sample_markers[all_sample_markers$symbol %in%
                                        single_sample_markers$symbol,]
write.table(overlap_markers, "~/cardiac_seq_analysis/data/ploidy/ploidy_overlap_DE_markers.csv", sep = ",",
            quote=F, row.names=F)

# pathway analysis on overlap data
marker_names <- overlap_markers$symbol
marker_names <- gsub("MT-", "", marker_names)
marker_names[15] <- "CYTB"
marker_names[10] <- "COX3"
marker_names[7] <- "COX2"
marker_names[6] <- "COX1"
marker_names <- select(hs,
                       keys = marker_names,
                       columns = c("ENTREZID"),
                       keytype = "SYMBOL")
marker_names <- na.omit(marker_names)
marker_names <- marker_names$ENTREZID

wp = enrichWP(marker_names, organism = "Homo sapiens") 
dotplot(wp, showCategory=20)
ggsave("~/cardiac_seq_analysis/images/ploidy/ploidy_overlap_pathways.png", width = 10,
       height = 8)
