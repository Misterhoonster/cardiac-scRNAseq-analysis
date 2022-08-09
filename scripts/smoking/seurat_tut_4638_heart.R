library(dplyr)
library(Seurat)
library(patchwork)

#Load the 4638 heart dataset
seurat.data <- Read10X(data.dir = "./filtered_feature_bc_matrix/")
seurat <- CreateSeuratObject(counts = seurat.data, project = "4638_heart", min.cells = 3, min.features = 200) 

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Visualize relationship between QC metrics
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# filter cells with features > 5k & mt > 10%
seurat <- subset(seurat, subset = nFeature_RNA < 5000 & percent.mt < 10)

# normalize the data
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize",
                        scale.factor = 10000)

# subset highly variable features
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)

# top 10 most variable genes
top10 <- head(VariableFeatures(seurat), 10)

# plot variable features with and wo labels
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0,
                     ynudge = 0)
plot1 + plot2

# scale the data
all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)

# run PCA
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))

# Examine and visualize PCA results in a few different ways
print(seurat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(seurat, dims = 1:2, reduction = "pca")
DimPlot(seurat, reduction = "pca")
DimHeatmap(seurat, dims = 1:20, cells = 100, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
seurat <- JackStraw(seurat, num.replicate = 100)
seurat <- ScoreJackStraw(seurat, dims = 1:20)

JackStrawPlot(seurat, dims = 1:20)
ElbowPlot(seurat)

# perform cell clustering
seurat <- FindNeighbors(seurat, dims = 1:10)
seurat <- FindClusters(seurat, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(seurat), 5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
seurat <- RunUMAP(seurat, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(seurat, reduction = "umap")

# run tSNE
seurat <- RunTSNE(seurat, dims = 1:10)
DimPlot(seurat, reduction = "tsne")

# save seurat object
saveRDS(seurat, file = "./seurat_tut_4638_heart.rds")

# find all markers of cluster 1
cluster1.markers <- FindMarkers(seurat, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers)

# final all markers distinguishing cluster 3 from clusters 5 to 8
cluster3.markers <- FindMarkers(seurat, ident.1 = 3, ident.2 = c(5, 8), min.pct = 0.25)
head(cluster3.markers)

# find markers for every cluster compared to all remaining cells, report
# only the positive ones
seurat.markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# visualizing marker expression patterns
## cardiomyocyte
VlnPlot(seurat, features = c("MYH6", "MYL7", "TNNT2"))
FeaturePlot(seurat, features = c("MYH6", "MYL7", "TNNT2"), reduction = "tsne",
            label = TRUE)
RidgePlot(seurat, features = c("MYH6", "MYL7", "TNNT2"))

## neuronal progenitor
VlnPlot(seurat, features = c("NEGR1", "NCAM2"))
FeaturePlot(seurat, features = c("NEGR1", "NCAM2"), reduction = "tsne",
            label = TRUE)
RidgePlot(seurat, features = c("NEGR1", "NCAM2"))

## fibroblast
VlnPlot(seurat, features = c("DDR2", "COL4A4", "COL6A3"))
FeaturePlot(seurat, features = c("DDR2", "COL4A4", "COL6A3"), reduction = "tsne",
            label = TRUE)
FeaturePlot(seurat, features = c("DDR2", "COL4A4", "COL6A3"), reduction = "umap",
            label = TRUE)
RidgePlot(seurat, features = c("DDR2", "COL4A4", "COL6A3"))

## pericyte
VlnPlot(seurat, features = c("PDGFRB"))
FeaturePlot(seurat, features = c("PDGFRB"), reduction = "tsne",
            label = TRUE)
FeaturePlot(seurat, features = c("PDGFRB"), reduction = "umap",
            label = TRUE)
RidgePlot(seurat, features = c("PDGFRB"))

## endothelium
VlnPlot(seurat, features = c("CDH5", "CD9", "VWF"))
FeaturePlot(seurat, features = c("CDH5", "CD9", "VWF"), reduction = "tsne",
            label = TRUE)
FeaturePlot(seurat, features = c("CDH5", "CD9", "VWF"), reduction = "umap",
            label = TRUE)
RidgePlot(seurat, features = c("CDH5", "CD9", "VWF"))

## adepocyte - hard to determine - n/a
VlnPlot(seurat, features = c("MEST", "PPARG", "FABP4", "GPAM", "FASN"))
FeaturePlot(seurat, features = c("MEST", "PPARG", "FABP4", "GPAM", "FASN"), reduction = "tsne",
            label = TRUE)
FeaturePlot(seurat, features = c("MEST", "PPARG", "FABP4", "GPAM", "FASN"), reduction = "umap",
            label = TRUE)
RidgePlot(seurat, features = c("MEST", "PPARG", "FABP4", "GPAM", "FASN"))

##lymphocyte
VlnPlot(seurat, features = c("CD53", "CD2", "CD4", "B2M", "CD96", "CD53"))
FeaturePlot(seurat, features = c("MEST", "PPARG", "FABP4", "GPAM", "FASN"), reduction = "tsne",
            label = TRUE)
FeaturePlot(seurat, features = c("MEST", "PPARG", "FABP4", "GPAM", "FASN"), reduction = "umap",
            label = TRUE)
RidgePlot(seurat, features = c("MEST", "PPARG", "FABP4", "GPAM", "FASN"))

seurat.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(seurat, features = top10$gene) + NoLegend()
