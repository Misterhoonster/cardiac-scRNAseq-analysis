library(dplyr)
library(Seurat)
library(patchwork)

### 4638-2n
#Load the 4638 heart dataset
setwd('./samples')
seurat.data <- Read10X(data.dir = "./4638-2n-1/")
seurat <- CreateSeuratObject(counts = seurat.data, project = "4638_2n_1", min.cells = 3, min.features = 200) 

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# QC
choudhuryCells <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mito")
VlnPlot(choudhuryCells, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

choudhuryCells <- subset(choudhuryCells, percent.mito<5)
choudhuryCells <- subset(choudhuryCells, nCount_RNA > 1000 & nCount_RNA < 15000)

# normalize
choudhuryCells <- SCTransform(choudhuryCells, vars.to.regress = "percent.mito", verbose = FALSE, conserve.memory=TRUE)

## PCA
choudhuryCells <- RunPCA(choudhuryCells, npcs=80, verbose = TRUE)
ElbowPlot(object = choudhuryCells, ndims=80) 

## UMAP
choudhuryCells <- RunUMAP(object = choudhuryCells, dims = 1:80, verbose=FALSE)

## clustering
choudhuryCells <- FindNeighbors(object = choudhuryCells, reduction='pca', dims = 1:80)
choudhuryCells <- FindClusters(object = choudhuryCells,
                               resolution = c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
Idents(choudhuryCells)<-choudhuryCells$SCT_snn_res.0.9

DimPlot(choudhuryCells, label=T)

# Define marker genes
marker_genes <- c("KCNJ8", "VWF", "DCN", "C1QC", "RYR2", "MYH11", "CD3E", "NRXN1",
                  "CCL21", "HAS1", "KIT", "PLIN1")

DotPlot(choudhuryCells, features = marker_genes) + RotatedAxis()

#Add Names Meta
current.cluster.ids<-levels(choudhuryCells) 

new.cluster.ids<-c('CM',
                   'Endo',
                   'FB',
                   'FB',
                   'PC',
                   'Endo',
                   'SM',
                   'Endo',
                   'Endo',
                   'Myeloid',
                   'CM',
                   'Neuron',
                   'Adipo',
                   'Endo')

Idents(choudhuryCells) <- plyr::mapvalues(x = Idents(choudhuryCells), from = current.cluster.ids, to = new.cluster.ids)


DimPlot(choudhuryCells, reduction='umap', label=T, cols = c('CM' = 'coral', 'Endo' = 'cornflowerblue', 'FB' = 'chartreuse3',
                                                            'PC' = 'darkolivegreen1', 'SM' = 'lavender', 'Myeloid' = 'orange', 'Neuron' = 'gray', 'Adipo' = 'pink'))

### 5828-heart
#Load the 5828 heart dataset
seurat.data <- Read10X(data.dir = "./TE/5828-heart/")
seurat <- CreateSeuratObject(counts = seurat.data, project = "5828_heart", min.cells = 3, min.features = 200) 

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# QC
choudhuryCells <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mito")
VlnPlot(choudhuryCells, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

choudhuryCells <- subset(choudhuryCells, percent.mito<10)
choudhuryCells <- subset(choudhuryCells, nCount_RNA > 100 & nCount_RNA < 5000)

# normalize
choudhuryCells <- SCTransform(choudhuryCells, vars.to.regress = "percent.mito", verbose = FALSE, conserve.memory=TRUE)

## PCA
choudhuryCells <- RunPCA(choudhuryCells, npcs=80, verbose = TRUE)
ElbowPlot(object = choudhuryCells, ndims=80) 

## UMAP
choudhuryCells <- RunUMAP(object = choudhuryCells, dims = 1:80, verbose=FALSE)

## clustering
choudhuryCells <- FindNeighbors(object = choudhuryCells, reduction='pca', dims = 1:80)
choudhuryCells <- FindClusters(object = choudhuryCells,
                               resolution = c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
Idents(choudhuryCells)<-choudhuryCells$SCT_snn_res.0.9

DimPlot(choudhuryCells, label=T)

# Define marker genes
marker_genes <- c("KCNJ8", "VWF", "DCN", "C1QC", "RYR2", "MYH11", "CD3E", "NRXN1",
                  "CCL21", "HAS1", "KIT", "PLIN1")

DotPlot(choudhuryCells, features = marker_genes) + RotatedAxis()

#Add Names Meta
current.cluster.ids<-levels(choudhuryCells)                  

new.cluster.ids<-c('FB',
                   'Endo',
                   'FB',
                   'CM',
                   'FB',
                   'Endo',
                   'PC')

Idents(choudhuryCells) <- plyr::mapvalues(x = Idents(choudhuryCells), from = current.cluster.ids, to = new.cluster.ids)


DimPlot(choudhuryCells, reduction='umap', label=T, cols = c('CM' = 'coral', 'Endo' = 'cornflowerblue', 'FB' = 'chartreuse3',
                                                            'PC' = 'darkolivegreen1'))

### 5828-4n
#Load the 5828-4n dataset
seurat.data <- Read10X(data.dir = "./5828-4n/")
seurat <- CreateSeuratObject(counts = seurat.data, project = "5828_4n", min.cells = 3, min.features = 200) 

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# QC
choudhuryCells <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mito")
VlnPlot(choudhuryCells, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

choudhuryCells <- subset(choudhuryCells, percent.mito<5)
choudhuryCells <- subset(choudhuryCells, nCount_RNA > 100 & nCount_RNA < 5000)

# normalize
choudhuryCells <- SCTransform(choudhuryCells, vars.to.regress = "percent.mito", verbose = FALSE, conserve.memory=TRUE)

## PCA
choudhuryCells <- RunPCA(choudhuryCells, npcs=80, verbose = TRUE)
ElbowPlot(object = choudhuryCells, ndims=80) 

## UMAP
choudhuryCells <- RunUMAP(object = choudhuryCells, dims = 1:80, verbose=FALSE)

## clustering
choudhuryCells <- FindNeighbors(object = choudhuryCells, reduction='pca', dims = 1:80)
choudhuryCells <- FindClusters(object = choudhuryCells,
                               resolution = c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
Idents(choudhuryCells)<-choudhuryCells$SCT_snn_res.0.9

DimPlot(choudhuryCells, label = T)

# Define marker genes
marker_genes <- c("KCNJ8", "VWF", "DCN", "C1QC", "RYR2", "MYH11", "CD3E", "NRXN1",
                  "CCL21", "HAS1", "KIT", "PLIN1")

DotPlot(choudhuryCells, features = marker_genes) + RotatedAxis()

#Add Names Meta
current.cluster.ids<-levels(choudhuryCells)                  

new.cluster.ids<-c('CM',
                   'CM',
                   'CM',
                   'CM',
                   'CM',
                   'CM',
                   'Endo',
                   'FB')

Idents(choudhuryCells) <- plyr::mapvalues(x = Idents(choudhuryCells), from = current.cluster.ids, to = new.cluster.ids)

DimPlot(choudhuryCells, reduction='umap', label=T, cols = c('CM' = 'coral', 'Endo' = 'cornflowerblue', 'FB' = 'chartreuse3'))

table(Idents(choudhuryCells))

