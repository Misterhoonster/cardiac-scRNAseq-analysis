library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(tibble)
library(sctransform)
library(purrr)

## CREDIT: https://github.com/alkoenig/Atlas_of_Human_Heart_Failure_Lavine/blob/main/DCM_Nuclei_Seurat.R

# set wd
setwd("~/cardiac_seq_analysis/samples/")

# gather data
dataFolders <- list("5828-4n")

sampleNames <- list("5828-4n")

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
# save(ploidy_cells, file="raw_ploidy_cells.Robj")

# load("raw_ploidy_cells.Robj")

# QC
ploidy_cells <- PercentageFeatureSet(ploidy_cells, pattern = "^MT-", col.name = "percent.mito")
VlnPlot(ploidy_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

ploidy_cells <- subset(ploidy_cells, percent.mito<5)
ploidy_cells <- subset(ploidy_cells, nCount_RNA > 1000 & nCount_RNA < 10000)

## saving
# save(ploidy_cells, file="ploidy_cells.Robj")

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
ggsave("~/cardiac_seq_analysis/images/5828-4n_Umap.png", width=12, height=10)

# Define marker genes
marker_genes <- c("KCNJ8", "VWF", "DCN", "C1QC", "RYR2", "MYH11", "CD3E", "NRXN1",
                  "CCL21", "HAS1", "KIT", "PLIN1")

DotPlot(ploidy_cells, features = marker_genes) + RotatedAxis()
ggsave("~/cardiac_seq_analysis/images/5828-4n_Dot.png", width=12, height=10)
