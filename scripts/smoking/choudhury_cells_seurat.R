library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(tibble)
library(sctransform)
library(purrr)

## CREDIT: https://github.com/alkoenig/Atlas_of_Human_Heart_Failure_Lavine/blob/main/DCM_Nuclei_Seurat.R

setwd("~/cardiac_seq_analysis/samples/")

# gather data
dataFolders<-as.matrix(read.csv("smoking_folders.txt", header=FALSE))
dataFolders<-as.list(dataFolders[,1])

sampleNames <-as.matrix(read.csv("smoking_names.txt", header=FALSE))
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
choudhuryCells<-reduce(objs, merge)
rm(objs)

# save Raw object
save(choudhuryCells, file="~/cardiac_seq_analysis/data/smoking/raw_choudhury_cells.Robj")

load("~/cardiac_seq_analysis/data/smoking/raw_choudhury_cells.Robj")

# QC
choudhuryCells <- PercentageFeatureSet(choudhuryCells, pattern = "^MT-", col.name = "percent.mito")
VlnPlot(choudhuryCells, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

choudhuryCells <- subset(choudhuryCells, percent.mito<5)
choudhuryCells <- subset(choudhuryCells, nCount_RNA > 1000 & nCount_RNA < 10000)

## saving
save(choudhuryCells, file="~/cardiac_seq_analysis/data/smoking/choudhury_cells.Robj")

### Choudhury cells seurat analysis

# load data
load("~/cardiac_seq_analysis/data/smoking/choudhury_cells.Robj")

# normalize
choudhuryCells <- SCTransform(choudhuryCells, vars.to.regress = "percent.mito", verbose = FALSE, conserve.memory=TRUE)

## PCA
choudhuryCells <- RunPCA(choudhuryCells, npcs=80, verbose = TRUE)
ElbowPlot(object = choudhuryCells, ndims=80) 

## TSNE
choudhuryCells <- RunTSNE(object = choudhuryCells, dims = 1:80)

## UMAP
choudhuryCells <- RunUMAP(object = choudhuryCells, dims = 1:80, verbose=FALSE)

## clustering
choudhuryCells <- FindNeighbors(object = choudhuryCells, reduction='pca', dims = 1:80)
choudhuryCells <- FindClusters(object = choudhuryCells,
                       resolution = c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
Idents(choudhuryCells)<-choudhuryCells$SCT_snn_res.0.3

# Define marker genes
marker_genes <- c("KCNJ8", "VWF", "DCN", "C1QC", "RYR2", "MYH11", "CD3E", "NRXN1",
                  "CCL21", "HAS1", "KIT", "PLIN1")

DotPlot(choudhuryCells, features = marker_genes) + RotatedAxis()

## Add Condition meta
choudhuryCells$Condition<- ifelse(choudhuryCells$orig.ident %in% c('1156',
                                                   '5111',
                                                   '5874'), "Smoking", "Control")

#Add Names Meta
current.cluster.ids<-levels(choudhuryCells)                  

new.cluster.ids<-c('Cardiomyocytes',
                   'Pericytes',
                   'Endothelium',
                   'Fibroblasts',
                   'Endothelium',
                   'Endothelium',
                   'Fibroblasts',
                   'Endothelium',
                   'Fibroblasts',
                   'Fibroblasts',
                   'Fibroblasts',
                   'Cardiomyocytes',
                   'Cardiomyocytes',
                   'Smooth Muscle',
                   'Cardiomyocytes',
                   'NK/T Cells',
                   'Cardiomyocytes',
                   'Myeloid',
                   'Neuron',
                   'Mast Cells',
                   'Lymphatic',
                   'B Cells',
                   'Neuron',
                   'Adipocytes')

Idents(choudhuryCells) <- plyr::mapvalues(x = Idents(choudhuryCells), from = current.cluster.ids, to = new.cluster.ids)

# Plot cell types
DimPlot(object = choudhuryCells, reduction = 'umap', label = T, split.by = "Condition")
ggsave("~/cardiac_seq_analysis/images/smoking/UMAP_Choudhury_Cells_0.3.png", width=12, height=6)

## finding and saving markers
markers <- FindAllMarkers(object = choudhuryCells,
                          only.pos = TRUE,
                          min.pct = 0.10,
                          thresh.use = 0.10)
write.table(markers, "~/cardiac_seq_analysis/data/smoking/choudhury_cell_type_markers_0.3.csv", sep=",", quote=F, row.names=F)

## saving
save(choudhuryCells, file="~/cardiac_seq_analysis/data/smoking/choudhury_cells.Robj")

# find smoking markers
smoking_markers <- FindMarkers(choudhuryCells,
                              ident.1 = "Smoking",
                              ident.2 = "Control",
                              group.by = "Condition",
                              min.pct = 0.25,
                              test.use = "wilcox")
smoking_markers <- tibble::rownames_to_column(smoking_markers, "symbol")
write.table(smoking_markers, "~/cardiac_seq_analysis/data/smoking/ch_DE_seurat.csv", sep=",", quote=F, row.names=F)
saveRDS(smoking_markers, '~/cardiac_seq_analysis/data/smoking/ch_DE_seurat.rds')

### wikipathways analysis w/ clusterprofiler
## download human genome database
library(org.Hs.eg.db)
hs <- org.Hs.eg.db

## enriched genes
# convert symbols to entrez id
marker_names <- smoking_markers$symbol
marker_names <- gsub("MT-", "", marker_names)
marker_names[114] <- "COX2"
marker_names[117] <- "COX3"
marker_names[121] <- "CYTB"
marker_names[155] <- "H1-10"
marker_names <- select(hs,
                       keys = marker_names,
                       columns = c("ENTREZID"),
                       keytype = "SYMBOL")
marker_names <- na.omit(marker_names)
marker_names <- marker_names$ENTREZID

# find enriched pathways and plot
wp = enrichWP(marker_names, organism = "Homo sapiens") 
dotplot(wp, showCategory=20)
ggsave("~/cardiac_seq_analysis/images/ploidy/2n_all_enriched_pathways.png", width = 10,
       height = 10)






