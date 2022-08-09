library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(tibble)
library(sctransform)
library(purrr)
library(patchwork)
library(plotly)

### Creating GENES-ONLY object
setwd("~/cardiac_seq_analysis/samples/TE")

# gather data
dataFolders<-list("4638-2n-1", "4638-heart", "5828-heart", "5087-heart-5prime",
                  "936-heart-5prime", "5919-2n-all", "1156-heart-5prime", 
                  "5874-heart-5prime", "5111-heart-5prime", "604-heart-5prime")

sampleNames <-list("4638-2n-1", "4638-heart", "5828-heart", "5087-heart-5prime",
                   "936-heart-5prime", "5919-2n-all", "1156-heart-5prime", 
                   "5874-heart-5prime", "5111-heart-5prime", "604-heart-5prime")

raw.data<-as.list(NULL)
for(i in 1:length(dataFolders)){
  print(paste("reading", dataFolders[[i]]))
  raw.data[i]<-Read10X(data.dir=dataFolders[[i]]) 
}

objs<-as.list(NULL)
for(i in 1:length(raw.data)){
  print(paste("creating object for", dataFolders[[i]]))
  objs[i] <- CreateSeuratObject(counts=raw.data[[i]], project = sampleNames[[i]])
  objs[[i]]$orig.ident<-sampleNames[[i]]                              
  objs[[i]]<-RenameCells(objs[[i]],add.cell.id = sampleNames[[i]])
}

rm(raw.data)

# merge objects
raw_cells <- purrr::reduce(objs, merge)
rm(objs)

# save raw cells
saveRDS(raw_cells, file = "~/cardiac_seq_analysis/data/TE/smoking_genes_raw.rds")
genes_only <- readRDS("~/cardiac_seq_analysis/data/TE/smoking_genes_raw.rds")
rm(raw_cells)

# QC
genes_only <- PercentageFeatureSet(genes_only, pattern = "^MT-", col.name = "percent.mito")
VlnPlot(genes_only, features = c("nCount_RNA"), y.max = 25000,
        pt.size = 0) + geom_boxplot(width = 0.1, color = "white")
VlnPlot(genes_only, features = c("percent.mito"),
        pt.size = 0) + geom_boxplot(width = 0.1, color = "white")
VlnPlot(genes_only, features = c("nFeature_RNA"),
        pt.size = 0) + geom_boxplot(width = 0.1, color = "white")

s_4638_2n_1 <- WhichCells(genes_only,
                       expression = nCount_RNA < 20000 & nFeature_RNA < 6000 & 
                         percent.mito < 5 &
                         orig.ident == "4638-2n-1")
s_4638_heart <- WhichCells(genes_only,
                       expression = nCount_RNA < 14000 & percent.mito < 5 &
                         nFeature_RNA < 5000 &
                         orig.ident == "4638-heart")
s_5919_2n_all <- WhichCells(genes_only,
                       expression = nCount_RNA < 9000 & percent.mito < 5 &
                         nFeature_RNA < 4500 &
                         orig.ident == "5919-2n-all")
s_5828_heart <- WhichCells(genes_only,
                        expression = nCount_RNA < 5000 & percent.mito < 40 &
                          nFeature_RNA < 3000 &
                          orig.ident == "5828-heart")
s_rest <- WhichCells(genes_only,
                     expression = nCount_RNA < 5000 & nFeature_RNA < 3000 & 
                       percent.mito < 5 &
                       orig.ident != "4638-2n-1" & orig.ident != "4638-heart" &
                       orig.ident != "5919-2n-all" & orig.ident != "5828-heart")

subset_cells <- subset(genes_only,
                       cells = intersect(WhichCells(genes_only),
                                         c(s_4638_2n_1, s_4638_heart,
                                           s_5919_2n_all, s_5828_heart,
                                           s_rest)))
genes_only <- subset_cells
rm(subset_cells)

# save QC cells
saveRDS(genes_only, file = "~/cardiac_seq_analysis/data/TE/smoking_genes_QC.rds")
genes_only <- readRDS("~/cardiac_seq_analysis/data/TE/smoking_genes_QC.rds")

# normalize
genes_only <- SCTransform(genes_only, vars.to.regress = "percent.mito",
                          verbose = F, conserve.memory = T)

## PCA
genes_only <- RunPCA(genes_only, npcs=80, verbose = T)
ElbowPlot(object = genes_only, ndims=80) 

## UMAP
genes_only <- RunUMAP(object = genes_only, dims = 1:80, verbose=F)

## clustering
genes_only <- FindNeighbors(object = genes_only, reduction='pca', dims = 1:80)
genes_only <- FindClusters(object = genes_only,
                             resolution = c(0.3))


Idents(genes_only) <- genes_only$SCT_snn_res.0.3

DimPlot(genes_only, reduction = "umap", label = T)

# add metadata
genes_only$Condition <- ifelse(genes_only$orig.ident %in% c('1156-heart-5prime', 
                                                            '5874-heart-5prime',
                                                            '5111-heart-5prime',
                                                            '604-heart-5prime'), "Smoking", "Control")
genes_only$Sex <- ifelse(genes_only$orig.ident %in% c('5087-heart-5prime',
                                                      '5919-2n-all',
                                                      '604-heart-5prime',
                                                      '5828-heart'), "Male", "Female")
genes_only$Age <- ifelse(genes_only$orig.ident %in% c('4638-2n-1',
                                                      '4638-heart'), "Young",
                         ifelse(genes_only$orig.ident %in% c('5087-heart-5prime',
                                                             '936-heart-5prime',
                                                             '1156-heart-5prime',
                                                             '5111-heart-5prime',
                                                             '604-heart-5prime'), "Middle", "Old"))

# Define marker genes
marker_genes <- c("KCNJ8", "VWF", "DCN", "C1QC", "RYR2", "MYH11", "CD3E", "NRXN1",
                  "CCL21", "HAS1", "KIT", "PLIN1")

DotPlot(genes_only, features = marker_genes) + RotatedAxis()

#Add Names Meta
current.cluster.ids<-levels(genes_only)                  

new.cluster.ids<-c('FB',
                   'Endo',
                   'PC',
                   'CM',
                   'Endo',
                   'CM',
                   'FB',
                   'FB',
                   'Endo',
                   'CM',
                   'T',
                   'SM',
                   'Endo',
                   'Myeloid',
                   'Neuron',
                   'PC',
                   'CM',
                   'CM',
                   'Mast',
                   'LEC',
                   'B',
                   'Adipo',
                   'Neuron')

Idents(genes_only) <- plyr::mapvalues(x = Idents(genes_only),
                                                 from = current.cluster.ids, to = new.cluster.ids)
# save genes only object
saveRDS(genes_only, file = "~/cardiac_seq_analysis/data/TE/smoking_genes_clustered.rds")
genes_only <- readRDS("~/cardiac_seq_analysis/data/TE/smoking_genes_clustered.rds")

# Find all markers
# DE (finding and saving markers)
markers <- FindAllMarkers(object = genes_only,
                          only.pos = T,
                          min.pct = 0.10,
                          thresh.use = 0.10)
write.csv(markers, "~/cardiac_seq_analysis/data/TE/genes_only_markers.csv")

### TE analysis
# set wd
setwd("~/cardiac_seq_analysis/samples/scTE_outs/")

# gather data
dataFolders<-list("4638-2n-1", "4638-heart", "5828-heart", "5087-heart-5prime",
                  "936-heart-5prime", "5919-2n-all", "1156-heart-5prime", 
                  "5874-heart-5prime", "5111-heart-5prime", "604-heart-5prime")

sampleNames <-list("4638-2n-1", "4638-heart", "5828-heart", "5087-heart-5prime",
                       "936-heart-5prime", "5919-2n-all", "1156-heart-5prime", 
                       "5874-heart-5prime", "5111-heart-5prime", "604-heart-5prime")

raw.data<-as.list(NULL)
for(i in 1:length(dataFolders)){
  print(paste("reading", dataFolders[[i]]))
  raw.data[[i]] <- read.table(paste0(dataFolders[[i]], ".csv"), sep = ",",
                            header = T)
  rownames(raw.data[[i]]) <- raw.data[[i]]$barcodes
  raw.data[[i]]$barcodes <- NULL
  raw.data[[i]] <- as.data.frame(t(raw.data[[i]]))
}


objs<-as.list(NULL)
for(i in 1:length(raw.data)){
  print(paste("creating object for", dataFolders[[i]]))
  objs[i] <- CreateSeuratObject(counts=raw.data[[i]], project = sampleNames[[i]])
  objs[[i]]$orig.ident<-sampleNames[[i]]                              
  objs[[i]]<-RenameCells(objs[[i]], add.cell.id = sampleNames[[i]])
}

rm(raw.data)

# merge objects
raw_cells <- reduce(objs, merge)
rm(objs)

# save raw cells
saveRDS(raw_cells, file = "~/cardiac_seq_analysis/data/TE/smoking_tes_raw.rds")
te <- readRDS("~/cardiac_seq_analysis/data/TE/smoking_tes_raw.rds")
rm(raw_cells)

## add sex + age metadata
# age column
Age <- c("Young", "Young", "Old", "Middle", "Middle", "Old", "Middle", "Old",
         "Middle", "Middle")

# sex column
Sex <- c("F", "F", "F", "M", "F", "M", "F", "F", "F", "M")

# add to metadata
te$Sex <- plyr::mapvalues(x = Idents(te),
                          from = unique(te$orig.ident), to = Sex)
te$Age <- plyr::mapvalues(x = Idents(te),
                          from = unique(te$orig.ident), to = Age)

# normalize
te <- PercentageFeatureSet(te, pattern = "^MT-", col.name = "percent.mito")
te <- SCTransform(te, vars.to.regress = c("percent.mito", "Sex", "Age"),
                          verbose = F, conserve.memory = T)

# rename cells to match genes_only
te <- RenameCells(te, new.names = paste0(Cells(te), "-1"))

# QC using genes only object
te <- subset(te, cells = intersect(WhichCells(te), WhichCells(genes_only)))

saveRDS(te, file = "~/cardiac_seq_analysis/data/TE/smoking_tes_qc.rds")
te <- readRDS("~/cardiac_seq_analysis/data/TE/smoking_tes_qc.rds")

# add metadata from genes only object
te <- AddMetaData(te, metadata = genes_only[[]] %>% select(SCT_snn_res.0.3, Condition))

## PCA
te <- RunPCA(te, npcs=80, verbose = T)
ElbowPlot(object = te, ndims=80) 

## UMAP
te <- RunUMAP(object = te, dims = 1:80, verbose=F)

DimPlot(te, label = T, group.by = "SCT_snn_res.0.3")

# Define marker genes
marker_genes <- c("KCNJ8", "VWF", "DCN", "C1QC", "RYR2", "MYH11", "CD3E", "NRXN1",
                  "CCL21", "HAS1", "KIT", "PLIN1")

DotPlot(te, features = marker_genes, group.by = "SCT_snn_res.0.3") + RotatedAxis()

# add cell type info
Idents(te) <- te$SCT_snn_res.0.3

current.cluster.ids<-levels(te)                

new.cluster.ids<-c('FB',
                   'Endo',
                   'PC',
                   'CM',
                   'Endo',
                   'CM',
                   'FB',
                   'FB',
                   'Endo',
                   'CM',
                   'T',
                   'SM',
                   'Endo',
                   'Myeloid',
                   'Neuron',
                   'PC',
                   'CM',
                   'CM',
                   'Mast',
                   'LEC',
                   'B',
                   'Adipo',
                   'Neuron')

Idents(te) <- plyr::mapvalues(x = Idents(te),
                                      from = current.cluster.ids, to = new.cluster.ids)
te$cell_type <- Idents(te)

# plot UMAP for genes-only and TE datasets
te_plot <- DimPlot(te, label = T) + labs(title = "Heart Cell Types (genes + TEs)") + theme_minimal()
genes_plot <- DimPlot(genes_only, label = T) + labs(title = "Heart Cell Types (only genes)") + theme_minimal()
patchwork::plot_layout(te_plot + genes_plot)
plot_grid(te_plot, genes_plot)
ggsave("~/cardiac_seq_analysis/images/TE/umap_te_genes.png",
       width = 20, height = 10)

# save TE seurat object
saveRDS(te, file = "~/cardiac_seq_analysis/data/TE/smoking_tes_clustered.rds")
te <- readRDS("~/cardiac_seq_analysis/data/TE/smoking_tes_clustered.rds")

# Confirm B cluster is B-cells
b_markers <- FindMarkers(object = te,
                       ident.1 = "B",
                          only.pos = T,
                          min.pct = 0.10,
                          thresh.use = 0.10)

# load and save list of TEs
te_list <- read.table("~/cardiac_seq_analysis/data/TE/TE_list.csv", sep = ",",
           header = T)[2:3]
saveRDS(te_list, "~/cardiac_seq_analysis/data/TE/TE_list.rds")

## run FindMarkers using just TEs

# find TEs in sample data
sample_tes <- rownames(te)
sample_tes <- sample_tes[sample_tes %in% te_list$final_TEs]

## find markers
# smoking markers
smoking_markers <- FindMarkers(object = te,
                               assay = "SCT",
                       features = sample_tes,
                       ident.1 = "Smoking",
                       ident.2 = "Control",
                       group.by = "Condition")
write.csv(smoking_markers, file = "~/cardiac_seq_analysis/data/TE/smoking_te_markers.csv")

# age markers
# extract age samples
age_te <- subset(te, subset = orig.ident == "4638-2n-1" | orig.ident == "4638-heart" |
                   orig.ident == "5919-2n-all" | orig.ident == "5087-heart-5prime" | orig.ident == "936-heart-5prime" |
                   orig.ident == "5828-heart")
saveRDS(age_te, file = "~/cardiac_seq_analysis/data/TE/age_te.rds")

young_age_markers <- FindMarkers(object = age_te,
                          features = sample_tes,
                          ident.1 = "Young",
                          group.by = "Age",
                          only.pos = T)
write.csv(young_age_markers, file = "~/cardiac_seq_analysis/data/TE/young_age_te_markers.csv")

middle_age_markers <- FindMarkers(object = age_te,
                                 features = sample_tes,
                                 ident.1 = "Middle",
                                 group.by = "Age",
                                 only.pos = T)
write.csv(middle_age_markers, file = "~/cardiac_seq_analysis/data/TE/middle_age_te_markers.csv")

old_age_markers <- FindMarkers(object = age_te,
                                 features = sample_tes,
                                ident.1 = "Old",
                                 group.by = "Age",
                                 only.pos = T)
write.csv(old_age_markers, file = "~/cardiac_seq_analysis/data/TE/old_age_te_markers.csv")

## Plots of top smoking TE DEs

setwd("~/cardiac_seq_analysis/images/TE/seurat/")

# down-regulated in smoking
FeaturePlot(te, features = c("AluSc", "AluSp", "AluY"), label = T, split.by = "Condition",
            min.cutoff = "q50")

DotPlot(te, features = c("AluSc", "AluSp", "AluY", "AluSc8",
                         "AluSq2", "AluSx3"), col.min = -1)
ggsave("dot_smoking_te_down.png",
       width = 10, height = 10)

plots <- VlnPlot(te,
                 features = c("AluSc", "AluSp", "AluY", "AluSc8", "AluSq2", "AluSx3"),
                 group.by = "Condition",
                 pt.size = 0,
                 ncol = 2,
                 combine = FALSE)

for(i in 1:length(plots)) {
  plots[[i]] <- plots[[i]] + geom_boxplot(width = 0.1, fill = "white")
}

vln_down_list <- Reduce( `+`, plots) +
  patchwork::plot_layout( ncol = 2 )

ggsave("vln_smoking_te_down.png", plot = vln_down_list,
       width = 7, height = 10)

alu_sc <- VlnPlot(te, features = c("AluSc"), group.by = "Condition", pt.size = 0) +
  geom_boxplot(width = 0.1, fill = "white")
ggsave("vln_alu_sc.png", plot = alu_sc,
       width = 5, height = 5)

ridge_down <- RidgePlot(te, features = c("AluSc", "AluSp", "AluY", "AluSc8",
                                         "AluSq2", "AluSx3"), group.by = "Condition",
          stack = T)
ggsave("ridge_smoking_te_down.png",
       plot = ridge_down, width = 8, height = 2)

# up-regulated in smoking
plots <- VlnPlot(te,
                 features = c("AluYb9", "MER5B"),
                 group.by = "Condition",
                 pt.size = 0,
                 ncol = 2,
                 combine = FALSE)

for(i in 1:length(plots)) {
  plots[[i]] <- plots[[i]] + geom_boxplot(width = 0.1, fill = "white")
}

vln_up_list <- Reduce( `+`, plots) +
  patchwork::plot_layout( ncol = 2 )

ggsave("vln_smoking_te_up.png", plot = vln_up_list, width = 10, height = 5)

DotPlot(te, features = c("AluYb9", "MER5B"), col.min = -1)
ggsave("dot_smoking_te_up.png",
       width = 10, height = 10)

ridge_up <- RidgePlot(te, features = c("AluYb9", "MER5B"), group.by = "Condition",
                        stack = T)
ggsave("ridge_smoking_te_up.png",
       plot = ridge_up, width = 4, height = 2)

# up-regulated in Young
DotPlot(age_te, features = c("AluYb8", "AluSp", "AluSq", "AluSc8",
                             "AluSx4", "AluSc"), col.min = -1)
ggsave("dot_smoking_te_young.png",
       width = 10, height = 10)

plots <- VlnPlot(age_te,
                 cols = c("#D7EAF3", "#77B5D9", "#14397D"),
                 features = c("AluYb8", "AluSp", "AluSq", "AluSc8",
                              "AluSx4", "AluSc"),
                 group.by = "Age",
                 pt.size = 0, 
                 ncol = 2,
                 combine = FALSE)

for(i in 1:length(plots)) {
  plots[[i]] <- plots[[i]] + geom_boxplot(width = 0.1, fill = "white")
}

vln_young_list <- Reduce( `+`, plots) +
  patchwork::plot_layout( ncol = 2 )

ggsave("vln_smoking_te_young.png", plot = vln_young_list,
       width = 7, height = 10)

vln_alu_yb8 <- VlnPlot(age_te, features = c("AluYb8"), group.by = "Age", pt.size = 0) +
  geom_boxplot(width = 0.1, fill = "white")
ggsave("vln_alu_yb8_age.png", plot = vln_alu_yb8,
       width = 7, height = 10)

ridge_young <- RidgePlot(age_te, features = c("AluYb8", "AluSp", "AluSq", "AluSc8",
                                          "AluSx4", "AluSc"), group.by = "Age",
                        stack = T)
ggsave("ridge_smoking_te_young.png",
       plot = ridge_young, width = 8, height = 2)

# up-regulated in Middle
DotPlot(age_te, features = c("MER5B", "AluYb9"), col.min = -1)
ggsave("dot_smoking_te_middle.png",
       width = 10, height = 10)

plots <- VlnPlot(age_te,
                 cols = c("#D7EAF3", "#77B5D9", "#14397D"),
                 features = c("MER5B", "AluYb9"),
                 group.by = "Age",
                 pt.size = 0,
                 ncol = 2,
                 combine = FALSE)

for(i in 1:length(plots)) {
  plots[[i]] <- plots[[i]] + geom_boxplot(width = 0.1, fill = "white")
}

vln_middle_list <- Reduce( `+`, plots) +
  patchwork::plot_layout( ncol = 2 )

ggsave("vln_smoking_te_middle.png", plot = vln_middle_list, width = 10, height = 5)

ridge_middle <- RidgePlot(age_te, features = c("MER5B", "AluYb9"),
                         group.by = "Age",
                         stack = T)
ggsave("ridge_smoking_te_middle.png",
       plot = ridge_middle, width = 8, height = 2)

# up-regulated in Old
DotPlot(age_te, features = c("Tigger9a", "AluSx1", "AluY", "AluSx",
                             "AluSz", "AluJb"), col.min = -1)
ggsave("dot_smoking_te_old.png",
       width = 10, height = 10)

plots <- VlnPlot(age_te,
                 cols = c("#D7EAF3", "#77B5D9", "#14397D"),
                 features = c("Tigger9a", "AluSx1", "AluY", "AluSx",
                              "AluSz", "AluJb"),
                 group.by = "Age",
                 pt.size = 0, 
                 ncol = 2,
                 combine = FALSE)

for(i in 1:length(plots)) {
  plots[[i]] <- plots[[i]] + geom_boxplot(width = 0.1, fill = "white")
}

vln_old_list <- Reduce( `+`, plots) +
  patchwork::plot_layout( ncol = 2 )

ggsave("vln_smoking_te_old.png", plot = vln_old_list,
       width = 7, height = 10)

ridge_old <- RidgePlot(age_te,
                          features = c("Tigger9a", "AluSx1", "AluY", "AluSx",
                                               "AluSz", "AluJb"),
                          group.by = "Age",
                          stack = T)
ggsave("ridge_smoking_te_old.png",
       plot = ridge_old, width = 8, height = 2)

### TE type analysis

# change TE names to rownames
rownames(te_list) <- te_list$final_TEs
te_list$final_TEs <- NULL

## Smoking
te_smoking <- subset(te, subset = Condition == "Smoking")
te_control <- subset(te, subset = Condition == "Control")

# calculate counts for smoking (use normalized vals)
counts <- rowSums(te_smoking[["SCT"]])[rownames(te_smoking[["SCT"]]) %in% rownames(te_list)]
smoking_counts <- as.data.frame(counts)
smoking_counts <- merge(te_list, smoking_counts, all=TRUE, by='row.names')
smoking_counts[is.na(smoking_counts$counts),]$counts = 0
smoking_counts <- smoking_counts %>% group_by(type) %>% summarise(counts = sum(counts)) %>% as.data.frame()
smoking_counts <- smoking_counts %>% mutate(percents = counts/sum(counts) * 100)
smoking_counts$Condition <- "Smoking"

# calculate counts for control (use normalized vals)
counts <- rowSums(te_control[["SCT"]])[rownames(te_control[["SCT"]]) %in% rownames(te_list)]
control_counts <- as.data.frame(counts)
control_counts <- merge(te_list, control_counts, all=TRUE, by='row.names')
control_counts[is.na(control_counts$counts),]$counts = 0
control_counts <- control_counts %>% group_by(type) %>% summarise(counts = sum(counts)) %>% as.data.frame()
control_counts <- control_counts %>% mutate(percents = counts/sum(counts) * 100)
control_counts$Condition <- "Control"

# combine dfs
combined_counts <- rbind(smoking_counts, control_counts)


p <- ggplot(data=combined_counts, aes(x=reorder(type, -percents), y=percents, fill=Condition)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal() + scale_fill_brewer(palette='Reds') +
  labs(title = "TEs by type (FB)", x = "Type", y = "% total")


# graph pie charts
m <- list(
  l = 50,
  r = 50,
  b = 100,
  t = 100,
  pad = 4
)

fig_smoking <- plot_ly(smoking_counts, labels = ~type, values = ~freq, type = 'pie', textinfo = 'label+percent')
fig_smoking <- fig_smoking %>% layout(title = 'TE Type Breakdown in Smokers', margin = m,
                      xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                      yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
fig_smoking

fig_control <- plot_ly(control_counts, labels = ~type, values = ~counts, type = 'pie', textinfo = 'label+percent')
fig_control <- fig_control %>% layout(title = 'TE Type Breakdown in Non-Smokers', margin = m,
                      xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                      yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
fig_control


## Age
## Smoking
te_young <- subset(te, subset = Age == "Young")
te_middle <- subset(te, subset = Age == "Middle")
te_old <- subset(te, subset = Age == "Old")

# calculate counts for young (use normalized vals)
counts <- rowSums(te_young[["SCT"]])[rownames(te_young[["SCT"]]) %in% rownames(te_list)]
young_counts <- as.data.frame(counts)
young_counts <- merge(te_list, young_counts, all=TRUE, by='row.names')
young_counts[is.na(young_counts$counts),]$counts = 0

# calculate counts for middle (use normalized vals)
counts <- rowSums(te_middle[["SCT"]])[rownames(te_middle[["SCT"]]) %in% rownames(te_list)]
middle_counts <- as.data.frame(counts)
middle_counts <- merge(te_list, middle_counts, all=TRUE, by='row.names')
middle_counts[is.na(middle_counts$counts),]$counts = 0

# calculate counts for old (use normalized vals)
counts <- rowSums(te_old[["SCT"]])[rownames(te_old[["SCT"]]) %in% rownames(te_list)]
old_counts <- as.data.frame(counts)
old_counts <- merge(te_list, old_counts, all=TRUE, by='row.names')
old_counts[is.na(old_counts$counts),]$counts = 0

rm(te_young)
rm(te_middle)
rm(te_old)

# graph pie charts
fig_young <- plot_ly(young_counts, labels = ~type, values = ~counts, type = 'pie', textinfo = 'label+percent')
fig_young <- fig_young %>% layout(title = 'TE Type Breakdown in Young Adults', margin = m,
                                      xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                      yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
fig_young

fig_middle <- plot_ly(middle_counts, labels = ~type, values = ~counts, type = 'pie', textinfo = 'label+percent')
fig_middle <- fig_middle %>% layout(title = 'TE Type Breakdown in Middle-Aged Adults', margin = m,
                                      xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                      yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
fig_middle

fig_old <- plot_ly(old_counts, labels = ~type, values = ~counts, type = 'pie', textinfo = 'label+percent')
fig_old <- fig_old %>% layout(title = 'TE Type Breakdown in Old Adults', margin = m,
                                      xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                      yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
fig_old

# for each cell type
setwd("~/cardiac_seq_analysis/images/TE/seurat/celltypes/")

smoking_plot_pie <- function(cell) {
  ## Smoking
  te_smoking <- subset(te, subset = Condition == "Smoking" & cell_type == cell)
  te_control <- subset(te, subset = Condition == "Control" & cell_type == cell)
  
  print(levels(te_smoking))
  
  # calculate counts for smoking (use normalized vals)
  counts <- rowSums(te_smoking[["SCT"]])[rownames(te_smoking[["SCT"]]) %in% rownames(te_list)]
  smoking_counts <- as.data.frame(counts)
  smoking_counts <- merge(te_list, smoking_counts, all=TRUE, by='row.names')
  smoking_counts[is.na(smoking_counts$counts),]$counts = 0
  smoking_counts <- smoking_counts %>% group_by(type) %>% summarise(counts = sum(counts)) %>% as.data.frame()
  smoking_counts <- smoking_counts %>% mutate(percents = counts/sum(counts) * 100)
  smoking_counts$Condition <- "Smoking"
  
  # calculate counts for control (use normalized vals)
  counts <- rowSums(te_control[["SCT"]])[rownames(te_control[["SCT"]]) %in% rownames(te_list)]
  control_counts <- as.data.frame(counts)
  control_counts <- merge(te_list, control_counts, all=TRUE, by='row.names')
  control_counts[is.na(control_counts$counts),]$counts = 0
  control_counts <- control_counts %>% group_by(type) %>% summarise(counts = sum(counts)) %>% as.data.frame()
  control_counts <- control_counts %>% mutate(percents = counts/sum(counts) * 100)
  control_counts$Condition <- "Control"
  
  # combine dfs
  combined_counts <- rbind(smoking_counts, control_counts)
  print(combined_counts)
  
  p <- ggplot(data=combined_counts, aes(x=reorder(type, -percents), y=percents, fill=Condition)) +
    geom_bar(stat="identity", color="black", position=position_dodge())+
    theme_minimal() + scale_fill_brewer(palette='Reds') +
    labs(title = paste0("TEs by type (", cell_type, ")"), x = "Type", y = "% total")
  
  ggsave(paste0("smoking_", cell_type, ".png"),
         plot = p,
         width = 8,
         height = 5)
}

aging_plot_pie <- function(cell) {
  ## Smoking
  te_young <- subset(te, subset = Age == "Young" & cell_type == cell)
  te_middle <- subset(te, subset = Age == "Middle" & cell_type == cell)
  te_old <- subset(te, subset = Age == "Old" & cell_type == cell)
  
  # calculate counts for young (use normalized vals)
  counts <- rowSums(te_young[["SCT"]])[rownames(te_young[["SCT"]]) %in% rownames(te_list)]
  young_counts <- as.data.frame(counts)
  young_counts <- merge(te_list, young_counts, all=TRUE, by='row.names')
  young_counts[is.na(young_counts$counts),]$counts = 0
  young_counts <- young_counts %>% group_by(type) %>% summarise(counts = sum(counts)) %>% as.data.frame()
  young_counts <- young_counts %>% mutate(percents = counts/sum(counts) * 100)
  young_counts$Age <- "Young"
  
  # calculate counts for middle (use normalized vals)
  counts <- rowSums(te_control[["SCT"]])[rownames(te_middle[["SCT"]]) %in% rownames(te_list)]
  middle_counts <- as.data.frame(counts)
  middle_counts <- merge(te_list, middle_counts, all=TRUE, by='row.names')
  middle_counts[is.na(middle_counts$counts),]$counts = 0
  middle_counts <- middle_counts %>% group_by(type) %>% summarise(counts = sum(counts)) %>% as.data.frame()
  middle_counts <- middle_counts %>% mutate(percents = counts/sum(counts) * 100)
  middle_counts$Age <- "Middle"
  
  # calculate counts for middle (use normalized vals)
  counts <- rowSums(te_control[["SCT"]])[rownames(te_old[["SCT"]]) %in% rownames(te_list)]
  old_counts <- as.data.frame(counts)
  old_counts <- merge(te_list, old_counts, all=TRUE, by='row.names')
  old_counts[is.na(old_counts$counts),]$counts = 0
  old_counts <- old_counts %>% group_by(type) %>% summarise(counts = sum(counts)) %>% as.data.frame()
  old_counts <- old_counts %>% mutate(percents = counts/sum(counts) * 100)
  old_counts$Age <- "Old"
  
  # combine dfs
  combined_counts <- rbind(young_counts, middle_counts, old_counts)
  combined_counts$Age <- factor(combined_counts$Age, levels = c("Young", "Middle", "Old"))
  print(combined_counts)
  
  p <- ggplot(data=combined_counts, aes(x=reorder(type, -percents), y=percents, fill=Age)) +
    geom_bar(stat="identity", color="black", position=position_dodge())+
    theme_minimal() + scale_fill_brewer(palette='Greens') +
    labs(title = paste0("TEs by type (", cell_type, ")"), x = "Type", y = "% total")
  
  ggsave(paste0("age_", cell_type, ".png"),
         plot = p,
         width = 8,
         height = 5)
}


for (cell_type in levels(te)) {
  print(paste0("Plotting ", cell_type))
  smoking_plot_pie(cell_type)
}

for (cell_type in levels(te)) {
  print(paste0("Plotting ", cell_type))
  aging_plot_pie(cell_type)
}
