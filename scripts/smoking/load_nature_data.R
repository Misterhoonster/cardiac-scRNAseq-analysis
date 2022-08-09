library(Seurat)
library(dplyr)
library(ggplot2)
library(sctransform)
library(patchwork)
library(SeuratDisk)

# load integrated cell + nuclei data from Nature paper
load('~/cardiac_seq_analysis/data/smoking/Koenig_seq_data')

# function to get age, sex, condition of each sample
get_sample_info <- function(name) {
  age <- RefMerge[["Age_Group_Tertile"]]$Age_Group_Tertile[which(RefMerge[["orig.ident"]]$orig.ident == name, useNames = F)[1]]
  sex <- RefMerge[["Sex"]]$Sex[which(RefMerge[["orig.ident"]]$orig.ident == name, useNames = F)[1]]
  condition <- RefMerge[["condition"]]$condition[which(RefMerge[["orig.ident"]]$orig.ident == name, useNames = F)[1]]
  return(c(age, sex, condition))
}

# retrieve all unique samples
all_samples <- unique(RefMerge[["orig.ident"]]$orig.ident)

# get info on all samples
sample_data <- sapply(all_samples, get_sample_info)

# convert matrix to df
sample_data <- t(sample_data)
colnames(sample_data) <- c("Age", "Sex", "Condition")
sample_data <- as.data.frame(sample_data)

# filter for F, middle + old, donor samples
target_samples <- filter(sample_data, Age %in% c("Middle", "Old") & Sex == "Female"
       & Condition == "Donor") %>% arrange(Age)
target_samples

target_subset <- subset(RefMerge, subset = (orig.ident == "TWCM-11-42" |
                                              orig.ident == "TWCM-13-1" |
                                              orig.ident == "TWCM-13-101" |
                                              orig.ident == "TWCM-13-104" |
                                              orig.ident == "TWCM-13-132" |
                                              orig.ident == "TWCM-10-68" |
                                              orig.ident == "TWCM-11-104"))

trimmed_subset <- CreateSeuratObject(target_subset[["RNA"]]@counts,
                                     meta.data = target_subset[["orig.ident"]],
                                     min.cells = 3,
                                     min.features = 200)

# target_subset <- FindVariableFeatures(target_subset)
# target_subset <- RunPCA(target_subset, npcs=80, verbose = TRUE)
# ElbowPlot(object = target_subset, ndims=80)

save(trimmed_subset, file = "~/cardiac_seq_analysis/data/smoking/nature_data.Robj")
