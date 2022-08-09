library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(tibble)
library(sctransform)
library(purrr)
library(patchwork)
library(RColorBrewer)

te <- readRDS("~/cardiac_seq_analysis/data/TE/smoking_tes_clustered.rds")
te_list <- readRDS("~/cardiac_seq_analysis/data/TE/TE_list.rds")

sample_tes <- rownames(te)
sample_tes <- sample_tes[sample_tes %in% te_list$final_TEs]

cell_type_markers <- FindAllMarkers(te,
                                    features = sample_tes,
                                    only.pos = T)
saveRDS(cell_type_markers, "~/cardiac_seq_analysis/data/TE/cell_type_markers.rds")

markers_abridged <- cell_type_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

write.csv(markers_abridged, "~/cardiac_seq_analysis/data/TE/cell_type_markers.csv")

l1hs <- FeaturePlot(te, features = "L1HS", label = T, min.cutoff = "q50")
ggsave("~/cardiac_seq_analysis/images/TE/seurat/CM_marker_L1HS.png", plot = l1hs,
       width = 10, height = 10)

l1pa2 <- FeaturePlot(te, features = "L1PA2", label = T, min.cutoff = "q50")
ggsave("~/cardiac_seq_analysis/images/TE/seurat/CM_marker_L1PA2.png", plot = l1pa2,
       width = 10, height = 10)

FeaturePlot(te, features = "L1ME4a", label = T, min.cutoff = "q50")
FeaturePlot(te, features = "AluYc", label = T, min.cutoff = "q0")
