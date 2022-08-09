library(Seurat)
library(patchwork)

# load samples
external_data <- load("~/cardiac_seq_analysis/data/smoking/nature_data.Robj")
choudhury_data <- load("~/cardiac_seq_analysis/data/smoking/choudhury_cells.Robj")

all_data <- lapply(c(external_data, choudhury_data), get)

# prep for integration
all_data <- lapply(X = all_data, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = all_data, nfeatures = 3000)
all_data <- PrepSCTIntegration(object.list = all_data, anchor.features = features)

# integrate
anchors <- FindIntegrationAnchors(object.list = all_data, normalization.method = "SCT",
                                         anchor.features = features)
cardiac_combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# save
save(cardiac_combined, file="~/cardiac_seq_analysis/data/smoking/cardiac_combined.Robj")
