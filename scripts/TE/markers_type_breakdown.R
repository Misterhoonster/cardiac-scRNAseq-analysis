library(Seurat)
library(ggplot2)
library(plotly)
library(RColorBrewer)
library(patchwork)
library(scCustomize)

# load TE list
te_list <- read.table("~/cardiac_seq_analysis/data/TE/TE_list.csv", sep = ",",
                      header = T)[2:3]

# load all seurat markers
setwd("~/cardiac_seq_analysis/data/TE/results/")

y <- read.csv("young_age_te_markers.csv", header = T, row.names = "X")
y <- rownames(y)

m <- read.csv("middle_age_te_markers.csv", header = T, row.names = "X")
m <- rownames(m)

o <- read.csv("old_age_te_markers.csv", header = T, row.names = "X")
o <- rownames(o)

smoking_up <- read.csv("smoking_te_markers.csv", header = T, row.names = "X") %>%
  filter(avg_log2FC > 0)
smoking_up <- rownames(smoking_up)

smoking_down <- read.csv("smoking_te_markers.csv", header = T, row.names = "X") %>%
  filter(avg_log2FC < 0)
smoking_down <- rownames(smoking_down)

plots = list()

plot_pie <- function(markers, plot_title, filename) {
  counts <- te_list %>% filter(final_TEs %in% markers)
  counts <- counts %>% group_by(type) %>% count() %>% as.data.frame()
  
  custom <- list(
    l = 50,
    r = 50,
    b = 100,
    t = 50,
    pad = 4
  )
  
  fig <- plot_ly(counts, labels = ~type, values = ~n, type = 'pie', textinfo = 'label+percent')
  fig <- fig %>% layout(title = paste("TE types in", plot_title, "markers"), margin = custom,
                                        xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
  print(fig)
  Sys.sleep(10)
}

# plot pie chart for all conditions
conditions <- list(list(smoking_up, "smoking up", "smoking_up_markers_pie.png"),
                   list(smoking_down, "smoking down", "smoking_down_markers_pie.png"),
                   list(y, "young", "smoking_young_markers_pie.png"),
                   list(m, "middle", "smoking_middle_markers_pie.png"),
                   list(o, "old", "smoking_old_markers_pie.png"))

for (condition in conditions) {
  plot_pie(condition[[1]], condition[[2]], condition[[3]])
}


# plot pie chart for all markers
counts <- te_list %>% count(type)

fig <- plot_ly(counts, labels = ~type, values = ~n, type = 'pie', textinfo = 'label+percent')
fig <- fig %>% layout(title = paste("TE types in", "all", "markers"), margin = custom,
                      xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                      yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
print(fig)

plot_pie(te_list, "all", "all_markers_pie.png")

# plot features for each condition
conditions <- list(list(smoking_up, "smoking.up", "smoking_up_markers_features.png", "Up-reg in smoking"),
                   list(smoking_down, "smoking.down", "smoking_down_markers_features.png", "Down-reg in smoking"),
                   list(y, "young", "smoking_young_markers_pie_features.png", "Up-reg in young"),
                   list(m, "middle", "smoking_middle_markers_pie_features.png", "Up-reg in middle"),
                   list(o, "old", "smoking_old_markers_pie_features.png", "Up-reg in old"))


plot_features <- function(markers, plot_title, file_name, y_label) {
  print(paste("plotting", plot_title))
  scores <- AddModuleScore(te, features = markers, name = plot_title)
  
  if (plot_title == "smoking.up" | plot_title == "smoking.down") {
    plot <- FeaturePlot(scores, features = paste0(plot_title, "1"), label = T,
                        split.by = "Condition", min.cutoff = "q10") & 
      theme_minimal()
      
    ggsave(file_name, plot, width = 10, height = 5)
  }
  
  else {
    plot <- FeaturePlot(scores, features = paste0(plot_title, "1"), label = T,
                        split.by = "Age", min.cutoff = "q10") & 
      theme_minimal()
    ggsave(file_name, plot, width = 15, height = 5)
  }
}

for (condition in conditions) {
  plot_features(condition[[1]], condition[[2]], condition[[3]], condition[[4]])
}







