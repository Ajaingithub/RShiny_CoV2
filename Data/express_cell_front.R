featureplot_front <- function(obj, feature, reduction = "umap", x = "UMAP_1", y = "UMAP_2", size = 0.2, color = "blue") {
  p <- FeaturePlot(obj, feature, reduction = reduction)
  feature <- gsub("-", ".", feature)
  p2 <- ggplot(p$data[order(p$data[feature][, 1]), ], aes_string(x = x, y = y, color = feature)) +
    geom_point(size = size) +
    scale_colour_gradientn(
      colours = c("lightgrey", color),
      limits = c(0, max(p$data[feature])),
      oob = scales::squish
    ) +
    ggtitle(paste(feature, "expressed front", sep = " ")) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 10),
      panel.background = element_rect(fill = "white", colour = "black"),
      legend.text = element_text(size = 10),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    )
  return(p2)
}
