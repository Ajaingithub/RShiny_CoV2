FeaturePlotSingleADT <- function(obj, feature, metadata_column, ...){
  all_cells<- colnames(obj)
  groups<- unique(obj@meta.data[,"orig.ident"])
  # the minimal and maximal of the value to make the legend scale the same. 
  minimal<- min(obj[['IADT']]@data[feature, ])
  maximal<- max(obj[['IADT']]@data[feature, ])
  ps<- list()
  for (group in groups) {
    subset_indx<- obj@meta.data[, metadata_column] == group
    subset_cells<- all_cells[subset_indx]
    p<- FeaturePlot(obj, features = feature, cells= subset_cells, reduction = "wnn.umap", ...) +
      scale_colour_gradientn(
        # colours = c("lightgrey", "darkgreen"),
        colours = c("#000436","#021EA9","#1632FB","#6E34FC","#C732D5","#FD619D","#FF9965","#FFD32B","#FFFC5A"),
        limits=c(minimal,maximal),
        oob = scales::squish) +
      ggtitle(group) +
      theme(plot.title = element_text(size = 10, face = "bold"))
    ps[[group]]<- p
  }
  return(ps)
}