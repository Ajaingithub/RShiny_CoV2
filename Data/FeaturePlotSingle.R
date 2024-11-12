FeaturePlotSingle <- function(obj, feature, metadata_column, ...){
  all_cells<- colnames(obj)
  groups<- unique(obj@meta.data[,"orig.ident"])
  # the minimal and maximal of the value to make the legend scale the same. 
  minimal<- min(obj[['RNA']]@data[feature, ])
  maximal<- max(obj[['RNA']]@data[feature, ])
  ps<- list()
  for (group in groups) {
    subset_indx<- obj@meta.data[, metadata_column] == group
    subset_cells<- all_cells[subset_indx]
    p<- FeaturePlot(obj, features = feature, cells= subset_cells, reduction = "wnn.umap", ...) +
      scale_colour_gradientn(
        # colours = c("lightgrey", "blue"),
        colours = c("#3361A5","#248AF3","#14B3FF","#88CEEF","#C1D5DC","#EAD397","#FDB31A","#E42A2A","#A31D1D"),
        limits=c(minimal,maximal),
        oob = scales::squish) +
      ggtitle(group) +
      theme(plot.title = element_text(size = 10, face = "bold"))
    ps[[group]]<- p
  }
  return(ps)
}