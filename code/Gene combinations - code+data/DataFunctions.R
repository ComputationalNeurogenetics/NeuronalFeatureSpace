 

# Function takes 1 cluster and sets the cells to be in group A and all the rest as group B. Returns a dataframe gathered from these to be used in combiroc analysis.
cluster.dfs <- function(obj, cluster,top.markers_rna){
  
  clust <- top.markers_rna[which(top.markers_rna$group == cluster),]
  clust.features <- clust$feature_symbol
  #message_parallel(print(cluster))
  #Data from a cluster's cells and 25 best markers
  clust <- FetchData(object = obj, vars = clust.features, cells = colnames(obj)[Idents(obj) == cluster])
  clust <- clust[, order(colnames(clust))]
  
  #Creating cell_id and Class columns, arranging data into right format
  clust <- rownames_to_column(clust, var = "cell_id")
  clust$Class <- rep(c("A"), time = nrow(clust))
  clust <- dplyr::select(clust, cell_id, Class, 2:(ncol(clust)-1))
  
  #Getting all other cluster cells and assigning them as group B
  non.clust <- FetchData(object = obj, vars = clust.features, cells = colnames(obj)[Idents(obj) != cluster])
  non.clust <- non.clust[, order(colnames(non.clust))]
  
  #Creating cell_id and Class columns, arranging data into right format
  non.clust <- rownames_to_column(non.clust, var = "cell_id")
  non.clust$Class <- rep(c("B"), time = nrow(non.clust))
  non.clust <- dplyr::select(non.clust, cell_id, Class, 2:(ncol(non.clust)-1))
  
  #Combining group A and group B for further combiROC analysis
  data.clust <- rbind(clust, non.clust)
  return(data.clust)
  
}



#Function to create RNA_name assay and getting gene names from ensmusg IDs. 
# !!! THIS IS FOR NEW scRNA DATA !!!
RenameGenesSeurat <- function(obj, newnames) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  obj[['RNA_name']] <- obj[['RNA']]
  RNA <- obj@assays$RNA_name
  
  #tmp.conv <- tibble(id=RNA@counts@Dimnames[[1]], symbol=newnames)
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
    #if (length(RNA@scale.data)) dimnames(RNA@scale.data)[[1]]    <- tmp.conv$symbol[match(dimnames(RNA@scale.data)[[1]],tmp.conv$id)]
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA_name <- RNA
  return(obj)
}
