
gene_set_prepare2 <- function(path_to_db_file, cell_type){
  
  cell_markers = openxlsx::read.xlsx(path_to_db_file)
  cell_markers = cell_markers[cell_markers$tissueType == cell_type,] 
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  # correct gene symbols from the given DB (up-genes)
  cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    paste0(markers_all, collapse=",")
  })
  
  # correct gene symbols from the given DB (down-genes)
  cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    paste0(markers_all, collapse=",")
  })
  
  cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
  
  list(gs_positive = gs, gs_negative = gs2)
}

convert_feature_identity <- function(object, assay, features, feature.format = "symbol") {
  
  #'
  #' Converts ENS ID -> gene symbol and vice versa
  #' Returns a vector of length(features) of either matches or NAs, in corresponding indices.
  #'
  #' Assumes libraries dplyr and seuratObject.
  #' Moreover, requires seuratObject[["assay"]] to contain df/tbl of ENS-symbol correspondences.
  #'
  # Protective tests
  if (!any(feature.format %in% c("ens", "symbol"))) {
    stop("Feature format should be a sting: either 'symbol' or 'ens'")
  }
  if (length(features) == 0) {
    stop("No features found in argument 'features'")
  }
  if (feature.format == "ens" & !all(grepl("*ENS", features))) {
    message("Warning: Found non-ENS ID for argument feature format 'ens'")
  }
  if (feature.format == "symbol" & any(grepl("*ENS", features))) {
    message("Warning: Found ENS ID for argument feature format 'symbol'")
  }
  # Diverging execution: case if provided features are ENSEBL IDs => conversion to symbols
  if (feature.format == "ens") {
    
    object.features <- object[[assay]][[]] %>%
      rownames_to_column(var = "gene_id") %>%
      as_tibble() %>%
      dplyr::select("gene_id", "feature_symbol") %>%
      dplyr::filter(gene_id %in% features)
    match.index <- match(features, object.features$gene_id, nomatch = NA)
    v.out <- sapply(match.index, function (i) { ifelse(is.na(i), NA, object.features$feature_symbol[i])})
    
    sprintf("Instance: Found matching for %d features out of total %d provided features", sum(!is.na(v.out)), length(features)) %>%
      print()
    
    return (v.out)
  }
  
  # Case: otherwise provided symbols => conversion to ENS IDs
  object.features <- object[[assay]][[]] %>%
    rownames_to_column(var = "gene_id") %>%
    as_tibble() %>%
    dplyr::select("gene_id", "feature_symbol") %>%
    dplyr::filter(feature_symbol %in% features)
  
  match.index <- match(features, object.features$feature_symbol, nomatch = NA)
  v.out <- sapply(match.index, function (i) { ifelse(is.na(i), NA, object.features$gene_id[i])})
  
  sprintf("Instance: Found matching for %d features out of total %d provided features", sum(!is.na(v.out)), length(features)) %>%
    print()
  
  return (v.out)
  }
  
create_dt <- function(x){
    DT::datatable(x,
                  extensions = 'Buttons',
                  options = list(dom = 'Blfrtip',
                                 buttons = c('copy', 'csv', 'excel', 'pdf'),
                                 lengthMenu = list(c(10,25,50,-1),
                                                   c(10,25,50,"All"))))
}
  
# Some of following code has been sourced from Stackoverflow forums as freely provided examples, contributions belong to original authors
  
prune_cutree_to_dendlist <- function(dend, k) {
  clusters <- cutree(dend,k, order_clusters_as_data = FALSE)
  dends <- vector("list", k)
  for(i in 1:k) { 
    leves_to_prune <- labels(dend)[clusters != i]
    dends[[i]] <- prune(dend, leves_to_prune)
  }
  class(dends) <- "dendlist"
  dends
}

get_attribute_counts <- function(seurat.data, clusters, all.clusters=FALSE, attribute=NULL){
  if (all.clusters){
    i <- TRUE
  } else {
    i <- which(seurat.data@meta.data[,"seurat_clusters"] %in% clusters)
  }
  return(table(seurat.data@meta.data[i,attribute]))
}

dendro_data_k <- function(hc, k) {
  
  hcdata    <-  ggdendro::dendro_data(hc, type = "rectangle")
  seg       <-  hcdata$segments
  labclust  <-  cutree(hc, k)[hc$order]
  segclust  <-  rep(0L, nrow(seg))
  heights   <-  sort(hc$height, decreasing = TRUE)
  height    <-  mean(c(heights[k], heights[k - 1L]), na.rm = TRUE)
  
  for (i in 1:k) {
    xi      <-  hcdata$labels$x[labclust == i]
    idx1    <-  seg$x    >= min(xi) & seg$x    <= max(xi)
    idx2    <-  seg$xend >= min(xi) & seg$xend <= max(xi)
    idx3    <-  seg$yend < height
    idx     <-  idx1 & idx2 & idx3
    segclust[idx] <- i
  }
  
  idx                    <-  which(segclust == 0L)
  segclust[idx]          <-  segclust[idx + 1L]
  hcdata$segments$clust  <-  segclust
  hcdata$segments$line   <-  as.integer(segclust < 1L)
  hcdata$labels$clust    <-  labclust
  
  hcdata
}

plot_ggdendro <- function(hcdata,
                          direction   = c("lr", "rl", "tb", "bt"),
                          fan         = FALSE,
                          scale.color = NULL,
                          branch.size = 1,
                          label.size  = 3,
                          nudge.label = 0.01,
                          expand.y    = 0.1) {
  
  direction <- match.arg(direction) # if fan = FALSE
  ybreaks   <- pretty(segment(hcdata)$y, n = 5)
  ymax      <- max(segment(hcdata)$y)
  
  ## branches
  p <- ggplot() +
    geom_segment(data         =  segment(hcdata),
                 aes(x        =  x,
                     y        =  y,
                     xend     =  xend,
                     yend     =  yend,
                     linetype =  factor(line),
                     colour   =  factor(clust)),
                 lineend      =  "round",
                 show.legend  =  FALSE,
                 size         =  branch.size)
  
  ## orientation
  if (fan) {
    p <- p +
      coord_polar(direction = -1) +
      scale_x_continuous(breaks = NULL,
                         limits = c(0, nrow(label(hcdata)))) +
      scale_y_reverse(breaks = ybreaks)
  } else {
    p <- p + scale_x_continuous(breaks = NULL)
    if (direction %in% c("rl", "lr")) {
      p <- p + coord_flip()
    }
    if (direction %in% c("bt", "lr")) {
      p <- p + scale_y_reverse(breaks = ybreaks)
    } else {
      p <- p + scale_y_continuous(breaks = ybreaks)
      nudge.label <- -(nudge.label)
    }
  }
  
  # labels
  labelParams <- set_labels_params(nrow(hcdata$labels), direction, fan)
  hcdata$labels$angle <- labelParams$angle
  
  p <- p +
    geom_text(data        =  label(hcdata),
              aes(x       =  x,
                  y       =  y,
                  label   =  label,
                  colour  =  factor(clust),
                  angle   =  angle),
              vjust       =  labelParams$vjust,
              hjust       =  labelParams$hjust,
              nudge_y     =  ymax * nudge.label,
              size        =  label.size,
              show.legend =  FALSE)
  
  # colors and limits
  if (!is.null(scale.color)) {
    p <- p + scale_color_manual(values = scale.color)
  }
  
  ylim <- -round(ymax * expand.y, 1)
  p    <- p + expand_limits(y = ylim)
  
  p
}

set_labels_params <- function(nbLabels,
                              direction = c("tb", "bt", "lr", "rl"),
                              fan       = FALSE) {
  if (fan) {
    angle       <-  360 / nbLabels * 1:nbLabels + 90
    idx         <-  angle >= 90 & angle <= 270
    angle[idx]  <-  angle[idx] + 180
    hjust       <-  rep(0, nbLabels)
    hjust[idx]  <-  1
  } else {
    angle       <-  rep(0, nbLabels)
    hjust       <-  0
    if (direction %in% c("tb", "bt")) { angle <- angle + 45 }
    if (direction %in% c("tb", "rl")) { hjust <- 1 }
  }
  list(angle = angle, hjust = hjust, vjust = 0.5)
}

construct_ident_vector <- function(tree.branches, old.idents){
  new.ident <- rep(NA, length(old.idents))
  branches.clusters <- as_tibble(tree.branches) %>% dplyr::select(prunned_tree, clusters) %>% distinct()
  sapply(1:nrow(branches.clusters), function(r){
    ii <- which(old.idents %in% branches.clusters[r,"clusters"][[1]][[1]])
    new.ident[ii] <<-r
    })
  return(new.ident)
}

prepare_comp_circos_data <- function(da_features=da_features, id.1, id.2, avg_log2FC.limit=0.75, s.data){
  
  feat.1 <- DA_feats(da.features = da_features, id.1=id.1, id.2=id.2, avg_log2FC.limit = avg_log2FC.limit)
  feat.2 <- DA_feats(da.features = da_features, id.1=id.2, id.2=id.1, avg_log2FC.limit = avg_log2FC.limit)
  
  dat.1 <- map(feat.1, function(x) {
    tmp <- unlist(str_split(x, pattern = "-"))
    tibble(chr = tmp[1], start=tmp[2], end=tmp[3],side=id.1) 
  }) %>% bind_rows() %>% mutate_at(c("start", "end"), .funs=as.numeric)
  
  dat.2 <- map(feat.2, function(x) {
    tmp <- unlist(str_split(x, pattern = "-"))
    tibble(chr = tmp[1], start=tmp[2], end=tmp[3],side=id.2) 
  }) %>% bind_rows() %>% mutate_at(c("start", "end"), .funs=as.numeric)
  
  dat <- bind_rows(dat.1, dat.2)
  
  closest.gene.1.genes <- as_tibble(ClosestFeature(object=s.data, regions=feat.1)) %>% pull(gene_name)
  closest.gene.1.dist <- as_tibble(ClosestFeature(object=s.data, regions=feat.1)) %>% pull(distance)
  
  closest.gene.2.genes <- as_tibble(ClosestFeature(object=s.data, regions=feat.2)) %>% pull(gene_name)
  closest.gene.2.dist <- as_tibble(ClosestFeature(object=s.data, regions=feat.2)) %>% pull(distance)
  
  dat$label <- c(closest.gene.1.genes,closest.gene.2.genes)
  dat$distance <- c(closest.gene.1.dist,closest.gene.2.dist)
  dat$color <- ifelse(dat$side==id.1, "red","blue")
  return(dat)
}

DA_feats <- function(da.features, id.1, id.2, avg_log2FC.limit=0.75, p_val_adj_limit = 0.05){
  da_features[[id.1]][[id.2]] %>% dplyr::filter((avg_log2FC >= avg_log2FC.limit) & p_val_adj <= p_val_adj_limit) %>% pull(feature)
}



RenameGenesSeurat <- function(obj, newnames) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  obj[['RNA_name']] <- obj[['RNA']]
  RNA <- obj@assays$RNA_name
  if (length(RNA@scale.data) > 0){
    tmp.conv <- tibble(id=RNA@counts@Dimnames[[1]], symbol=newnames)
  }
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts) >0 & class(RNA@data)[1]=="dgCMatrix") {RNA@counts@Dimnames[[1]]            <- newnames}
    if (length(RNA@data) >0 ){ RNA@data@Dimnames[[1]]                <- newnames}
    if (length(RNA@scale.data) > 0 & !is.matrix(RNA@scale.data)){RNA@scale.data@Dimnames[[1]]    <- newnames}
    if (length(RNA@scale.data) > 0 & is.matrix(RNA@scale.data)){rownames(RNA@scale.data)    <- tmp.conv$symbol[match(rownames(RNA@scale.data),tmp.conv$id)]}
    #if (length(RNA@scale.data)) dimnames(RNA@scale.data)[[1]]    <- tmp.conv$symbol[match(dimnames(RNA@scale.data)[[1]],tmp.conv$id)]
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA_name <- RNA
  return(obj)
}

doRNAintegration <- function(scATAC.data, scRNA.data, results.path, cores=1, run.date=run.date, max.lsi.dim){
  # RNA activity estimation
  genomic.metadata <- mcols(Annotation(s.data[['peaks']]))
  # Generate conversion table from gene_name to gene_id
  gene_name2gene_id <- as_tibble(genomic.metadata[,c("gene_name","gene_id")])
  
  # Calculate gene activity estimate from scATAC reads based on the scATAC, by using gene_names as function does not support any other id
  gene.activities <- GeneActivity(s.data, assay="peaks")
  
  # Store gene_names
  gene.names <- rownames(gene.activities)
  
  # Switch sparse matrix to use ensmusg id
  ensmusg.ids <- gene_name2gene_id[match(gene.names,pull(gene_name2gene_id,"gene_name")),] %>% pull("gene_id")
  gene_names <- gene_name2gene_id[match(gene.names,pull(gene_name2gene_id,"gene_name")),] %>% pull("gene_name")
  
  # Dropping NAs
  non.na.i <- !is.na(ensmusg.ids)
  gene.activities.gene_id <- gene.activities[non.na.i,]
  rownames(gene.activities.gene_id) <- ensmusg.ids[non.na.i]
  
  # Add the gene activity matrix to the Seurat object as a new assay
  s.data[['Activity']] <- CreateAssayObject(counts = gene.activities.gene_id)
  s.data <- NormalizeData(
    object = s.data,
    assay = 'Activity',
    normalization.method = 'LogNormalize',
    scale.factor = median(s.data$nCount_Activity)
  )
  
  # Add gene_name to gene_id mapping into the s.data[['Activity']] assays metadata
  s.data[['Activity']]<- AddMetaData(s.data[['Activity']], col.name = "feature_symbol", metadata = gene_names[non.na.i])
  
  #Perform label transfer
  
  s.data_RNA@meta.data$CellType<-s.data_RNA@meta.data$seurat_cluster
  
  # Finding transfer anchors
  transfer.anchors <- FindTransferAnchors(
    reference = s.data_RNA,
    query = s.data,
    reduction = 'cca',
    reference.assay="RNA",
    query.assay = "Activity",
    features = VariableFeatures(object=s.data_RNA)
  )
  
  predicted.labels <- TransferData(   
    anchorset = transfer.anchors,
    refdata = s.data_RNA@meta.data$CellType,
    weight.reduction = s.data[['lsi']],
    dims = 2:max.lsi.dim
  )
  
  predicted.NT.labels <- TransferData(
    anchorset = transfer.anchors,
    refdata = s.data_RNA@meta.data$NT.type,
    weight.reduction = s.data[['lsi']],
    dims = 2:max.lsi.dim
  )
  
  s.data <- AddMetaData(object = s.data, metadata = apply(predicted.NT.labels,1,max), col.name ="RNA.predicted.NT")
  s.data <- AddMetaData(s.data, metadata = predicted.labels)
  
  #Perform scRNA data imputation
  DefaultAssay(s.data_RNA) <- "RNA"
  refdata <- GetAssayData(s.data_RNA, assay = "RNA", slot = "data")
  
  s.data_RNA@meta.data$tech<-"scRNA"
  s.data@meta.data$tech<-"scATAC"
  
  imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = s.data[["lsi"]], dims = 2:max.lsi.dim)
  
  s.data[["RNA"]] <- imputation
  coembed <- merge(x = s.data_RNA, y = s.data)
  
  # Copy feature metadata from s.data_rna to s.data
  s.data_rna.feature.metadata <- s.data_RNA[["RNA"]][[]]
  s.data[["RNA"]] <- AddMetaData(s.data[["RNA"]], metadata = s.data_rna.feature.metadata[rownames(s.data[["RNA"]]),"feature_symbol"], col.name = "feature_symbol")
  
  # Find variable features
  coembed <- FindVariableFeatures(coembed)
  
  # Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both datasets
  coembed <- ScaleData(coembed, do.scale = FALSE)
  coembed <- RunPCA(coembed, verbose = FALSE)
  coembed <- RunUMAP(coembed, dims = 2:30)
  coembed@meta.data$CellType <- ifelse(!is.na(coembed@meta.data$CellType), coembed@meta.data$CellType, coembed@meta.data$predicted.id)
  
  qsave(s.data, file=paste(results.path,"/E14_s.data.integrated.",run.date,".qs",sep=""), nthreads = cores)
  qsave(coembed, file=paste(results.path,"/E14_coembed.",run.date,".qs",sep=""), nthreads = cores)
  return(s.data)
}

