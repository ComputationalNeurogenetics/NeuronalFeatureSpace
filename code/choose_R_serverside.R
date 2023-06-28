.libPaths(c("/projappl/project_2001539/r411_packages", .libPaths()))
libpath <- .libPaths()[1]

require(parallel)
require(Seurat)
require(Signac)
require(qs)

source("./R/pipeline.R")

#s.data <- qread("../NeuronalFeatureSpace/scATAC_data/s.data.merged.for.chooseR_q0.qs", nthreads=4)

s.data <- qread("../NeuronalFeatureSpace/scRNA_data/E14_scRNA_for_chooseR_280623.qs", nthreads=6)

#std.proportion.lsi<- sapply(2:length(s.data@reductions$lsi@stdev),function(x){s.data@reductions$lsi@stdev[x]/s.data@reductions$lsi@stdev[2]})

std.proportion.pca <- sapply(2:length(s.data@reductions$pca@stdev),function(x){s.data@reductions$pca@stdev[x]/s.data@reductions$pca@stdev[2]})

max.dim <- which.min(abs(std.proportion.pca-0.15))
npcs <- 2:max.dim

resolutions <- seq(1,20,by=1)
#resolutions <- c(0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
#resolutions.cont <- c(11,12,13,14,15,16,17,18,19,20)
assay <- "RNA"
reduction <- "pca"
results_path <- "results_leiden_RNA/"
algorithm <- 4
random.seed <- 2020

mclapply(resolutions, function(res) {
  message(paste0("Clustering ", res, "..."))
  message("\tFinding ground truth...")

  # "Truths" will be stored at glue::glue("{reduction}.{assay}_res.{res}")
  s.data <- find_clusters(
                       s.data,
                       reduction = reduction,
                       assay = assay,
			algorithm = algorithm,
                       resolution = res,
			random.seed = random.seed,
                       npcs=npcs,
			verbose=FALSE
  )
  clusters <- s.data[[glue::glue("S.{reduction}.{assay}_res.{res}")]]

  # Now perform iterative, sub-sampled clusters
  results <- multiple_cluster(
                              s.data,
                              n = 100,
                              size = 0.8,
                              npcs = npcs,
                              res = res,
				random.seed = random.seed,
                              reduction = reduction,
                              assay = assay,
				algorithm = algorithm
  )

  # Now calculate the co-clustering frequencies
  message(paste0("Tallying ", res, "..."))
  # This is the more time efficient vectorisation
  # However, it exhausts vector memory for (nearly) all datasets
  # matches <- purrr::map(columns, find_matches, df = results)
  # matches <- purrr::reduce(matches, `+`)
  columns <- colnames(dplyr::select(results, -cell))
  mtchs <- matrix(0, nrow = dim(results)[1], ncol = dim(results)[1])
  i <- 1 # Counter
  for (col in columns) {
    message(paste0("\tRound ", i, "..."))
    mtchs <- Reduce("+", list(
      mtchs,
      find_matches(col, df = results)
    ))
    i <- i + 1
  }

  message(paste0("Scoring ", res, "..."))
  mtchs <- dplyr::mutate_all(
    dplyr::as_tibble(mtchs),
    function(x) dplyr::if_else(Re(x) > 0, percent_match(x), 0)
  )

  # Now calculate silhouette scores
  message(paste0("Silhouette ", res, "..."))
  sil <- cluster::silhouette(
    x = as.numeric(as.character(unlist(clusters))),
    dmatrix = (1 - as.matrix(mtchs))
  )
  saveRDS(sil, paste0(results_path, "silhouette_", res, ".rds"))

# Store cells per cluster table per resolution
saveRDS(table(Idents(s.data)), paste0(results_path, "cells_per_cluster__", res, ".rds"))
  # Finally, calculate grouped metrics
  message(paste0("Grouping ", res, "..."))
  grp <- group_scores(mtchs, unlist(clusters))
  saveRDS(grp, paste0(results_path, "frequency_grouped_", res, ".rds"))
  sil <- group_sil(sil, res)
  saveRDS(sil, paste0(results_path, "silhouette_grouped_", res, ".rds"))
  }, mc.cores=10)
