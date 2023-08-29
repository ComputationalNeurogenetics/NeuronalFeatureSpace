#.libPaths(c("/projappl/project_2006402/r422_packages", .libPaths()))
#libpath <- .libPaths()[1]

setwd("/home/kilpinen/NeuronalFeatureSpace/chooseR/")

source("R/pipeline.R")
library(Seurat)
library(writexl)
library(parallel)
library(qs)
`%>%` <- magrittr::`%>%`

# Load Seurat object
# Make sure this has already been normalised and that PCA has been performed
# The data included here is the Ding, et al., Nature Biotechnology, 2020
# human PBMC Smart-Seq data from our manuscript
#obj <- readRDS("data/ding_smartseq_pbmc_preprocessed.rds")
obj <- qread("../analysis/E14_scATAC_integrated.280823.qs",nthreads = 20)

# Define the number of PCs to use, and which assay and reduction to use.
# We recommend testing a broad range of resolutions
# For more on picking the correct number of PCs, see:
# https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
npcs <- 2:59
resolutions <- seq(1,20,by=1)
assay <- "peaks"
reduction <- "lsi"
results_path <- "../analysis/results_E14_280823_scATAC_leiden/"

# Run pipeline
mclapply(resolutions,function(res){
  message(paste0("Clustering ", res, "..."))
  message("\tFinding ground truth...")

  # "Truths" will be stored at glue::glue("{reduction}.{assay}_res.{res}")
  obj <- find_clusters(
                       obj,
                       reduction = reduction,
                       assay = assay,
                       resolution = res,
                       npcs = npcs,
		       method="igraph",
                       algorithm=4
  )
  clusters <- obj[[glue::glue("{reduction}.{assay}_res.{res}")]]

  # Now perform iterative, sub-sampled clusters
  results <- multiple_cluster(
                              obj,
                              n = 100,
                              size = 0.8,
                              npcs = npcs,
                              res = res,
			      method="igraph",
                              reduction = reduction,
                              assay = assay,
                              algorithm=4
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

  # Finally, calculate grouped metrics
  message(paste0("Grouping ", res, "..."))
  grp <- group_scores(mtchs, unlist(clusters))
  saveRDS(grp, paste0(results_path, "frequency_grouped_", res, ".rds"))
  sil <- group_sil(sil, res)
  saveRDS(sil, paste0(results_path, "silhouette_grouped_", res, ".rds"))
  
  # Store cells per cluster table per resolution
  saveRDS(table(Idents(obj)), paste0(results_path, "cells_per_cluster__", res, ".rds"))
}, mc.cores=3)

# Save original data, with ground truth labels
saveRDS(obj, paste0(results_path, "clustered_data.rds"))

