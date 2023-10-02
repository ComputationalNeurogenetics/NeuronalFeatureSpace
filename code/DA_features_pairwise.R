library(Seurat)
library(tidyverse)
library(qs)
library(Signac)
library(parallel)

setwd("/scratch/project_2001539/SAMI/NeuronalFeatureSpace/scATAC_data")
s.data.e14 <- qread("s.data.res16.qs", nthreads=20)


DefaultAssay(s.data.e14) <- "peaks"
all.cluster.idents <- levels(Idents(s.data.e14))

da_features <- lapply(all.cluster.idents, function(id.1){
  print(paste("Calculating cluster",id.1,sep=" "))
  da_counts <- mclapply(all.cluster.idents, function(id.2){
    if (id.1!=id.2){
      tmp_da_peaks_count <- FindMarkers(
        object = s.data.e14,
        ident.1 = id.1,
        ident.2 = id.2,
        min.pct = 0.05,
        test.use = 'LR',
        latent.vars = 'peak_region_fragments'
      ) %>% rownames_to_column(var="feature") %>% as_tibble()# %>% dplyr::filter((avg_log2FC > 0.75 | avg_log2FC < -0.75) & p_val_adj < 0.05) %>% nrow()
    } else {
      tmp_da_peaks_count <- 0
    }
    return(tmp_da_peaks_count)
  }, mc.cores=40)
  #return(unlist(da_counts))
})

qsave(da_features,"da_features_16.qs", nthreads = 20)
