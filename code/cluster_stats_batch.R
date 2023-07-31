.libPaths(c("/projappl/project_2001539/r411_packages", .libPaths()))
libpath <- .libPaths()[1]

library(dplyr)
library(dendextend)
library(fpc)
library(qs)
library(Seurat)
library(Signac)
library(parallel)

run.date<-"260723"

neuronal.merged.filt <- qread("../scATAC_data/neuronal.merged.clade.phase.filt.260723.qs",nthreads=6)
hclust_cells <- qread("../scATAC_data/hclust.260723.qs",nthreads=6)
cell.dist <- qread("../scATAC_data/cell.dist.260723.qs",nthreads=6)

heights <- seq(4,10,by=.5)

cluster.stats.list <- mclapply(heights,function(h){
	cells_tree_cut = cutree(hclust_cells, h=h)
	lsi_cells = dplyr::tibble(barcode = Cells(neuronal.merged.filt), cells_tree_cut = cells_tree_cut)
	tmp.stats <- cluster.stats(d=cell.dist, clustering=lsi_cells$cells_tree_cut)
	return(tmp.stats)
},mc.cores=6)

qsave(cluster.stats.list,file=paste("../analysis/cluster.stats.list.",run.date,".qs",sep=""), nthreads=6)
