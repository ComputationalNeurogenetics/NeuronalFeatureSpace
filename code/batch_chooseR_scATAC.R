.libPaths(c("/projappl/project_2006402/r422_packages", .libPaths()))
libpath <- .libPaths()[1]
library(qs)
library(Seurat)
library(tidyverse)
source("AdditionalFunctions.R")
run.date<-"280823"
s.data <- qread(paste("../analysis/E14_scATAC_integrated.",run.date,".qs",sep=""), nthreads = 6)
doChooseR(obj=s.data,npcs=2:59,resolutions=seq(1,3,by=1),assay="peaks",reduction="lsi",results_path=createDir("../analysis/results_E14_280823_scATAC_leiden_test/"),cores=3,method="igraph",algorithm=4,iterations=3)
