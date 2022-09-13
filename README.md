# NeuronalFeatureSpace

This is subproject of defining joint feature space for mouse CNS E12-E14, more specifically from E12R1, E14DI, E14MB and E14R1 non-multimodal scATAC processed samples.

Open topics include how to validate clades biologically, or is this necessary as such. However, some type of optimization for cut height needs to be done.

Main steps to reproduce analysis goes as described below (Conda env: r411_291021 for step 1, and then r421 for steps after that):

1) Run JointFeatureSpace_220622.Rmd to create joint featurespace.
2) Run R14_scRNA_300822.Rmd to merge and process scRNA samples together.
3) Run E14_clade_peaks_joint_space_290822.Rmd to read scATAC data from fragements files based on joint featurespace, replicates are first merged and then brain regions. Data is QC filtered and binarized.
4) Run e14_counts_to_downstream_290822.Rmd to perfom LSI, RNA integration, clustering resolution optimization and cluster marker detection for each modality.
