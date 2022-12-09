# Neuron types in the developing mouse CNS can be divided into several epigenomic and transcriptomic classes

This is repository for code used in the manuscript "Neuron types in the developing mouse CNS can be divided into several epigenomic and transcriptomic classes".

Main steps to reproduce analysis goes as described below by using Conda envs r411_291021 for step 1, and then r421 for steps after that. ViaEnv for VIA analysis.

1) Run JointFeatureSpace_220622.Rmd to create joint featurespace.
2) Run R14_scRNA_071222.Rmd to merge and process scRNA samples together.
3) Run E14_clade_peaks_joint_space_290822.Rmd to read scATAC data from fragements files based on joint featurespace, replicates are first merged and then brain regions. Data is QC filtered and binarized.
4) Run e14_counts_to_downstream_290822.Rmd to perfom LSI, RNA integration, clustering resolution optimization and cluster marker detection for each modality.

Essential steps in e14_counts_to_downstream_071222.Rmd are as follows:
1.	Read in both scATAC and scRNA data objects containing merged E14 data
2.	Run LSI, find top (95%) features
3.	First singular component was dropped due to correlating with read depth
4.	Create gene activity matrix based on estimated RNA levels
5.	Find transfer anchors between scATAC and scRNA samples
6.	Perform label transfer and RNA imputation
7.	Filter scATAC cells to >0.5 scRNA integration prediction score
8.	Perform identification of neurotransmitter type per cell
9.	Calculate UMAP and cluster detection optimization (separate HPC run code for optimization)
10.	Pairwise differential chromatin accessibility calculations (Separate HPC run code)
11.	Adding motif information to the object
12.	Run ChromVar
13.	Find Closest features for ATAC features
14.	Identification of cluster markers (not actively used but included to the object for future research)
a.	Expression markers
b.	scATAC feature markers
c.	TF motif accessibility
d.	TF motif enrichments

Additional analyses are coded in ClusterCircosPlots.Rmd
