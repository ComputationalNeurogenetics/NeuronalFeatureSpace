# NeuronalFeatureSpace

This is subproject of defining joint feature space for mouse CNS E12-E14, more specifically from E12R1, E14Di, E14vMB and E14R1 non-multimodal scATAC processed samples.

Open topics include how to validate clades biologically, or is this necessary as such. However, some type of optimization for cut height needs to be done.


Creation of joint scATAC feature space data goes along these steps.

Conda env: r411_221021

First phase:

Per each sample location (E12R1, E14DI, E14vMB, E14vR1)

1) Read cell-level information and fragments files for each replicate

2) Create feature bin x cell matrix per replicate
  - Binning mm10 to 5kbp bins and counting reads in them per cell
  - Returning sparse feature matrix 

3) Create chromatin assay object per replicate

4) Create Seurat object per replicate

5) Merging replicates

7) Storing resulting Seurat object as Rds file


Second phase:

1) Read Seurat objects with merged data from previous phase for all sample locations

2) Creating new merged Seurat object from all of samples

3) Checking the size, and presence of fragment file links of fully merged object

4) Add gene annotations from Ensembl

5) Calculate nucleosome signal

6) Plotting Nucleosome signal

7) Calculate TSS enrichment

8) Plot biological QC data

9) Filter based on QC values (peak_region_fragments, pct_reads_in_peaks, blacklist_ratio, nucleosome_signal, TSS.enrichment)

10) Filtering based on number of cells per feature or vice versa

11) Binarize data, calculate TFIDF and SVD

12) Correlation between depth and reduced dimension components

13) Take approximate right singular vectors, matrix multiplication with singular values and transpose, dropping first component as it correlates with read depth

14) Perform hierarchical clustering of cells (cosine, ward.D2) and cut at specific level to get clades with about balanced distribution of cells

15) Run MACS for clade level peak detection (effective.genome.size = 1.87e9, group.by = "clade", combine.peaks = TRUE, other defaults for Signac::CallPeaks)

16) Store resulting Granges object as Rds file, this is end result of the entire analysis. Granges object containind coordinates of features to be actually used in the processing of these samples in their separate analyses.
