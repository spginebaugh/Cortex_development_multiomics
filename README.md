# Cortex development multiomics
Work in progress.

In this repository, I reanalyze human corticogenesis data from [Trevino et al., 2021](https://doi.org/10.1016/j.cell.2021.07.039). This dataset contains single-cell RNAseq, single-cell ATAC seq, and single-cell multiomics data (simultaneous RNA+ATACseq), making it an ideal dataset for testing multiomic methods.

The original manuscript utilized the Seurat implementation of Canonical Correlation Analysis to merge the different sequencing modalities. However, since the initial submission of this manuscript, there have been numerous advances in multi-omics analysis with programs/algorithms like [MultiVI](https://doi.org/10.1038/s41592-023-01909-9), [Muon](https://doi.org/10.1186/s13059-021-02577-8), and [MOFA+](https://doi.org/10.1186/s13059-020-02015-1) (see my analysis of muscle aging using a MOFA based method [here](https://github.com/spginebaugh/muscle_aging_ML)). Additionally, I wanted to test out some of the python-based ATAC-seq methods, like SnapATAC2, as well as some of the python-based multiomic downstream analysis methods, like Scenic+. 

Here, we explore different AI/ML approaches to the analysis of multi-omic data, and identify transcription factors involved in corticogenesis, with applications to understanding Autism Spectrum Disorder and improved differentiation protocols for neural organoids.
