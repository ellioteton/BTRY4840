Author: Elliot Eton
Due Date: December 17, 2022
Course: BTRY4840

Document title: README.md
Description: Here, I outline step-by-step how to recapitulate the findings in my final project report.

I performed all analysis on a dataset from a single patient (M12) with a JAK2-mutant-driven myeloproliferative neoplasm. To generate the dataset, the single cell method DOGMA-seq was performed. The ATAC (chromatin accessibility) and RNA (gene expression) libraries were sequenced on an Illumina 10X machine. FASTQs were constructed and loaded into the Cellranger pipeline, which processes cell-barcoded reads and performs STAR transcriptome alignment, unique molecular identifier (UMI) correction, and fragment/read counting. The resultant ATAC and RNA matrices were loaded into a Seurat object and subsequently filtered through a range of quality control measures. 

All ATAC fragments and gene expression (GEX) reads were aligned to the human reference genome GRCh38 using Cellranger (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger). 


This dataset includes chromatin accessibility data (ATAC data) and gene expression data (RNA data).

Step 1: 
