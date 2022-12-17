#Following multivelo tutorial with adaptations

#####import packages
import os
import multivelo as mv
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt

#####set settings
scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_rows', 200)
np.set_printoptions(suppress=True)

#####Reading in unspliced and spliced counts
adata_rna = scv.read('velocyto/M12_dogma_velocyto_output.loom', cache=True)
adata_rna.obs_names = [x.split(':')[1][:-1] + '-1' for x in adata_rna.obs_names]
adata_rna.var_names_make_unique()

# Top 1000 variable genes are used for downstream analyses.
scv.pp.filter_and_normalize(adata_rna, min_shared_counts=10, n_top_genes=1000)

# Load cell annotations
cell_annot = pd.read_csv('02_ee_azimuth_pred.tsv', sep='\t', index_col=0)

adata_rna = adata_rna[cell_annot.index,:]
adata_rna.obs['celltype'] = cell_annot['predicted.celltype.l2']

######Preprocessing the ATAC counts
adata_atac = sc.read_10x_mtx('outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=True, gex_only=False)
adata_atac = adata_atac[:,adata_atac.var['feature_types'] == "Peaks"]

# We aggregate peaks around each gene as well as those that have high correlations with promoter peak or gene expression.
# Peak annotation contains the metadata for all peaks.
# Feature linkage contains pairs of correlated genomic features.
adata_atac = mv.aggregate_peaks_10x(adata_atac, 
                                    'outs/atac_peak_annotation.tsv', 
                                    'outs/analysis/feature_linkage/feature_linkage.bedpe', 
                                    verbose=True)

# Let's examine the total count distribution and remove outliers.
plt.hist(adata_atac.X.sum(1), bins=100, range=(0, 100000));

# We normalize aggregated peaks with TF-IDF.
mv.tfidf_norm(adata_atac)

#####Finding shared barcodes and features between RNA and ATAC

shared_cells = pd.Index(np.intersect1d(adata_rna.obs_names, adata_atac.obs_names))
shared_genes = pd.Index(np.intersect1d(adata_rna.var_names, adata_atac.var_names))
print("Length shared_cells")
print(len(shared_cells))
print("Length shared_genes")
print(len(shared_genes))

# We reload in the raw data and continue with a subset of cells.
adata_rna = scv.read("velocyto/M12_dogma_velocyto_output.loom", cache=True)
adata_rna.obs_names = [x.split(':')[1][:-1] + '-1' for x in adata_rna.obs_names]
adata_rna.var_names_make_unique()

adata_rna = adata_rna[shared_cells, shared_genes]
adata_atac = adata_atac[shared_cells, shared_genes]

scv.pp.normalize_per_cell(adata_rna)
scv.pp.log1p(adata_rna)
scv.pp.moments(adata_rna, n_pcs=30, n_neighbors=50)

adata_rna.obs['celltype'] = cell_annot.loc[adata_rna.obs_names, 'predicted.celltype.l2']
adata_rna.obs['celltype'] = adata_rna.obs['celltype'].astype('category')

scv.tl.umap(adata_rna)

scv.pl.umap(adata_rna, 
            color='celltype',
           legend_loc = "right margin")

#####Smoothing gene aggregagted peaks by neighbors

# Write out filtered cells and prepare to run Seurat WNN --> R script can be found on Github.
adata_rna.obs_names.to_frame().to_csv('filtered_cells.txt', header=False, index=False)

# Read in Seurat WNN neighbors.
nn_idx = np.loadtxt("seurat_wnn/nn_idx.txt", delimiter=',')
nn_dist = np.loadtxt("seurat_wnn/nn_dist.txt", delimiter=',')
nn_cells = pd.Index(pd.read_csv("seurat_wnn/nn_cells.txt", header=None)[0])

for element in adata_atac.obs_names:
    if element not in nn_cells:
        print(element)

shared_cells = pd.Index(np.intersect1d(adata_rna.obs_names, nn_cells))
len(shared_cells)

adata_rna = adata_rna[shared_cells, shared_genes]
adata_atac = adata_atac[shared_cells, shared_genes]


# Make sure cell names match.
np.all(nn_cells == adata_atac.obs_names)

mv.knn_smooth_chrom(adata_atac, nn_idx, nn_dist)

adata_atac

#####Running multi-omic dynamical model

# This will take a while. Parallelization is high recommended.
adata_result = mv.recover_dynamics_chrom(adata_rna, 
                                         adata_atac, 
                                         max_iter=5, 
                                         init_mode="invert", 
                                         verbose=False, 
                                         parallel=True, 
                                         save_plot=False,
                                         rna_only=False,
                                         fit=True,
                                         n_anchors=500, 
                                         extra_color_key='celltype'
                                        )

# Save the result for use later on
adata_result.write("multivelo_result.h5ad")
