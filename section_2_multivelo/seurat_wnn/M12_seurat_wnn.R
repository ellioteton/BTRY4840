# MultiVelo Seurat WNN Demo
# The procedure mostly follows Seurat tutorial: https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html
# Note that we do not claim these preprocessing steps to be the best, as there aren't any. Feel free to make any changes you deem necessary.
# Please use libblas 3.9.1 and liblapack 3.9.1 for reproducing the 10X mouse M12 demo, or use supplied WNN files on GitHub.

library(Seurat)
library(Signac)
# 
# # read in expression and accessbility data
# M12.data <- Read10X(data.dir = "../outs/filtered_feature_bc_matrix/")
# 
# # subset for the same cells in the jointly filtered anndata object
# barcodes <- read.delim("../filtered_cells.txt", header = F, stringsAsFactors = F)$V1
# 
# # preprocess RNA
# M12 <- CreateSeuratObject(counts = M12.data$`Gene Expression`[,barcodes])
# M12 <- NormalizeData(M12)
# M12 <- FindVariableFeatures(M12)
# M12 <- ScaleData(M12, do.scale = F) # not scaled for consistency with scVelo (optionally, use SCTransform)
# M12 <- RunPCA(M12, verbose = FALSE)
# M12 <- RunUMAP(M12, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_') # optional
# 
# # preprocess ATAC
# M12[["ATAC"]] <- CreateAssayObject(counts = M12.data$`Peaks`[,barcodes], min.cells = 1)
# DefaultAssay(M12) <- "ATAC"
# M12 <- RunTFIDF(M12)
# M12 <- FindTopFeatures(M12, min.cutoff = 'q0')
# M12 <- RunSVD(M12)
# M12 <- RunUMAP(M12, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_") # optional
# 

M12<-readRDS("~/M12/M12_dogma_with_dictionary_mapping.rds")

M12 <- RenameCells(
  M12,
  new.names = gsub("M12_dogma_","",M12$cell)
)

DefaultAssay(M12)<-'ATAC'
M12_w_neighbor_object <- FindNeighbors(M12, 
                     reduction = 'integrated_lsi',
                     dims = 1:50,
                     return.neighbor = TRUE)

# # find weighted nearest neighbors
# M12 <- FindMultiModalNeighbors(M12, 
#                                reduction.list = list("integrated_lsi", "umap"), 
#                                dims.list = list(1:50, 1:2), 
#                                k.nn = 50)
# M12 <- RunUMAP(M12, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") # optional

# extract neighborhood graph
nn_idx <- M12_w_neighbor_object@neighbors$ATAC.nn@nn.idx
nn_dist <- M12_w_neighbor_object@neighbors$ATAC.nn@nn.dist
nn_cells <- M12_w_neighbor_object@neighbors$ATAC.nn@cell.names

# save neighborhood graph
write.table(nn_idx, "/gpfs/commons/home/eeton/M12/06_multivelo/seurat_wnn/nn_idx.txt", sep = ',', row.names = F, col.names = F, quote = F)
write.table(nn_dist, "/gpfs/commons/home/eeton/M12/06_multivelo/seurat_wnn/nn_dist.txt", sep = ',', row.names = F, col.names = F, quote = F)
write.table(nn_cells, "/gpfs/commons/home/eeton/M12/06_multivelo/seurat_wnn/nn_cells.txt", sep = ',', row.names = F, col.names = F, quote = F)

# save sessionInfo for reproducibility
writeLines(capture.output(sessionInfo()), "/gpfs/commons/home/eeton/M12/06_multivelo/seurat_wnn/sessionInfo.txt")

cell_annotations <- M12$predicted.l2
cell_annotations<-data.frame(cell_annotations)
colnames(cell_annotations)<-"predicted.celltype.l2"
cell_annotations$bc <- rownames(cell_annotations)
cell_annotations <- cell_annotations[,c('bc','predicted.celltype.l2')]
write.table(cell_annotations,"/gpfs/commons/home/eeton/M12/06_multivelo/02_ee_dict_pred.tsv",sep='\t',row.names = F,quote=F)
