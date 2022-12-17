sample_path<-"~/05_M12/07_cp/03_M12_w_DORC_scores.rds"
sample_name<-"M12"
output_path<-"~/05_M12/07_cp/"

message(Sys.time())

message("Loading packages")
library(Seurat)
library(Signac)
library(stringr)
library(plyr)
library(dplyr)
library(dbscan)
library(ggplot2)

set.seed(10021)

message("Loading sample ",sample_name,"...")
sample <- readRDS(sample_path)
head(sample[["DORC"]])
sample <- NormalizeData(sample, assay="DORC",normalization.method = "RC", scale.factor = 1e6)
sample <- FindVariableFeatures(sample, assay = "DORC",selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sample[["DORC"]])
sample <- ScaleData(sample,assay="DORC",features = all.genes)

message("Running PCA on DORC and RNA data...")
sample<-RunPCA(object = sample,assay = "DORC",reduction.name = "dorc_pca",reduction.key = "dorcpca_")
a <- ElbowPlot(sample,ndims=50,reduction="dorc_pca")+theme_bw()+ggtitle("Elbow plot for PCA on domains of regulatory chromatin")
ggsave(paste(output_path,"dorc_pca_elbow_plot.jpeg"),plot=a,device="jpeg")
sample<-RunPCA(object = sample,assay = "RNA",reduction.name = "rna_pca",reduction.key = "rnapca_")
b<-ElbowPlot(sample,ndims=50,reduction="rna_pca")+theme_bw()+ggtitle("Elbow plot for PCA on RNA")
ggsave(paste(output_path,"rna_pca_elbow_plot.jpeg"),plot=b,device="jpeg")

message("Define function to construct knn graph using several distance metrics, find clusters accordingly, and plot first two PCs with appropriate labels")
construct_nn_cluster_graph <- function(sample,assay_name,reduction_name,distance_metric,graph_name,output_path,output){
  DefaultAssay(sample)<-assay_name
  if (assay_name == "DORC"){
    sample<-FindNeighbors(object = sample,reduction = reduction_name,dims = 1:15,k.param = 30,annoy.metric = distance_metric,graph.name = graph_name)
  }
  else{
    sample<-FindNeighbors(object = sample,reduction = reduction_name,dims = 1:30,k.param = 30,annoy.metric = distance_metric,graph.name = graph_name)
  }
  sample<-FindClusters(object = sample,graph.name = graph_name,resolution = 1)
  c <- DimPlot(object = sample,dims = c(1,2),group.by="predicted.l2.plot",reduction = reduction_name,label = T,repel = T,label.box = T)+
    ggtitle(paste0("PCA on ",assay_name,": clusters defined through kNN with ",distance_metric," distance"))+NoLegend()
  ggsave(paste0(output_path,output),plot=c,device="jpeg",height=8,width=10)
  return(sample)
}

message("Find neighbors in chromatin space...")
message("Find Euc neighbors in chromatin space...")
sample<-construct_nn_cluster_graph(sample,"DORC","dorc_pca","euclidean","dorc_pca_euc","/gpfs/commons/home/eeton/05_M12/07_cp/","dorc_pca_euc_dimplot.jpeg")
message("Find Cos neighbors in chromatin space...")
sample<-construct_nn_cluster_graph(sample,"DORC","dorc_pca","cosine","dorc_pca_cos","/gpfs/commons/home/eeton/05_M12/07_cp/","dorc_pca_cos_dimplot.jpeg")
message("Find Manhattan neighbors in chromatin space...")
sample<-construct_nn_cluster_graph(sample,"DORC","dorc_pca","manhattan","dorc_pca_man","/gpfs/commons/home/eeton/05_M12/07_cp/","dorc_pca_man_dimplot.jpeg")
message("Find hamming neighbors in chromatin space...")
sample<-construct_nn_cluster_graph(sample,"DORC","dorc_pca","hamming","dorc_pca_ham","/gpfs/commons/home/eeton/05_M12/07_cp/","dorc_pca_ham_dimplot.jpeg")

message("Find neighbors in RNA space...")
sample<-construct_nn_cluster_graph(sample,"RNA","rna_pca","euclidean","rna_pca_euc","/gpfs/commons/home/eeton/05_M12/07_cp/","rna_pca_euc_dimplot.jpeg")
sample<-construct_nn_cluster_graph(sample,"RNA","rna_pca","cosine","rna_pca_cos","/gpfs/commons/home/eeton/05_M12/07_cp/","rna_pca_cos_dimplot.jpeg")
sample<-construct_nn_cluster_graph(sample,"RNA","rna_pca","manhattan","rna_pca_man","/gpfs/commons/home/eeton/05_M12/07_cp/","rna_pca_man_dimplot.jpeg")
sample<-construct_nn_cluster_graph(sample,"RNA","rna_pca","hamming","rna_pca_ham","/gpfs/commons/home/eeton/05_M12/07_cp/","rna_pca_ham_dimplot.jpeg")

message("Find neighbors across chromatin and RNA space...")
sample<-FindMultiModalNeighbors(object = sample,
                                reduction.list = list("dorc_pca","rna_pca"),
                                dims.list=list(1:15, 1:30),
                                k.nn = 10,
                                knn.graph.name = "dorc_rna_wknn",
                                snn.graph.name = "dorc_rna_wsnn",
                                weighted.nn.name = "dorc_rna_weighted.nn")
sample <- RunUMAP(sample, nn.name = "dorc_rna_weighted.nn", reduction.name = "dorc_rna_wnn.umap", reduction.key = "dorc_rna_wnnUMAP_")
sample <- FindClusters(sample, graph.name = "dorc_rna_wsnn", algorithm = 3, resolution = 1, verbose = FALSE)
p1 <- DimPlot(sample, reduction = 'dorc_rna_wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p2 <- DimPlot(sample, reduction = 'dorc_rna_wnn.umap', group.by = 'predicted.celltype.l2', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
ggsave(paste0(output_path,"weighted_multimodal_umap.jpeg"),plot=p1+p2,device="jpeg",width=10,height=5)

saveRDS(sample,paste0(output_path,"04_completed_knn.rds"))
message(Sys.time())