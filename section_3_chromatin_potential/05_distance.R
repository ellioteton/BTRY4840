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

sample_path<-"~/05_M12/07_cp/04_completed_knn.rds"
sample_name<-"M12"
output_path<-"~/05_M12/07_cp/"

message("Loading sample ",sample_name,"...")
sample <- readRDS(sample_path)

calculate_distance<-function(distance_metric){
  #Construct distance matrix
  d<-dist(sample@reductions[["wnn.umap"]]@cell.embeddings, method = distance_metric)
  #For each chromatin profile, calculate the 10 RNA nearest neighbors
  sa<-dbscan::kNNdist(x = d,k=10,all=T)
  distance_df <- data.frame(sa)
  #Calculate average of the 10 distances
  distance_df$average <- apply(distance_df,1,mean)
  #Min max normalize the averages: the output is the chromatin potential
  min_average <- min(distance_df$average)
  max_average <- max(distance_df$average)
  distance_df$norm_average <- (distance_df$average - min_average) / (max_average-min_average)
  #Add chromatin potential to Seurat object
  chromatin_potential <- distance_df[,12,drop=F]
  colnames(chromatin_potential)<-"chromatin_potential"
  sample<-AddMetaData(sample,chromatin_potential)
  #Visualize chromatin potential
  a <- VlnPlot(sample,features="chromatin_potential",group.by="predicted.l2.plot")
  ggsave(paste0(output_path,distance_metric,"_vln_plot.jpeg"),plot = a,device = "jpeg",width = 10,height=7)
  b <- FeaturePlot(sample,features="chromatin_potential",reduction = "dorc_rna_wnn.umap")
  ggsave(paste0(output_path,distance_metric,"_umap.jpeg"),plot = b,device = "jpeg",width = 10,height=7)
}
  
calculate_distance("euclidean")

message(Sys.time())