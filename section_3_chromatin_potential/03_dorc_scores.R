sample_path<-"~/05_M12/07_cp/M12_linked.rds"
sample_name<-"M12"
DORCs_path<-"~/05_M12/07_cp/02_DORCs.csv"
links_table_path<-"~/05_M12/07_cp/02_DORCs_links_df.rds"
output_path<-"~/05_M12/07_cp/"


calculate_DORC_scores <- function(sample_path,sample_name,DORCs_path,links_table_path,output_path){
  message("Loading packages")
  library(Seurat)
  library(Signac)
  library(stringr)
  library(plyr)
  library(dplyr)
  library(ggplot2)
  
  set.seed(10021)
  
  message("Loading sample ",sample_name,"...")
  sample <- readRDS(sample_path)
  
  message("Loading list of DORCs...")
  DORCs<-read.csv(DORCs_path)
  
  message("Loading DORCs links table...")
  DORCs_links <- readRDS(links_table_path)
  DORC_peaks<-unique(DORCs_links$links.peak)
  
  message("Output path is ",output_path)
  
  library(Matrix)
  message("First: normalize peak counts by the total number of unique fragments in peaks per cell")
  peaks <- GetAssayData(object = sample,assay = 'Peaks', slot = 'data')
  peaks@x <- peaks@x / rep.int(colSums(peaks), diff(peaks@p))
  
  message("Second: subset complete peaks matrix to contain only peaks significantly correlated with given DORCs")
  subset_peaks <- peaks[rownames(peaks) %in% DORC_peaks,]
  message("Convert peak names to gene names")
  rownames(subset_peaks) <- mapvalues(rownames(subset_peaks),DORCs_links$links.peak,DORCs_links$links.gene)
  message("Calculate raw DORC score per gene per cell")
  library(Matrix.utils)
  DORC_scores <- aggregate.Matrix(subset_peaks,groupings=rownames(subset_peaks),fun="sum")
  
  message("Add DORC scores matrix to Seurat object")
  DORC_scores_assay <- CreateAssayObject(counts = DORC_scores)
  sample[["DORC"]]<-DORC_scores_assay
  saveRDS(sample,paste0(output_path,"03_M12_w_DORC_scores.rds"))
  
  message("Calculate an average DORC score by cell type")
  #df_DORC_scores_by_cell <- data.frame(DORC_scores)
  #colnames(df_DORC_scores_by_cell)<-colnames(DORC_scores)
  #df_DORC_scores_by_cell <- t(df_DORC_scores_by_cell)

  #sample<-AddMetaData(object = sample,metadata = df_DORC_scores_by_cell)
  
  clusters <- sort(as.numeric(as.character(unique(sample$wsnn_res.0.8))))
  cluster<-0

  df <- data.frame(matrix(nrow = nrow(DORC_scores)))
  rownames(df) <- rownames(DORC_scores)
  DORC_scores <- GetAssayData(object = sample,slot = "data",assay = "DORC")
  
  for (cluster in clusters){
    #Need to get barcodes corresponding to cluster from Seurat object
    cluster_cells <- Cells(sample)[sample$wsnn_res.0.8 == cluster]
    #Subset DORC_scores matrix by these barcodes
    cluster_DORC_scores <- DORC_scores[,colnames(DORC_scores) %in% cluster_cells]
    #Average across each row
    cluster_average <- rowMeans(cluster_DORC_scores)
    #Save averages to new dataframe, where header = cluster, rownames = genes
    df[[as.character(cluster)]]<-cluster_average
  }

  clusters_DORC_scores <- df[,2:ncol(df)]

  # custom function to implement min max scaling
  minMax <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }

  #normalise data using custom function
  clusters_DORC_scores_minmax_norm <- data.frame(lapply(clusters_DORC_scores, minMax))
  rownames(clusters_DORC_scores_minmax_norm)<-rownames(clusters_DORC_scores)
  clusters_DORC_scores_minmax_norm$gene <- rownames(clusters_DORC_scores_minmax_norm)

  library(reshape2)
  melt_DORC_scores <- melt(clusters_DORC_scores_minmax_norm)
  colnames(melt_DORC_scores)<-c("Gene","Cluster","Score")

  DORC_scores_by_cell_type<-ggplot(melt_DORC_scores,aes(x=Gene,y=Cluster,fill=Score))+
    geom_tile() +
    guides (fill = guide_colourbar(title="DORC score")) +
    xlab('DORC') +
    ylab('ATAC cluster') +
    ggtitle("DORC scores by cluster")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_fill_gradient(low = "white",
                        high = "darkgreen",
                        guide = "colorbar")

  ggsave(paste0(output_path,"03_DORC_scores.jpeg"),plot=DORC_scores_by_cell_type,device="jpeg",width=20,height=10)

  #subset by most representative DORCs
  top_DORCs <- na.exclude(c(colnames(t(clusters_DORC_scores_minmax_norm))[max.col(t(clusters_DORC_scores_minmax_norm),ties.method="first")]))
  top_DORCs <- unique(top_DORCs)

  DORC_scores_by_cell_type_top_reps<-ggplot(melt_DORC_scores[melt_DORC_scores$Gene %in% top_DORCs,],aes(x=Gene,y=Cluster,fill=Score))+
    geom_tile() +
    guides (fill = guide_colourbar(title="DORC score")) +
    xlab('DORC') +
    ylab('ATAC cluster') +
    ggtitle("DORC scores by cluster")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_fill_gradient(low = "white",
                        high = "darkgreen",
                        guide = "colorbar")
  ggsave(paste0(output_path,"03_top_DORC_scores.jpeg"),plot=DORC_scores_by_cell_type_top_reps,device="jpeg",width=10,height=10)
}

message(Sys.time())
message("##########################Starting DORC score calculation##########################")
calculate_DORC_scores(
  "~/05_M12/07_cp/M12_linked.rds",
  "M12",
  "~/05_M12/07_cp/02_DORCs.csv",
  "~/05_M12/07_cp/02_DORCs_links_table.rds",
  "~/05_M12/07_cp/")
message("##########################Ending DORC score calculation##########################")
message(Sys.time())