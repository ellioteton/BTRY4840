id_dorcs <- function(sample_path,output_path){
  message("Loading packages")
  library(ggrepel)
  library(Seurat)
  library(Signac)
  library(stringr)
  library(dplyr)
  library(ggplot2)
  
  set.seed(10021)
  
  message("Read in links dataframe")
  links <- readRDS(sample_path)
  
  message("For peaks associated with multiple genes, keep only the peak-gene association with the smallest p-value.")
  links <- links[order(links$pvalue, decreasing=FALSE),]
  linksKeepSmallest <- links[!duplicated(links$peak),]
  links <- linksKeepSmallest
  
  message("Saving dataframe of kept peaks with each linked gene")
  links_df <- data.frame(links$score,links$gene,links$peak,links$zscore,links$pvalue)
  saveRDS(object = links_df,file = paste0(output_path,"02_DORCs_links_df.rds"))
  
  message("Rank order genes by number correlated peaks")
  links_table <- data.frame(as.matrix(table(links_df$links.gene)))
  #Sort by descending
  colnames(links_table) <- "num"
  links_table$genes <- rownames(links_table)
  links_table <- links_table[order(-links_table$num),]
  #rank order
  links_table$rank <- c(seq(1,length(links_table$genes)))
  saveRDS(object = links_df,file = paste0(output_path,"02_DORCs_links_table.rds"))
  
  message("Identify and save list of DORCs")
  DORCgenes<-unique(links_table$genes)
  write.csv(x = DORCgenes,file = paste0(output_path,"02_DORCs.csv"))
  
  message("Visualize DORCs")
  DORCs_plot<-ggplot(data=links_table, aes(x=rank, y=num, label=genes)) +
    geom_line()+
    geom_point()+
    xlab("Ranked genes")+
    ylab("Number of correlated peaks")+
    scale_x_reverse()+
    geom_hline(yintercept=10, linetype="dashed", color = "gray")+
    geom_point(data=links_table[links_table$num>=10, ], aes(x=rank, y=num), colour="red")+
    theme_classic()+
    geom_text_repel(aes(label=ifelse(num>20,as.character(genes),'')),
                    max.overlaps=10,
                    color="blue",
                    min.segment.length = 0,
                    hjust=0)+
    ggtitle("DORCs identified in study of hematopoiesis in M12 patient sample")
  ggsave(filename = paste0(output_path,"02_DORCs_plot.jpeg"),plot = DORCs_plot,device = "jpeg")
}

message(Sys.time())
message("##########################Starting to ID DORCs##########################")
id_dorcs("~/05_M12/07_cp/M12link_df.rds",
         "~/05_M12/07_cp/")
message("##########################Ending ID of DORCs##########################")
message(Sys.time())
