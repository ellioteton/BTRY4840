link_peaks <- function(sample_path,sample_name,output_path){
  message("Loading packages")
  library(Seurat)
  library(Signac)
  library(stringr)
  library(dplyr)
  library(BSgenome.Hsapiens.UCSC.hg38)
  
  set.seed(10021)
  
  message("Output path is ",output_path)
  message("Loading sample...")
  sample <- readRDS(sample_path)
  message("Calculating region stats...")
  sample <- RegionStats(sample, BSgenome.Hsapiens.UCSC.hg38, assay = "Peaks", verbose = TRUE)
  message("Linking peaks...")
  sample <- LinkPeaks(object = sample,
                      peak.assay = "Peaks",
                      expression.assay = "RNA"
                      )
  message("Saving object...")
  saveRDS(object = sample,file = paste0(output_path,sample_name,"_linked.rds"))
  
  #### Save links matrix ####
  message("Saving links matrix...")
  links <- GetAssayData(object = sample,assay = 'Peaks', slot = 'links')
  saveRDS(object = links,file = paste0(output_path,sample_name,"link_df.rds"))
}

message(Sys.time())
message("##########################Starting linking##########################")
link_peaks("~/05_M12/06_integrated/M12_final_ATAC_RNA_Peaks.rds",
           "M12",
           "~/05_M12/07_cp/")
message("##########################Ending linking##########################")
message(Sys.time())