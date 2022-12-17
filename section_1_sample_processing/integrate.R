#As noted in the README.md file in the GitHub repository, I drew from two Seurat vignettes
#for this integration of ATAC and RNA data. The vignettes are noted in the repository.
#I made many customizations to the vignettes, although the general pseudocode is consistent.

message("Loading packages")
library(Seurat)
library(Signac)
library(stringr)
library(dplyr)
library(ggplot2)
library(EnsDb.Hsapiens.v86)

set.seed(10021)

output_path <- "~/05_M12/06_integrated/"
message("Output path is ",output_path)

####Load M12_dogma RNA data from cellranger
message("Loading data from cellranger")
counts <- Read10X_h5("/gpfs/commons/home/eeton/05_M12/00_cellranger/M12dogma/filtered_feature_bc_matrix.h5")
#extract RNA data
rna_counts <- counts$`Gene Expression`
#extract ATAC data
atac_counts <- counts$Peaks

####Create Seurat object
#Add RNA data
M12 <- CreateSeuratObject(counts=rna_counts)
M12[["percent.mt"]]<-PercentageFeatureSet(M12,pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotations),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotations) <- ucsc.levels

frag.file <- "~/05_M12/00_cellranger/M12dogma/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)

M12[["ATAC"]] <- chrom_assay

per_barcode_metrics <- read.csv("~/05_M12/00_cellranger/M12dogma/per_barcode_metrics.csv")
rownames(per_barcode_metrics)<-per_barcode_metrics$barcode
M12<-AddMetaData(M12,per_barcode_metrics)
DefaultAssay(M12)<-"ATAC"
M12 <- NucleosomeSignal(object = M12)
M12 <- TSSEnrichment(object = M12, fast = FALSE)

####Evaluate quality of the ATAC, RNA, and filter full object based on the quality of both
message("ATAC and RNA QC and filtering")
sample.list <- list(M12sample = M12)
sample<-M12
#### QC and filtering
sample.list.filtered <- lapply(sample.list, function(sample) {

  #Overall QC
  p1 <- VlnPlot(sample, features = c("nCount_ATAC","nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 4)+NoLegend()
  ggsave(paste0(output_path,"qc_pre_filtering.jpeg"),device="jpeg",plot=p1,height=5,width=10)
  p2 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "percent.mt")+theme_bw()+NoLegend()
  ggsave(paste0(output_path,"qc_pre_filtering_nCount_v_mt.jpeg"),device="jpeg",plot=p2,height=5,width=7)
  p3 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+theme_bw()+NoLegend()
  ggsave(paste0(output_path,"qc_pre_filtering_nCount_v_nFeature.jpeg"),device="jpeg",plot=p3,height=5,width=7)
  p3_2 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nCount_ATAC")+theme_bw()+NoLegend()
  ggsave(paste0(output_path,"qc_pre_filtering_nRNA_v_nATAC.jpeg"),device="jpeg",plot=p3_2,height=5,width=7)

  #ATAC QC and filtering
  DefaultAssay(sample)<-"ATAC"

  sample$pct_reads_in_peaks <- sample$atac_peak_region_fragments / sample$atac_fragments * 100
  pct_reads_in_peaks_low <- median(sample$pct_reads_in_peaks) - 2*mad(sample$pct_reads_in_peaks)
  pct_reads_in_peaks_hi <- median(sample$pct_reads_in_peaks) + 2*mad(sample$pct_reads_in_peaks)

  atac_peak_region_fragments_low <- median(sample$atac_peak_region_fragments) - 2*mad(sample$atac_peak_region_fragments)
  atac_peak_region_fragments_hi <- median(sample$atac_peak_region_fragments) + 2*mad(sample$atac_peak_region_fragments)

  TSS_enrichment_low <- median(sample$TSS.enrichment) - 2*mad(sample$TSS.enrichment)
  TSS_enrichment_hi <- median(sample$TSS.enrichment) + 2*mad(sample$TSS.enrichment)

  sample$high.tss <- ifelse(sample$TSS.enrichment > TSS_enrichment_hi,"High","Low")
  p4 <- TSSPlot(sample, group.by='high.tss') + NoLegend() +theme_bw()
  ggsave(paste0(output_path,"qc_pre_filtering_tss_plot.jpeg"),device="jpeg",plot=p4,height=5,width=7)

  sample$blacklist_fraction <- FractionCountsInRegion(object=sample,assay="ATAC",regions=blacklist_hg38)
  blacklist_fraction_low <- median(sample$blacklist_fraction) - 2*mad(sample$blacklist_fraction)
  blacklist_fraction_hi <- median(sample$blacklist_fraction) + 2*mad(sample$blacklist_fraction)

  nucleosome_signal_low <- median(sample$nucleosome_signal) - 2*mad(sample$nucleosome_signal)
  nucleosome_signal_hi <- median(sample$nucleosome_signal) + 2*mad(sample$nucleosome_signal)
  sample$nucleosome_group <- ifelse(sample$nucleosome_signal > nucleosome_signal_hi, paste0("NS > ",nucleosome_signal_hi), paste0("NS < ",nucleosome_signal_hi))
  p5 <- FragmentHistogram(sample, group.by='nucleosome_group') + NoLegend() +theme_bw()
  ggsave(paste0(output_path,"qc_pre_filtering_nucleosome_group.jpeg"),device="jpeg",plot=p5,height=5,width=7)

  p6 <- VlnPlot(
    object = sample,
    features = c('pct_reads_in_peaks', 'atac_peak_region_fragments',
                 'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal'),
    pt.size = 0.1,
    ncol = 5
  )
  ggsave(paste0(output_path,"qc_pre_filtering_peak_info.jpeg"),device="jpeg",plot=p6,height=5,width=15)

  message("nCount_ATAC summary")
  print(summary(sample$nCount_ATAC))
  boxplot(sample$nCount_ATAC)

  nATAC_low <- median(sample$nCount_ATAC)-2*mad(sample$nCount_ATAC)
  nATAC_high <- median(sample$nCount_ATAC)+2*mad(sample$nCount_ATAC)

  message("nFeature_RNA summary")
  print(summary(sample$nFeature_RNA))
  hist(sample$nFeature_RNA)
  boxplot(sample$nFeature_RNA)
  nGene.low <- median(sample$nFeature_RNA)-2*mad(sample$nFeature_RNA)
  nGene.high <- median(sample$nFeature_RNA)+2*mad(sample$nFeature_RNA)

  message("nCount_RNA summary")
  nUMI.low <- median(sample$nCount_RNA)-2*mad(sample$nCount_RNA)
  nUMI.high <-  median(sample$nCount_RNA)+2*mad(sample$nCount_RNA)

  message("Percent.mito summary")
  mito.low <- 0
  mito.high <- 10

  sample.filtered <- subset(sample,
                            subset =
                              pct_reads_in_peaks > pct_reads_in_peaks_low &
                              pct_reads_in_peaks < pct_reads_in_peaks_hi &
                              atac_peak_region_fragments > atac_peak_region_fragments_low &
                              atac_peak_region_fragments < atac_peak_region_fragments_hi &
                              TSS.enrichment > TSS_enrichment_low &
                              TSS.enrichment < TSS_enrichment_hi &
                              blacklist_fraction > blacklist_fraction_low &
                              blacklist_fraction < blacklist_fraction_hi &
                              nucleosome_signal > nucleosome_signal_low &
                              nucleosome_signal < nucleosome_signal_hi &
                              nCount_ATAC > nATAC_low &
                              nCount_ATAC < nATAC_high &
                              nFeature_RNA>nGene.low &
                              nFeature_RNA<nGene.high &
                              nCount_RNA > nUMI.low &
                              nCount_RNA < nUMI.high &
                              percent.mt>mito.low &
                              percent.mt<mito.high)

  p1 <- VlnPlot(sample.filtered, features = c("nCount_ATAC","nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 4)+NoLegend()
  ggsave(paste0(output_path,"qc_post_filtering.jpeg"),device="jpeg",plot=p1,height=5,width=10)
  p2 <- FeatureScatter(sample.filtered, feature1 = "nCount_RNA", feature2 = "percent.mt")+theme_bw()+NoLegend()
  ggsave(paste0(output_path,"qc_post_filtering_nCount_v_mt.jpeg"),device="jpeg",plot=p2,height=5,width=7)
  p3 <- FeatureScatter(sample.filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+theme_bw()+NoLegend()
  ggsave(paste0(output_path,"qc_post_filtering_nCount_v_nFeature.jpeg"),device="jpeg",plot=p3,height=5,width=7)
  p3_2 <- FeatureScatter(sample.filtered, feature1 = "nCount_RNA", feature2 = "nCount_ATAC")+theme_bw()+NoLegend()
  ggsave(paste0(output_path,"qc_post_filtering_nRNA_v_nATAC.jpeg"),device="jpeg",plot=p3_2,height=5,width=7)
  p4 <- TSSPlot(sample.filtered, group.by='high.tss') + NoLegend() +theme_bw()
  ggsave(paste0(output_path,"qc_post_filtering_tss_plot.jpeg"),device="jpeg",plot=p4,height=5,width=7)
  p5 <- FragmentHistogram(sample.filtered, group.by='nucleosome_group') + NoLegend() +theme_bw()
  ggsave(paste0(output_path,"qc_post_filtering_nucleosome_group.jpeg"),device="jpeg",plot=p5,height=5,width=7)
  p6 <- VlnPlot(
    object = sample.filtered,
    features = c('pct_reads_in_peaks', 'atac_peak_region_fragments',
                 'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal'),
    pt.size = 0.1,
    ncol = 5
  )
  ggsave(paste0(output_path,"qc_post_filtering_peak_info.jpeg"),device="jpeg",plot=p6,height=5,width=10)

  return(sample.filtered)
})
#Save filtered object under new key
M12f <- sample.list.filtered$M12sample

#### Azimuth ####
#Save RNA counts matrix
rna.matrix <- GetAssayData(object = M12f,slot = 'counts',assay = 'RNA')
saveRDS(rna.matrix,paste0(output_path,"RNA_data_for_azimuth.rds"))

#Run Azimuth via web app, bone marrow reference#
#Load Azimuth predictions
predictions <- read.delim(paste0(output_path,'azimuth/azimuth_pred.tsv'), row.names = 1)
M12f <- AddMetaData(
  object = M12f,
  metadata = predictions)

####GET PEAKS-------------------------
peaks <- CallPeaks(object = M12f,
                   group.by = "predicted.celltype.l2",
                   macs2.path = "/gpfs/commons/home/eeton/.local/share/r-miniconda/envs/PeakCalls_analysis/bin/macs2",
                   assay = "ATAC")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

DefaultAssay(M12f) <- "ATAC"
# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(M12f),
  features = peaks,
  cells = colnames(M12f)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
M12f[["Peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = frag.file,
  annotation = annotations
)
saveRDS(object = M12f,file = paste0(output_path,"M12f_pre_analysis.rds"))

#### RNA analysis ####
DefaultAssay(M12f) <- "RNA"
Idents(M12f)<-"orig.ident"

###Normalize the RNA data
message("Normalizing RNA")
M12f <- NormalizeData(M12f,normalization.method="LogNormalize",scale.factor=10000)
M12f <- FindVariableFeatures(M12f, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(M12f), 10)

###Find and plot variable features with and without labels
message("Variable Features")
plot1 <- VariableFeaturePlot(M12f)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)+theme_bw()
ggsave(paste0(output_path,"highvarfx.jpeg"),plot = plot2,device="jpeg",height=5,width=5)

### Scaling the data
message("Scaling data")
all.genes <- rownames(M12f)
M12f <- ScaleData(M12f,features=all.genes)

###Perform linear dimensional reduction
message("Running PCA")
M12f <- RunPCA(M12f,
               features=VariableFeatures(object=M12f),
               reduction.name = "pca")

###Examine and visualize PCA results in a few different ways
message("Visualizing PCA")
print(M12f[["pca"]],dims=1:20,nfeatures=5)
plot3 <- VizDimLoadings(M12f,dims=1:5,reduction="pca")
ggsave(paste0(output_path,"pca1.jpeg"),plot = plot3,device="jpeg",height=5,width=10)

plot4 <- DimPlot(M12f,reduction="pca")+NoLegend()
ggsave(paste0(output_path,"pca2.jpeg"),plot = plot4,device="jpeg",height=5,width=7)

jpeg(file=paste0(output_path,"pca3.jpeg"))
plot5 <- DimHeatmap(M12f,dims=1:20,cells=500,balanced=TRUE,reduction="pca")
dev.off()

###Determine the dimensionality of the dataset
message("JackStraw")
M12f<-JackStraw(M12f,num.replicate=100,dims=50)
M12f<-ScoreJackStraw(M12f,dims=1:50)
plot6<-JackStrawPlot(M12f,dims=1:50)+theme_bw()
ggsave(paste0(output_path,"jackstraw.jpeg"),plot = plot6,device="jpeg",height=5,width=7)
plot7<-ElbowPlot(M12f,ndims=50)+theme_bw()
ggsave(paste0(output_path,"elbowplot.jpeg"),plot = plot7,device="jpeg",height=5,width=7)
saveRDS(M12f,"intermediate_M12f.rds")

M12f <- readRDS("~/05_M12/06_integrated/intermediate_M12f.rds")
DefaultAssay(M12f)<-'RNA'

### Cluster the cells
message("RNA clustering")
final_dim <- min(which(M12f@reductions$pca@jackstraw$overall.p.values[,2] > 0.05))-1
M12f <- FindNeighbors(M12f,dims=1:final_dim)
M12f <- FindClusters(M12f,resolution = 1)

### Run non-linear dimensional reduction (UMAP)
message("RNA UMAPping")
M12f <- RunUMAP(M12f,dims=1:final_dim,reduction.name="rna_umap")
plot8<-DimPlot(M12f,reduction="rna_umap")
ggsave(paste0(output_path,"RNA_umap.jpeg"),plot = plot8,device="jpeg",height=5,width=7)

#### ATAC analysis ####
message("Starting ATAC analysis")
DefaultAssay(M12f)<-"Peaks"
M12f <- RunTFIDF(M12f)
M12f <- FindTopFeatures(M12f,min.cutoff = "q0")
M12f <- RunSVD(M12f)
plot10<-DepthCor(M12f,n = 30)
ggsave(paste0(output_path,"atac_depth_cor.jpeg"),plot = plot10,device="jpeg",height=5,width=7)
plot10_5 <- ElbowPlot(M12f,reduction='lsi',ndims = 30)
ggsave(paste0(output_path,"atac_stdev_elbow_plot.jpeg"),plot = plot10_5,device="jpeg",height=5,width=7)
message("ATAC UMAPping")

M12f <- RunUMAP(M12f, reduction = "lsi", dims = 2:15, reduction.name = "atac_umap", reduction.key="atacUMAP_")
M12f <- FindNeighbors(object = M12f, reduction = 'lsi')
M12f <- FindClusters(object=M12f,verbose=FALSE,algorithm=3,resolution = 1)

plot11<-DimPlot(object = M12f,label = T,reduction="atac_umap")+NoLegend()
ggsave(paste0(output_path,"atac_umap.jpeg"),plot = plot11,device="jpeg",height=5,width=7)

#### WNN vignette ####
message("Computing weighted nearest neighbors")
M12f <- FindMultiModalNeighbors(M12f, reduction.list = list("pca", "lsi"), dims.list = list(1:38, 2:15))
M12f <- RunUMAP(M12f, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
M12f <- FindClusters(M12f, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
saveRDS(M12f,paste0(output_path,"almost_final.rds"))
print(colnames(M12f@meta.data))
plot12<-DimPlot(M12f,reduction = "wnn.umap",group.by='wsnn_res.0.8')
ggsave(paste0(output_path,"wnn_umap.jpeg"),plot = plot12,device="jpeg",height=5,width=7)

####Visualize cells on UMAP, labeled by Azimuth####
#UMAP
plot13 <- DimPlot(M12f,reduction="wnn.umap",group.by='predicted.celltype.l2',label=T,shuffle=T,repel=T)
ggsave(paste0(output_path,"wnn_umap_azi_labels.jpeg"),plot = plot13,device="jpeg",height=5,width=7)

#Barplot frequency of celltypes
celltypes <- data.frame(table(M12f$predicted.celltype.l2),row.names = 'Var1')
plot14 <- ggplot(data=celltypes,aes(x=rownames(celltypes), y=Freq))+
  geom_bar(stat='identity')+
  coord_flip()+
  theme_bw()
ggsave(paste0(output_path,"celltypes_barplot.jpeg"),plot = plot14,device="jpeg",height=5,width=7)

#Display UMAP of seurat with filtered azimuth labels
to_remove_celltypes <- M12f@meta.data %>% 
  count(predicted.celltype.l2) %>% 
  dplyr::filter(n < 50) %>% 
  pull(predicted.celltype.l2)
M12f$predicted.l2.plot <- if_else(M12f$predicted.celltype.l2 %in% 
                                    to_remove_celltypes, "Other", M12f$predicted.celltype.l2)
M12f_filtered_azimuth<-DimPlot(M12f, 
                               reduction = 'wnn.umap',
                               group.by = "predicted.l2.plot", 
                               shuffle = T, 
                               label = T, 
                               label.box = T,
                               repel=T) + NoLegend()
ggsave(paste0(output_path,"azimuth_umap_filtered.jpeg"),plot=M12f_filtered_azimuth,device="jpeg",width=7,height=5)

####Load M12_dogma with ATAC,GeneScores,Peaks
message("Loading M12 dogma ATAC")
old_M12<-readRDS("~/05_M12/01_subset/M12_dogma.rds")
mito_calls_colnames <- grepl("threshold_",colnames(old_M12@meta.data))
mito_calls <- old_M12@meta.data[,mito_calls_colnames]
rownames(mito_calls)<-gsub("M12_dogma_","",rownames(mito_calls))
M12f<-AddMetaData(M12f,mito_calls)

M12f_mutant_plot<-DimPlot(M12f, 
                     reduction = 'wnn.umap',
                     group.by = "threshold_10_reads", 
                     shuffle = T, 
                     label = F, 
                     label.box = T,
                     repel=T)
ggsave(paste0(output_path,"threshold_10_reads.jpeg"),plot=M12f_mutant_plot,device="jpeg",width=7,height=5)

M12f$Sample <- rep("M12",length(Cells(M12f)))
M12f$Method <- rep("dogma",length(Cells(M12f)))

#### Save 
message("Saving RNA+ATAC+Peaks object")

saveRDS(M12f,paste0(output_path,"M12_final_ATAC_RNA_Peaks.rds"))
