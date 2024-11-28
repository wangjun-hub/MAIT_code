####scRNA-seq
library(Seurat)
library(SeuratData)
library(paletteer) 
library(scales)
library(ggplot2)
library(harmony) 
library(celldex)
library(dplyr)
library(patchwork)
library(patchwork)
library(rtracklayer)

data_total<-list()

new_counts1 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD2_LC1_5/filtered_feature_bc_matrix")
new_counts2 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD2_LC2_5/filtered_feature_bc_matrix")
new_counts3 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD1_LC1_5/filtered_feature_bc_matrix")
new_counts4 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD1_LC2_5/filtered_feature_bc_matrix")
new_counts5 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD1_LC3_5/filtered_feature_bc_matrix")
new_counts6 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD4_LC1_5/filtered_feature_bc_matrix")
new_counts7 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD4_LC2_5/filtered_feature_bc_matrix")
new_counts8 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD5_LC1_5/filtered_feature_bc_matrix")
new_counts9 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD5_LC2_5/filtered_feature_bc_matrix")
new_counts10 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD9_lc1/filtered_feature_bc_matrix")
new_counts11 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD9_lc2/filtered_feature_bc_matrix")
new_counts12 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD8_LC1_5/filtered_feature_bc_matrix")
new_counts13 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD8_LC2_5/filtered_feature_bc_matrix")
new_counts14 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD14_lc1/filtered_feature_bc_matrix")
new_counts15 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD14_lc2/filtered_feature_bc_matrix")

new_counts16 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD16_lc1/filtered_feature_bc_matrix")
new_counts17 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD16_lc2/filtered_feature_bc_matrix")


data_total[[1]] <- CreateSeuratObject(counts = new_counts1, min.cells = 3,min.features = 200, project = "T1")
data_total[[2]] <- CreateSeuratObject(counts = new_counts2, min.cells = 3,min.features = 200, project = "T2")
data_total[[3]] <- CreateSeuratObject(counts = new_counts3, min.cells = 3,min.features = 200, project = "T3")
data_total[[4]] <- CreateSeuratObject(counts = new_counts4, min.cells = 3,min.features = 200, project = "T4")
data_total[[5]] <- CreateSeuratObject(counts = new_counts5, min.cells = 3,min.features = 200, project = "T5")
data_total[[6]] <- CreateSeuratObject(counts = new_counts6, min.cells = 3,min.features = 200, project = "T6")
data_total[[7]] <- CreateSeuratObject(counts = new_counts7, min.cells = 3,min.features = 200, project = "T7")
data_total[[8]] <- CreateSeuratObject(counts = new_counts8, min.cells = 3,min.features = 200, project = "T8")
data_total[[9]] <- CreateSeuratObject(counts = new_counts9, min.cells = 3,min.features = 200, project = "T9")
data_total[[10]] <- CreateSeuratObject(counts = new_counts10, min.cells = 3,min.features = 200, project = "T10")
data_total[[11]] <- CreateSeuratObject(counts = new_counts11, min.cells = 3,min.features = 200, project = "T11")
data_total[[12]] <- CreateSeuratObject(counts = new_counts12, min.cells = 3,min.features = 200, project = "T12")
data_total[[13]] <- CreateSeuratObject(counts = new_counts13, min.cells = 3,min.features = 200, project = "T13")
data_total[[14]] <- CreateSeuratObject(counts = new_counts14, min.cells = 3,min.features = 200, project = "T14")
data_total[[15]] <- CreateSeuratObject(counts = new_counts15, min.cells = 3,min.features = 200, project = "T15")

data_total[[16]] <- CreateSeuratObject(counts = new_counts16, min.cells = 3,min.features = 200, project = "T16")
data_total[[17]] <- CreateSeuratObject(counts = new_counts17, min.cells = 3,min.features = 200, project = "T17")



###Quality Control
#Homo_sapiens.GRCh38.103.chr.gtf.gz
#Loading Annotations
gtf <- import("/home/shuishan/software/genome/Homo_sapiens.GRCh38.104.chr.gtf.gz")

gtf <- gtf[!is.na(gtf$gene_name)]
gtf <- gtf[gtf$gene_name!=""]

## protein coding
protein <- 
  gtf$gene_name[gtf$transcript_biotype %in% 
                  c("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene", 
                    "IG_M_gene", "IG_V_gene", "IG_Z_gene", 
                    "nonsense_mediated_decay", "nontranslating_CDS", 
                    "non_stop_decay", "polymorphic_pseudogene", 
                    "protein_coding", "TR_C_gene", "TR_D_gene", "TR_gene", 
                    "TR_J_gene", "TR_V_gene")]

## mitochondrial genes
mito <- gtf$gene_name[as.character(seqnames(gtf)) %in% "MT"]

## long noncoding
lincRNA <- 
  gtf$gene_name[gtf$transcript_biotype %in% 
                  c("3prime_overlapping_ncrna", "ambiguous_orf", 
                    "antisense_RNA", "lincRNA", "ncrna_host", "non_coding", 
                    "processed_transcript", "retained_intron", 
                    "sense_intronic", "sense_overlapping")]

## short noncoding
sncRNA <- 
  gtf$gene_name[gtf$transcript_biotype %in% 
                  c("miRNA", "miRNA_pseudogene", "misc_RNA", 
                    "misc_RNA_pseudogene", "Mt_rRNA", "Mt_tRNA", 
                    "Mt_tRNA_pseudogene", "ncRNA", "pre_miRNA", 
                    "RNase_MRP_RNA", "RNase_P_RNA", "rRNA", "rRNA_pseudogene", 
                    "scRNA_pseudogene", "snlRNA", "snoRNA", 
                    "snRNA_pseudogene", "SRP_RNA", "tmRNA", "tRNA",
                    "tRNA_pseudogene", "ribozyme", "scaRNA", "sRNA")]

## pseudogene
pseudogene <- 
  gtf$gene_name[gtf$transcript_biotype %in% 
                  c("disrupted_domain", "IG_C_pseudogene", "IG_J_pseudogene", 
                    "IG_pseudogene", "IG_V_pseudogene", "processed_pseudogene", 
                    "pseudogene", "transcribed_processed_pseudogene",
                    "transcribed_unprocessed_pseudogene", 
                    "translated_processed_pseudogene", 
                    "translated_unprocessed_pseudogene", "TR_J_pseudogene", 
                    "TR_V_pseudogene", "unitary_pseudogene", 
                    "unprocessed_pseudogene")]

annotations <- list(protein=unique(protein), 
                    mito=unique(mito),
                    lincRNA=unique(lincRNA),
                    sncRNA=unique(sncRNA),
                    pseudogene=unique(pseudogene))


names(annotations )
annotations <- annotations[lengths(annotations)>0]

#Quality Control and Cell Selection

for (i in 1:17){

  data_total[[i]][["percent.mt"]] <- PercentageFeatureSet(data_total[[i]], pattern = "^MT-")
  
  # Visualize QC metrics as a violin plot
  plot=VlnPlot(data_total[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  #plot
  ggsave(paste0("./result_both/",i,"_nFeature_RNA & nCount_RNA & percent.mt.png"), plot, width=20 ,height=10)



  percent <- lapply(annotations, function(.ele){
    PercentageFeatureSet(data_total[[i]], features = rownames(data_total[[i]])[rownames(data_total[[i]]) %in% .ele])
  })

  
  for(j in seq_along(percent)){
    data_total[[i]][[paste0("percent.", names(percent)[j])]] <- percent[[j]]
  }

  
  plot1 <- VlnPlot(object = data_total[[i]], 
                   features = c("nFeature_RNA", "nCount_RNA",
                                paste0("percent.", names(percent))), 
                   ncol = 4)
  #plot1
  ggsave(paste0("./result_both/",i,"_annotations_information.png"), plot1, width=20 ,height=20)

  plot2 <- FeatureScatter(data_total[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot3 <- FeatureScatter(data_total[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot4 <- plot2 + plot3
  #plot4
  ggsave(paste0("./result_both/",i,"_correlation.png"), plot4, width=20 ,height=10)


  data_total[[i]] <- subset(x = data_total[[i]], 
                 subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &percent.protein > 80 & 
                 percent.mt < 20)
}



#Harmony Algorithm
scRNA<-merge(data_total[[1]],data_total[2:length(data_total)])
table(scRNA$orig.ident)

#Run PCA
scRNA <- NormalizeData(scRNA)

scRNA <- FindVariableFeatures(scRNA, selection.method = "vst") %>% ScaleData()

scRNA <- RunPCA(scRNA, verbose = F)


cellinfo<-subset(scRNA@meta.data,select=c("orig.ident","percent.mt","nFeature_RNA","nCount_RNA","percent.protein"))

scRNA<-CreateSeuratObject(scRNA@assays$RNA@counts,meta.data=cellinfo)

scRNA <- SCTransform(scRNA)

scRNA <- RunHarmony(scRNA,group.by.vars="orig.ident",assay.use="SCT",max.iter.harmony=20)


saveRDS(scRNA,file="./scRNA.rds")

