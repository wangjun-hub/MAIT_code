####Figure 5G
####pyscenic 
###-----python script-----

pyscenic grn --seed 21 --num_workers 40 --method grnboost2 --output mait_grn.tsv -t mait.csv hs_hgnc_tfs.txt

pyscenic ctx mait_grn.tsv hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather  hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather --annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl -t mait.csv  --mode "dask_multiprocessing" --output regulons.csv --num_workers 40

pyscenic aucell -t mait.csv regulons.csv -o auc_mtx.csv --num_workers 40


###-----R script-----
###Read the scenic results

data<-read.csv(file="/Users/wangjun/Downloads/auc_mtx.csv",row.names = 1)

data_t <- data.frame(t(data))
data_t[1:4,1:4]
colnames(mait)[1:4]

table(Idents(mait))

data_t$celltype <- Idents(mait)

library(pheatmap)

dim(data_t)

data_seu<-CreateSeuratObject(counts = data)


data_seu <- NormalizeData(data_seu,normalization.method="LogNormalize")

data_seu <- ScaleData(data_seu)
data_seu <- RunPCA(data_seu,features = rownames(data_seu))


data_seu <- AddMetaData(data_seu, metadata = as.character(Idents(mait)),col.name = "celltype")

Idents(data_seu) <- "celltype"
table(Idents(data_seu))

new.cluster.ids <- c("CXCR4+IFNG+ MAIT", "Other MAIT","Other MAIT","Other MAIT","Other MAIT")

#Identify the corresponding cell groups
names(new.cluster.ids) <- levels(data_seu)
data_seu <- RenameIdents(data_seu, new.cluster.ids)

data_seu$Cell_type <- Idents(data_seu)

library(Seurat)
####Calculate differential genes
scRNA.markers <- FindAllMarkers(object = data_seu, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)

library(dplyr)
marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.05,cluster=="CXCR4+IFNG+ MAIT") %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)


write.csv(marker_gene,"/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/all_DEG_gene.csv",row.names=F)

####Draw heatmap
marker_gene_heatmap <- scRNA.markers %>% filter(p_val_adj <= 0.05,cluster=="CXCR4+IFNG+ MAIT") %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)



markers <- marker_gene_heatmap$gene
markers <- as.data.frame(markers)
markerdata <- ScaleData(data_seu, features = as.character(unique(markers$markers)), assay = "RNA")


p2<-DoHeatmap(markerdata,
              features = as.character(unique(markers$markers)),
              group.by = "Cell_type",
              assay = 'RNA',
              angle=45,
              group.colors = pal,size = 0)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))

###save the figure

ggsave(p2,file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast_9_12/scenic.pdf")
ggsave(p2,file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast_9_12/scenic.png",dpi = 1000)

