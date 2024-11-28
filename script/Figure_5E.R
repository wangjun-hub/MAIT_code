####Figure 5E
####Calculate differential genes
scRNA.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)

library(dplyr)
marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.05) %>% group_by(cluster) %>% top_n(n = 2000, wt = p_val_adj)



####Draw heatmap
marker_gene_heatmap <- scRNA.markers %>% filter(p_val_adj <= 0.05) %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)



markers <- marker_gene_heatmap$gene
markers <- as.data.frame(markers)
markerdata <- ScaleData(pbmc, features = as.character(unique(markers$markers)), assay = "RNA")


p2<-DoHeatmap(markerdata,
          features = as.character(unique(markers$markers)),
          group.by = "Cluster",
          assay = 'RNA',
          angle=45,
          group.colors = pal,size = 0)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))

ggsave(p2,file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/MAIT_cluster_heatmap.pdf", width = 4.85, height = 4.39,device=cairo_pdf)
ggsave(p2,file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/MAIT_cluster_heatmap.png", width = 4.85, height = 4.39,dpi=1000)



####Pathway enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)


##GO
for (i in 0:4){
  print(i)
  marker_gene_cluster1<-marker_gene[marker_gene$cluster%in%i,]

  geneList1<-bitr(unique(marker_gene_cluster1$gene), fromType = "SYMBOL",
             toType = c( "ENTREZID"),
             OrgDb = org.Hs.eg.db)

  GO2<-enrichGO( geneList1$ENTREZID,
                 OrgDb = org.Hs.eg.db,
                 keyType = "ENTREZID",
                 ont = "ALL",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 readable = TRUE)

  df<-GO2@result

  df2<-df[,c(1,2,3,6,7,9,10)]

  write.table(df2,file = paste0("/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/Cluster_",i,"_GO.csv"),sep = ",",row.names=F)

}


# BiocManager::install("msigdbr")
library(msigdbr)
# Extract C6 library
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

for (i in 0:4){
  print(i)
  marker_gene_cluster1<-marker_gene[marker_gene$cluster%in%i,]

  geneList1<-bitr(unique(marker_gene_cluster1$gene), fromType = "SYMBOL",
             toType = c( "ENTREZID"),
             OrgDb = org.Hs.eg.db)

  em <- enricher(geneList1$ENTREZID, TERM2GENE=m_t2g,readable = TRUE)
  df<-em@result

  write.table(df,file = paste0("/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/Cluster_",i,"_Hallmark.csv"),sep = ",",row.names=F)

}
