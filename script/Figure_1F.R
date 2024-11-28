###Figure 1F
library(clusterProfiler)
library(org.Hs.eg.db)
library(devtools)


marker_gene <- read.csv(file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast/mast_6_25/DEG.csv",sep = ",",header = T,row.names = F)

geneList1<-bitr(unique(marker_gene$gene), fromType = "SYMBOL",
                toType = c( "ENTREZID"),
                OrgDb = org.Hs.eg.db)


marker_gene<-data.frame(marker_gene)

marker_gene_total<-merge(marker_gene,y=geneList1,by.x="gene",by.y="SYMBOL")
marker_gene_total[1:4,]

marker_gene_total<-marker_gene_total[!is.na(marker_gene_total$ENTREZID),]
dim(marker_gene_total)
dim(marker_gene_total[!is.na(marker_gene_total$ENTREZID),])
dim(marker_gene)
write.table(marker_gene_total,file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast/mast_6_25/DEG.csv",sep = ",",col.names = T,row.names = F)

##GO enrichment analysis
for (i in 1:9){
  print(i)
  
  cluster <- names(table(marker_gene$cluster))[i]
  
  marker_gene_cluster1<-marker_gene[marker_gene$cluster%in%cluster,]
  
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
  
  #colnames(df)
  write.table(df2,file = paste0("/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast/mast_6_25/",cluster,"_GO.csv"),sep = ",",row.names=F)
}


##KEGG enrichment analysis
for (i in 1:9){
  print(i)
  cluster <- names(table(marker_gene$cluster))[i]
  
  marker_gene_cluster1<-marker_gene[marker_gene$cluster%in%cluster,]
  
  kk2<-enrichKEGG( marker_gene_cluster1$ENTREZID,
                   organism = "hsa",
                   keyType = "kegg",
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH")
  
  df<-kk2@result
  
  df2<-df[,c(1:7,9)]
  #colnames(df)
  write.table(df2,file = paste0("/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast/mast_6_25/",cluster,"_KEGG.csv"),sep = ",",row.names=F)
}

