####Figure 1E

####Calculate differential genes
scRNA.markers <- FindAllMarkers(object = mast, only.pos = TRUE, 
                                min.pct = 0, thresh.use = 0)

library(dplyr)
marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.05) %>% group_by(cluster) %>% top_n(n = 2000, wt = avg_log2FC)


write.csv(marker_gene,"/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast/mast_6_25/DEG.csv",row.names=F)

####Draw heatmap
marker_gene_heatmap <- scRNA.markers %>% filter(p_val_adj <= 0.05) %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

markers <- marker_gene_heatmap$gene
markers <- as.data.frame(markers)
markerdata <- ScaleData(mast_subset, features = as.character(unique(markers$markers)), assay = "RNA")


p2<-DoHeatmap(markerdata,
          features = as.character(unique(markers$markers)),
          group.by = "RNA_snn_res.0.1",
          assay = 'RNA',
          angle=45,
          group.colors = pal[c(2,1,5,4,7)],size = 0)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))

ggsave(p2,file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast/MAST_cluster_heatmap.pdf", width = 4.85, height = 4.39,device=cairo_pdf)
ggsave(p2,file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast/MAST_cluster_heatmap.png", width = 4.85, height = 4.39,dpi=1000)

