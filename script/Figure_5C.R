####Figure 5C
pbmc<-readRDS(file="/Users/wangjun/Desktop/FluT_MAIT/maitdata.rds")


table(Idents(pbmc))

new.cluster.ids <- c("CXCR4+IFNG+ MAIT","NCR3+KLRG1+ MAIT","FCER1G+IFI30+ MAIT","CX3CR1+GNLY+ MAIT","CD24+HSPB1+ MAIT")

names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)


####marker projection 
marker_list<-c("CXCR4", "CD69","KLRB1", "NCR3", "KLRG1", "GZMK", "S1PR1", "IL7R", "FCER1G", "IFI30", "HLA-DRA", "CD74", "NKG7", "KLRD1", "GZMB", "CX3CR1", "GNLY", "HSPB1", "CD9", "CD24", "KLRK1")

for (i in marker_list){
  p9 <- FeaturePlot(pbmc, reduction = "tsne",features = i,cols = c("gray","red"),repel = TRUE,order = T)+ 
  RotatedAxis()+
    theme(
      plot.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 6),
      text = element_text(family = "sans"),
      aspect.ratio = 1,
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.5,
                               arrow = arrow(length = unit(0.1, "cm"), 
                                             ends = "last", type = "closed")),
      axis.title = element_text(hjust = 0, size = 6, face = "bold")) +
    guides(size = 3, x = axis, y = axis) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL)

  ggsave(p9,file=paste0("/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/FeaturePlot_new/MAIT_",i,"_RNA.pdf"),width=4.35,height=4.02,limitsize = FALSE)
  ggsave(p9,file=paste0("/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/FeaturePlot_new/MAIT_",i,"_RNA.png"),width=4.35,height=4.02,dpi=1000,limitsize = FALSE)

}
