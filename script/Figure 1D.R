##Figure 1D

mast<-readRDS(file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/lung_mast.rds")

mast_nc<-readRDS(file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/Mast_cells.rds")

Idents(mast_nc)<-"orig.ident"

mast_nc$Sample<-mast_nc$orig.ident

#mast_nc$Sample<-Idents(mast_nc)
mast_nc$Condition<-"Tumor"

mast_total<-merge(mast,y=mast_nc,add.cell.ids=NULL)


mast <- FindVariableFeatures(mast, selection.method = "vst")

f=VariableFeatures(object = mast)
mast <- ScaleData(mast, features = f,verbose = FALSE)
mast <- RunPCA(mast, npcs = 30, verbose = FALSE)

mast <- RunHarmony(mast,reduction = "pca",lambda=1.7,group.by.vars = "Condition",reduction.save = "harmony")

mast <- FindNeighbors(mast, dims = 1:30, reduction = "harmony")
mast <- FindClusters(mast, resolution = 0.3)

p1<-DotPlot(mast, features = c("IL18"), cols = c("green","red"))+scale_size_continuous(range = c(2, 6)) + 
RotatedAxis()


table(Idents(mast))

##进行重命名
#new.cluster.ids <- c("FCER1A+ID2+MCs","CD81+CCL18+MCs","C1q+IL18+MCs", "RPL+HDC+MCs","Septin2+SCL24A3+MCs","HSP90+AREG+MCs","EEF1G+MIF+ MCs","CD74+HLADR+MCs","CCL5+IL32+MCs")
new.cluster.ids <- c("C0: FCER1A+ID2+MCs","C1: CD81+CCL18+MCs","C2: C1q+IL18+MCs", "C3: RPL+HDC+MCs","C4: Septin2+SCL24A3+MCs","C5: HSP90+AREG+MCs","C6: EEF1G+MIF+ MCs","C7: CD74+HLADR+MCs","C8: CCL5+IL32+MCs")


#鉴定出相应的细胞类群后
names(new.cluster.ids) <- levels(mast)
mast <- RenameIdents(mast, new.cluster.ids)

mast$new_name<-Idents(mast)
saveRDS(mast,file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast_6_25.rds")



mast <- readRDS(file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast_6_25.rds")

IL18_mast <- subset(mast,idents = "C2: C1q+IL18+MCs")

saveRDS(IL18_mast,file="/Users/wangjun/Desktop/IL18_mast.rds")


Idents(mast)<-factor(Idents(mast),levels = c("C2: C1q+IL18+MCs","C0: FCER1A+ID2+MCs","C1: CD81+CCL18+MCs", "C3: RPL+HDC+MCs","C4: Septin2+SCL24A3+MCs","C5: HSP90+AREG+MCs","C6: EEF1G+MIF+ MCs","C7: CD74+HLADR+MCs","C8: CCL5+IL32+MCs"))

##绘制聚类图
p2 <- DimPlot(mast,reduction = "tsne",pt.size = 0.3)+
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
p2

ggsave(p2,file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast/mast_tsne_6_30.pdf",width = 5.81,height = 4.15) 
ggsave(p2,file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast/mast_tsne_6_30.png",width = 5.81,height = 4.15 ,dpi=1000)
