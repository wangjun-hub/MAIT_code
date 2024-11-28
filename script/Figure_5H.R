####Figure 5H

pbmc<-readRDS(file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/MAIT_anno_6_22.rds")

table(Idents(pbmc))

DefaultAssay(pbmc)<-"AbC"

marker_list1<-c("CCR4-Ab","CCR5-Ab","CCR6-Ab","CXCR3-Ab","CX3CR1-Ab","IL2RA-Ab","IL4R-Ab","ITGAE-Ab","CD69-Ab","ITGA1-Ab")
marker_list2<-c("NCR1-Ab","NCAM1-Ab","KLRB1-Ab","KLRG1-Ab","KLRK1-Ab","KLRD1-Ab","KIR2DL1-Ab","KIR2DL3-Ab","KIR3DL1-Ab","B3GAT1-Ab")
marker_list3<-c("PDCD1-Ab","TIGIT-Ab","CTLA4-Ab","LAG3-Ab","ENTPD1-Ab","NT5E-Ab","HLA-DRA-Ab","ICOS-Ab","TNFRSF9-Ab","TNFRSF4-Ab")
marker_list4<-c("PTPRC-CD45RA-Ab","PTPRC-CD45RO-Ab","CD44-Ab","SELL-Ab","SELP-Ab","FAS-Ab","DPP4-Ab","CD28-Ab","CD27-Ab","LAMP1-Ab")

rownames(pbmc)

marker_list_total<-c(marker_list1,marker_list2,marker_list3,marker_list4)

df <- data.frame(t(data.frame(GetAssayData(pbmc))[marker_list_total,]))

df[1:4,]
df$Cluster <- "test"

dim(df)
df$Cluster<-Idents(pbmc)
df[1:4,38:41]

grouped_means <- df %>%
  group_by(Cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

grouped_means <- data.frame(grouped_means)

write.table(grouped_means,sep = ",",row.names = F,col.names = T, file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast/mast_6_25/出图文件/protein_data.csv")

grouped_means_new <- grouped_means
for (i in 2:ncol(grouped_means)){
  grouped_means_new[,i] <- (grouped_means[,i]-min(grouped_means[,i]))/(max(grouped_means[,i])-min(grouped_means[,i]))

}
write.table(grouped_means_new,sep = ",",row.names = F,col.names = T, file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast/mast_6_25/出图文件/protein_data_normalization.csv")
