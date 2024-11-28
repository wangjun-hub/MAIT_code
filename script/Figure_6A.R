####Figure 6A
marker_list<-c("NLRP3", "NLRC4", "NLRP1", "AIM2", "MEFV", "PYCARD", "CASP1", "CASP4", "CASP5", "GSDMA", "GSDMB", "GSDMC", "GSDMD", "GSDME", "IL1B", "IL18")

DefaultAssay(pbmc)<-"RNA"
#DefaultAssay(pbmc)<-"AbC"

df <- data.frame(t(data.frame(GetAssayData(pbmc))[marker_list,]))

df$Cluster <- "test"

df$Cluster<-Idents(pbmc)


grouped_means <- df %>%
  group_by(Cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

grouped_means <- data.frame(grouped_means)

write.table(grouped_means,sep = ",",row.names = F,col.names = T, file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast/mast_6_25/出图文件/RNA_NLRP3_list_data.csv")

grouped_means_new <- grouped_means
for (i in 2:ncol(grouped_means)){
  grouped_means_new[,i] <- (grouped_means[,i]-min(grouped_means[,i]))/(max(grouped_means[,i])-min(grouped_means[,i]))
  
}
write.table(grouped_means_new,sep = ",",row.names = F,col.names = T, file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast/mast_6_25/出图文件/RNA_NLRP3_list_data_normalization.csv")

