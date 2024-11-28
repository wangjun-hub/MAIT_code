####Figure 2D

data_signature<-read.csv(file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast/inflammaseome_complex.csv",sep=",",header=T)

colnames(mast@meta.data)
mast@meta.data <- mast@meta.data[,1:130]


for (i in 1:ncol(data_signature)){
  print(i)
  data_geneset<- as.list(data_signature[,i][data_signature[,i]!=""])
  mast<-AddModuleScore(object=mast,features=data_geneset,ctrl = 100,name = colnames(data_signature)[i])
}

colnames(mast@meta.data)


mast$negative_regulation_of_NLRP3_inflammasome_complex_assembly<-(rowMeans(mast@meta.data[132:138])-min(rowMeans(mast@meta.data[132:138])))/(max(rowMeans(mast@meta.data[132:138]))-min(rowMeans(mast@meta.data[132:138])))

mast$positive_regulation_of_NLRP3_inflammasome_complex_assembly<-(rowMeans(mast@meta.data[139:147])-min(rowMeans(mast@meta.data[139:147])))/(max(rowMeans(mast@meta.data[139:147]))-min(rowMeans(mast@meta.data[139:147])))

mast$regulation_of_NLRP3_inflammasome_complex_assembly<-(rowMeans(mast@meta.data[148:164])-min(rowMeans(mast@meta.data[148:164])))/(max(rowMeans(mast@meta.data[148:164]))-min(rowMeans(mast@meta.data[148:164])))

mast$inflammasome_complex_score<-(rowMeans(mast@meta.data[165:180])-min(rowMeans(mast@meta.data[165:180])))/(max(rowMeans(mast@meta.data[165:180]))-min(rowMeans(mast@meta.data[165:180])))


mast_pathway_plot<-data.frame(Idents(mast),mast$negative_regulation_of_NLRP3_inflammasome_complex_assembly,mast$positive_regulation_of_NLRP3_inflammasome_complex_assembly,mast$regulation_of_NLRP3_inflammasome_complex_assembly,mast$inflammasome_complex_score)
colnames(mast_pathway_plot)<-c("Cluster","Pathway_1","Pathway_2","Pathway_3","inflammasome_pathway_score")

mast_pathway_plot[1:4,]
data <- mast_pathway_plot[,c("Cluster","Pathway_2")]

data[1:4,]

###Remove Outliers
# Calculate quartiles
Q1 <- quantile(data$Pathway_2, 0.25)
Q3 <- quantile(data$Pathway_2, 0.75)
IQR <- Q3 - Q1

# Calculate the upper and lower limits of outliers
lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR

# Mark Outliers
data$outlier <- data$Pathway_2 < lower_bound | data$Pathway_2 > upper_bound

filtered_data <- data[!data$outlier, ]
filtered_data[1:4,]

data<-filtered_data[,c("Cluster","Pathway_2")]


colnames(data)<-c("group","value")

table(data$bin)
data[1:4,]
# Create Bins
data$bin <- cut(data$value, breaks = seq(0, 1, by = 0.01), include.lowest = TRUE)

# Calculate the number of groups in each bin
data_count <- data %>%
  group_by(bin, group) %>%
  summarise(count = n()) %>%
  ungroup()

tail(data_count)

# Calculate the proportion of groups in each bin
data_percentage <- data_count %>%
  group_by(bin) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()


# Use ggplot2 to draw a stacked graph
p1<-ggplot(data_percentage, aes(x = bin, y = percentage, fill = group)) + 
  geom_bar(stat = "identity", position = "fill",width=1.01) +
  labs(x = "X Axis (binned)", y = "Percentage", title = "inflammasome complex assembly") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c(pal[c(8,3,5,9,4,1,6,7,2)]))+
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, title.theme = element_text(size = 12))) +
  
  theme(plot.title = element_text(hjust = 0.5),axis.ticks.x = element_blank(),panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),axis.text.x.bottom = element_text(angle=60,size=6, hjust = 1))+
  
  theme(legend.key = element_rect(fill = "transparent"),legend.background = element_rect(fill = "white"),axis.title.x = element_text(size = 0),axis.title.y = element_text(size = 12),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black",linewidth = 0),legend.box.background=element_rect(colour = "white",fill = "white"))

p1


##Read data2 data

mast_pathway_plot[1:4,]

mast_pathway_plot$IL18<-"test"
IL18_exp <- t(data.frame(GetAssayData(mast)))[,"IL18"]
IL18_exp[1:4]
mast_pathway_plot$IL18 <- IL18_exp


grouped_means <- mast_pathway_plot %>%
  group_by(Cluster) %>%
  summarise(
    mean_Pathway_2 = mean(Pathway_2, na.rm = TRUE),
    mean_inflammasome_score = mean(inflammasome_pathway_score, na.rm = TRUE),
    mean_IL18 = mean(IL18, na.rm = TRUE)
  )


grouped_means[1:4,]

p2<-ggplot(grouped_means, aes(x = mean_Pathway_2, y = mean_inflammasome_score, color = Cluster,size = mean_IL18)) + 
  geom_point() + 
  geom_text_repel(aes(label = Cluster),size=3,force = 3) + 
  xlim(0.18,0.3)+
  ylim(0.2,0.32)+
  guides(color = "none",size = "legend")+
  labs(x = "Positive_Pathway", y = "Inflammasome_component") +
  scale_color_manual(values = c(pal[c(8,3,5,9,4,1,6,7,2)]))+ 
  scale_size_continuous(name = "Mean_IL18") + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),axis.text.x.bottom = element_text(size=12, hjust = 0.5))+
  theme(legend.key = element_rect(fill = "transparent"),legend.background = element_rect(fill = "white"),axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black",linewidth = 0),legend.box.background=element_rect(colour = "white",fill = "white"))


library(patchwork)
p3<-p1/p2
p3


ggsave(p3,file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast/mast_6_25/Figure_5.pdf",width = 8.86,height = 4.93)
ggsave(p3,file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast/mast_6_25/Figure_5.png",dpi = 1000)

