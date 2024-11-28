####Figure 5J
####Statistically analyze the diversity and cloning of TCR

##Read data
data_multi <- read.csv(file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/材料补充/Multi_tcr_total_new.csv",sep = ",",header = T)

data_unique <- read.csv(file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/材料补充/Unique_tcr_total_new.csv",sep = ",",header = T)

data_multi_subset <- data_multi[data_multi$v_gene=="TRAV1-2",]

data_multi_subset <- data_multi[data_multi$SampleBarcode%in%data_multi_subset$SampleBarcode & data_multi$chain=="TRB",]
table(data_multi_subset$chain)


data_unique_subset <- data_unique[data_unique$TRA_v_gene=="TRAV1-2",]

dim(data_multi_subset)
dim(data_unique_subset)


#data_multi_subset<-data_multi
#data_unique_subset<-data_unique
data_multi_subset <- data_multi_subset[,c("SampleBarcode","cdr3","orig.ident","Patient","Cluster")]

data_unique_subset <- data_unique_subset[,c("SampleBarcode","TRB_cdr3","orig.ident","Patient","Cluster")]

colnames(data_multi_subset) <- c("SampleBarcode","TRB","orig.ident","Patient","Cluster")
colnames(data_unique_subset) <- c("SampleBarcode","TRB","orig.ident","Patient","Cluster")

data_subset_total <- rbind(data_multi_subset,data_unique_subset)
dim(data_subset_total)

data_subset_total <- unique(data_subset_total)

df <- data_subset_total


#####
df[df$Cluster!="CXCR4+IFNG+ MAIT",]$Cluster <- "Other MAIT"

table(df$Cluster)


###Statistical abundance
df_freq <- df %>%
  group_by(orig.ident, Cluster, raw_clonotype_id) %>%
  summarise(freq = n()) %>%
  ungroup()

# Divide the frequency into different levels 
# (assuming 1, 2,>=3 levels are used here)
df_freq <- df_freq %>%
  mutate(level = cut(freq, 
                     breaks = c(-Inf,  2, Inf), 
                     labels = c("1-2", ">=3")))

# Calculate the relative abundance of each level in each cluster
df_abundance <- df_freq %>%
  group_by(Cluster, freq) %>%
  summarise(weighted_count = sum(freq)) %>%  #Calculate the weighted total based on Freq
  mutate(rel_abundance = weighted_count / sum(weighted_count)) %>%  #Calculate relative abundance
  ungroup()



write.table(df_abundance,file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast_mait/df_abundance_new.csv",sep=",",col.names=T,row.names=F)


# Draw a stacked graph grouped by Cluster
p <- ggplot(df_abundance, aes(x = Cluster, y = rel_abundance, fill = level)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Cluster", y = "Relative Abundance", fill = "Level") +
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = c("1-2" = "#F39B7FFF", ">=3" = "#DC0000FF"))


# show img
print(p)

ggsave(p,file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast_mait/TCR_clonotype_new4.pdf")
ggsave(p,file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast_mait/TCR_clonoytype_new4.png",dpi = 1000)



###Figure 5K

library(vegan)

# Step 1: Group according to Orig.ident and Cluster
df_grouped <- df %>%
  group_by(orig.ident, Cluster)

# Step 2: Calculate the Shannon diversity index of TRB
df_shannon <- df_grouped %>%
  summarise(shannon_index = diversity(table(TRB), index = "shannon")) %>%
  ungroup()

# Step 3: Draw a box plot and add scatter points
p <- ggplot(df_shannon, aes(x = Cluster, y = shannon_index)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue") +  
  geom_jitter(width = 0.2, size = 2, color = "black") + 
  labs(x = "Cluster", y = "Shannon Diversity Index", title = "Shannon Diversity Index by Cluster") +
  theme_minimal()

ggsave(p,file="/Users/wangjun/Desktop/FluT_MAIT/MAIT_final/mast_mait/shannon.png",dpi = 1000)

