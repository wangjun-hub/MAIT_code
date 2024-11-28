####Figure 4B

####Visual receptor analysis
library(devtools)
library(ktplots)

setwd("/Users/wangjun/Desktop/out")

mypvals <- read.delim("pvalues.txt", check.names = FALSE)
mymeans <- read.delim("means.txt", check.names = FALSE)

colnames(mypvals)

mymeans[grepl("IL18",mymeans$interacting_pair),grepl("MCs",colnames(mymeans))]

mypvals[grepl("IL18",mypvals$interacting_pair),grepl("MCs",colnames(mypvals))]

mypvals[grepl("|C2: C1q",mypvals$interacting_pair),1:6]

kp = grep("MCs", colnames(mypvals))


table(kp)
pos = (1:ncol(mypvals))[kp] 
choose_pvalues <- mypvals[,c(c(1,2,5,6,8,9),pos  )]
choose_means <- mymeans[,c(c(1,2,5,6,8,9),pos)]

dim(choose_pvalues)


logi <- apply(choose_pvalues[,5:ncol(choose_pvalues)]<0.01, 1, sum) 
# Only retaining some cell specific interactions
choose_pvalues <- choose_pvalues[logi>=2,]

dim(choose_pvalues)

choose_pvalues[,1:4]



# Remove null values
logi1 <- choose_pvalues$gene_a != "ddd"
logi2 <- choose_pvalues$gene_b != "ddd"
logi <- logi1 & logi2
choose_pvalues <- choose_pvalues[logi,]

max(choose_means[,7:17])


# Retain choose-means under the same conditions
choose_means <- choose_means[choose_means$id_cp_interaction %in% choose_pvalues$id_cp_interaction,]

choose_means[,1:4]
dim(choose_means)
colnames(choose_pvalues)

max_values <- as.numeric(apply(choose_means[7:17], 1, max))>0.5


choose_means <- choose_means[,c(1:6,25,20:24,26,28,29,7,27)]
choose_pvalues <- choose_pvalues[,c(1:6,25,20:24,26,28,29,7,27)]

choose_means <- choose_means[max_values,]
choose_pvalues<-choose_pvalues[max_values,]

choose_pvalues[1:4,]

logi <- apply(choose_pvalues[,5:ncol(choose_pvalues)]<0.01, 1, sum) 
# Only retaining some cell specific interactions
choose_pvalues <- choose_pvalues[logi>=2,]
choose_means <- choose_means[logi>=2,]





dim(choose_means)
  

# Convert choose_pvalues and choose-means data
library(tidyverse)
meansdf <- choose_means %>% reshape2::melt()
meansdf <- data.frame(interacting_pair = meansdf$interacting_pair,
                      CC = meansdf$variable,
                      means = meansdf$value)
pvalsdf <- choose_pvalues %>% reshape2::melt()
pvalsdf <- data.frame(interacting_pair = pvalsdf$interacting_pair,
                      CC = pvalsdf$variable,
                      pvals = pvalsdf$value)

# Merge p-value and mean files
pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
meansdf$joinlab<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
pldf <- merge(pvalsdf,meansdf,by = "joinlab")

# Dotplot visualization
summary((filter(pldf,means >0))$means)
head(pldf)

pcc =  pldf%>% filter(means >0) %>% 
  ggplot(aes(CC.x,interacting_pair.x) )+ 
  geom_point(aes(color=means,size=-log10(pvals+0.0001)) ) +
  scale_size_continuous(range = c(1,5))+
  scale_color_gradient2(high="#E64B35FF",mid = "#F39B7FFF",low ="#4DBBD5FF",midpoint = 0.9  )+ 
  theme_bw()+ 
  # scale_color_manual(values = rainbow(100))+
  theme(axis.text.x = element_text(angle = -45,hjust = 0,vjust = 0.5))
pcc

ggsave(pcc,file="/Users/wangjun/Desktop/pcc.pdf",width = 7.98,height = 5.33,device=cairo_pdf)

write.table(pldf,file="/Users/wangjun/Desktop/cellphonedb_data.csv",sep = ",",row.names = F,col.names = T)


