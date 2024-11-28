####Figure 4A

count_net <- read.delim("../data/count_network.txt", check.names = FALSE)
inter_net <- read.delim("../data/interaction_count.txt", check.names = FALSE)
pvalues <- read.delim("../data/pvalues.txt", check.names = FALSE)
means <- read.delim("../data/means.txt", check.names = FALSE)
sig.means <- read.delim("../data/significant_means.txt", check.names = FALSE)
deconvoluted <- read.delim("../data/deconvoluted.txt", check.names = FALSE)

#Interactive Network Diagram
library(CellChat)
count_inter <- count_net
#count_inter$count <- count_inter$count/100
count_inter$count <- count_inter$count
library(tidyr)
count_inter<-spread(count_inter, TARGET, count)
rownames(count_inter) <- count_inter$SOURCE
count_inter <- count_inter[, -1]
count_inter <- as.matrix(count_inter)


count_inter <- count_inter[,c("CD8+ Exhausted T","CD56brightCD16- NK cells","Treg","CD56dimCD16+ NK cells","MAIT","gd T","C2: C1q+IL18+MCs","Plasma cells","CD8+ GZMK+ T","B cells","Naive T","CD8+ GZMB+ T")]

count_inter <- count_inter[c("CD8+ Exhausted T","CD56brightCD16- NK cells","Treg","CD56dimCD16+ NK cells","MAIT","gd T","C2: C1q+IL18+MCs","Plasma cells","CD8+ GZMK+ T","B cells","Naive T","CD8+ GZMB+ T"),]


node.colors <- c("CD8+ Exhausted T" = "pink" ,"CD56brightCD16- NK cells"="#8491B4FF","Treg"="#F39B7FFF","CD56dimCD16+ NK cells"="#3C5488FF","MAIT" = "#DC0000FF","gd T"="#E64B35FF","C2: C1q+IL18+MCs" = "black","Plasma cells"="#B09C85FF","CD8+ GZMK+ T"="#91D1C2FF","B cells"="#7E6148FF","Naive T"="#4DBBD5FF","CD8+ GZMB+ T"="#00A087FF")

#Set the default edge color to gray
edge.colors <- rep("purple", nrow(count_inter))
#Name each edge
names(edge.colors) <- rownames(count_inter)



library(igraph)
igraph.options(vertex.color = node.colors)

#Draw graphics using netVisual_circle
netVisual_circle(count_inter, 
                 weight.scale = TRUE, 
                 sources.use = "C2: C1q+IL18+MCs",
                 arrow.size = 0.2,
                 color.use = node.colors,
                 vertex.label.cex = 1,
                 vertex.label.color = "black",
                 edge.width.max = 2,
                 alpha.edge = 0.8,
                 edge.curved = TRUE,
                 shape = "circle",
                 layout = layout_in_circle)


pdf("network_7_21_no_number.pdf",height = 6,width = 6)
netVisual_circle(count_inter, 
                   weight.scale = T, 
                   sources.use = "C2: C1q+IL18+MCs",
                   #edge.weight.max = max(mat2), 
                   #label.edge = TRUE,
                   title.name = rownames(count_inter)[i],
                   arrow.size=0.2)

dev.off()
