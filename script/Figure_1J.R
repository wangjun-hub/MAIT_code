#####Figure 1J
##Read file
data <- read.csv(file = "/Users/wangjun/Desktop/mast_function.csv", sep = ",", row.names = 1,header = F, fill = T)
data[1:4,]

data <- data.frame(t(data.frame(data)))

colnames(mast@meta.data)

data[,i][data[,i]!=""]

for (i in 1:ncol(data)){
  print(i)
  data_geneset<- as.list(data[,i][data[,i]!=""])
  #print(data_geneset)
  mast<-AddModuleScore(object=mast,features=data_geneset,ctrl = 100,name = colnames(data)[i])
}

colnames(mast@meta.data)

mast$mast.cell.degranulation<-(rowMeans(mast@meta.data[131:138])-min(rowMeans(mast@meta.data[131:138])))/(max(rowMeans(mast@meta.data[131:138]))-min(rowMeans(mast@meta.data[131:138])))

mast$mast.cell.mediated.immunity<-(rowMeans(mast@meta.data[139:147])-min(rowMeans(mast@meta.data[139:147])))/(max(rowMeans(mast@meta.data[139:147]))-min(rowMeans(mast@meta.data[139:147])))

mast$positive.regulation.of.mast.cell.activation<-(rowMeans(mast@meta.data[148:154])-min(rowMeans(mast@meta.data[148:154])))/(max(rowMeans(mast@meta.data[148:154]))-min(rowMeans(mast@meta.data[148:154])))

mast$positive.regulation.of.mast.cell.activation.involved.in.immune.response<-(rowMeans(mast@meta.data[155:159])-min(rowMeans(mast@meta.data[155:159])))/(max(rowMeans(mast@meta.data[155:159]))-min(rowMeans(mast@meta.data[155:159])))

mast$positive.regulation.of.mast.cell.chemotaxis<-(rowMeans(mast@meta.data[160:164])-min(rowMeans(mast@meta.data[160:164])))/(max(rowMeans(mast@meta.data[160:164]))-min(rowMeans(mast@meta.data[160:164])))

mast$positive.regulation.of.mast.cell.degranulation<-(rowMeans(mast@meta.data[165:173])-min(rowMeans(mast@meta.data[165:173])))/(max(rowMeans(mast@meta.data[165:173]))-min(rowMeans(mast@meta.data[165:173])))

mast$regulation.of.histamine.secretion.by.mast.cell<-(rowMeans(mast@meta.data[174:180])-min(rowMeans(mast@meta.data[174:180])))/(max(rowMeans(mast@meta.data[174:180]))-min(rowMeans(mast@meta.data[174:180])))

mast$regulation.of.mast.cell.activation<-(rowMeans(mast@meta.data[181:190])-min(rowMeans(mast@meta.data[181:190])))/(max(rowMeans(mast@meta.data[181:190]))-min(rowMeans(mast@meta.data[181:190])))

mast$regulation.of.mast.cell.activation.involved.in.immune.response<-(rowMeans(mast@meta.data[191:201])-min(rowMeans(mast@meta.data[191:201])))/(max(rowMeans(mast@meta.data[191:201]))-min(rowMeans(mast@meta.data[191:201])))

mast$regulation.of.mast.cell.chemotaxis<-(rowMeans(mast@meta.data[202:207])-min(rowMeans(mast@meta.data[202:207])))/(max(rowMeans(mast@meta.data[202:207]))-min(rowMeans(mast@meta.data[202:207])))

mast$regulation.of.mast.cell.degranulation<-(rowMeans(mast@meta.data[208:223])-min(rowMeans(mast@meta.data[208:223])))/(max(rowMeans(mast@meta.data[208:223]))-min(rowMeans(mast@meta.data[208:223])))

mast$Fc.epsilon.RI.signaling.pathway<-(rowMeans(mast@meta.data[224:290])-min(rowMeans(mast@meta.data[224:290])))/(max(rowMeans(mast@meta.data[224:290]))-min(rowMeans(mast@meta.data[224:290])))

mast$VEGF.signaling.pathway<-(rowMeans(mast@meta.data[291:349])-min(rowMeans(mast@meta.data[291:349])))/(max(rowMeans(mast@meta.data[291:349]))-min(rowMeans(mast@meta.data[291:349])))

mast$Antigen.processing.and.presentation<-(rowMeans(mast@meta.data[350:418])-min(rowMeans(mast@meta.data[350:418])))/(max(rowMeans(mast@meta.data[350:418]))-min(rowMeans(mast@meta.data[350:418])))

mast$Apoptosis<-(rowMeans(mast@meta.data[419:559])-min(rowMeans(mast@meta.data[419:559])))/(max(rowMeans(mast@meta.data[419:559]))-min(rowMeans(mast@meta.data[419:559])))

mast$Cellular.senescence<-(rowMeans(mast@meta.data[560:715])-min(rowMeans(mast@meta.data[560:715])))/(max(rowMeans(mast@meta.data[560:715]))-min(rowMeans(mast@meta.data[560:715])))

mast$Chemokine.signaling.pathway<-(rowMeans(mast@meta.data[716:905])-min(rowMeans(mast@meta.data[716:905])))/(max(rowMeans(mast@meta.data[716:905]))-min(rowMeans(mast@meta.data[716:905])))

mast$Complement.and.coagulation.cascades<-(rowMeans(mast@meta.data[906:990])-min(rowMeans(mast@meta.data[906:990])))/(max(rowMeans(mast@meta.data[906:990]))-min(rowMeans(mast@meta.data[906:990])))

mast$NOD.like.receptor.signaling.pathway<-(rowMeans(mast@meta.data[991:1168])-min(rowMeans(mast@meta.data[991:1168])))/(max(rowMeans(mast@meta.data[991:1168]))-min(rowMeans(mast@meta.data[991:1168])))

mast$Toll.like.receptor.signaling.pathway<-(rowMeans(mast@meta.data[1169:1270])-min(rowMeans(mast@meta.data[1169:1270])))/(max(rowMeans(mast@meta.data[1169:1270]))-min(rowMeans(mast@meta.data[1169:1270])))


table(Idents(mast))



mast_pathway_plot<-data.frame(Idents(mast),mast$mast.cell.degranulation, mast$mast.cell.mediated.immunity,
                              mast$positive.regulation.of.mast.cell.activation, mast$positive.regulation.of.mast.cell.activation.involved.in.immune.response,mast$positive.regulation.of.mast.cell.chemotaxis,
                              mast$positive.regulation.of.mast.cell.degranulation, mast$regulation.of.histamine.secretion.by.mast.cell,
                              mast$regulation.of.mast.cell.activation, mast$regulation.of.mast.cell.activation.involved.in.immune.response,
                              mast$regulation.of.mast.cell.chemotaxis,mast$regulation.of.mast.cell.degranulation,mast$Fc.epsilon.RI.signaling.pathway,
                              mast$VEGF.signaling.pathway, mast$Antigen.processing.and.presentation, mast$Apoptosis, mast$Cellular.senescence,
                              mast$Chemokine.signaling.pathway, mast$Complement.and.coagulation.cascades, mast$NOD.like.receptor.signaling.pathway,
                              mast$Toll.like.receptor.signaling.pathway)
colnames(mast_pathway_plot)<-c("Cluster","Pathway_1","Pathway_2","Pathway_3","inflammasome_pathway_score")

colnames(mast_pathway_plot)

mast_pathway_plot[1:4,]


mean_values <- mast_pathway_plot %>%
  group_by(Idents.mast.) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

mean_values <- data.frame(mean_values)

write.table(mean_values, file = "/Users/wangjun/Desktop/mast_pathway.csv", sep = ",", col.names = T, row.names = F)

for (i in 2:ncol(mean_values)){
  mean_values[,i] <- (mean_values[,i]-min(mean_values[,i]))/(max(mean_values[,i])-min(mean_values[,i]))
}

write.table(mean_values, file = "/Users/wangjun/Desktop/mast_pathway_normalize.csv", sep = ",", col.names = T, row.names = F)
