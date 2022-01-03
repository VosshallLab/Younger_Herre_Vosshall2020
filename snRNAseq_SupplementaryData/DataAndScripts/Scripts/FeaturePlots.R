suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(sctransform)
  library(ComplexHeatmap)
  library(future)
  library(scales)
  library(data.table)
  library(ggplot2)
})

theme_set(theme_cowplot())

#projectFolder="C:/ "
setwd(projectFolder)

#Load one of the following datasets
load("AntennaNeurons.RData")
#load("AntennaAllCells.RData")
#load("MaxPalpNeurons.RData") 
#load("MaxPalpAllCells.RData")

#Name dataset of interest "dataset"
dataset <- AntennaNeurons
#dataset <- AntennaAllCells
#dataset <- MaxPalpNeurons
#dataset <- MaxPalpAllCells

############################# Pseudo bulk (may be required for code further on)

counts <- dataset@assays$RNA@counts                  #Pull raw counts from Seurat object
pseudo_bulk <- rowSums(counts)                                             #Sum values across cells

pseudo_bulk <- as.data.frame(pseudo_bulk)                                  #Convert to dataframe
pseudo_bulk <- tibble::rownames_to_column(pseudo_bulk, "Gene")             #Pull out rownames as column

pseudo_bulk <- pseudo_bulk[!grepl("LOC*", pseudo_bulk$Gene),]              #get rid of LOC & MC genes (for just chemoreceptor analysis)                     #
pseudo_bulk <- pseudo_bulk[!grepl("MT*", pseudo_bulk$Gene),]                    
pseudo_bulk <- pseudo_bulk[order(pseudo_bulk[,2],decreasing=TRUE),]        #Sort genes in decreasing order 

pseudo_bulk.ORs <- pseudo_bulk[grepl("Or*", pseudo_bulk$Gene),]            #Separate ORs, IRs, GRs into three variables
pseudo_bulk.IRs <- pseudo_bulk[grepl("Ir*", pseudo_bulk$Gene),]
pseudo_bulk.GRs <- pseudo_bulk[grepl("Gr*", pseudo_bulk$Gene),]

rm(counts)

################################################################### Feature Plots (Figure S9H, Figure S13F)

geneName1 <- toString("LOC110678282") #repo, glial marker
geneName2 <- toString("LOC5564305") #grh, epithelial marker

#Neural markers
geneName3 <- toString("LOC5565901") #syt1
geneName4 <- toString("LOC5575210") #nompC, mechanosensory channel
geneName5 <- toString("LOC5564848") #CadN
geneName6 <- toString("LOC5570381") #brp
geneName7 <- toString("LOC5570204") #elav


FeaturePlot(dataset,  reduction = 'tsne', features = c(geneName1, geneName2, geneName3, geneName4) , pt.size = .1, slot = "data", order=TRUE, ncol = 2)
genes <- paste(toString(geneName1),toString(geneName2),toString(geneName3),toString(geneName4), sep="_")
Vln.filename <- paste("FeaturePlot_",toString(genes),".pdf", sep="")
ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 2000, height = 2000)


###################### Feature Plots (Figure S9H, Figure S14H, )

nrow(pseudo_bulk.ORs) 
nrow(pseudo_bulk.IRs) 
nrow(pseudo_bulk.GRs) 

#Feature Plots for ORs
for (i in 1:20) {
  #i <- 1
  c <- (i*6)-5
  geneName1 <- toString(pseudo_bulk.ORs$Gene[c])
  geneName2 <- toString(pseudo_bulk.ORs$Gene[c+1])
  geneName3 <- toString(pseudo_bulk.ORs$Gene[c+2])
  geneName4 <- toString(pseudo_bulk.ORs$Gene[c+3])
  geneName5 <- toString(pseudo_bulk.ORs$Gene[c+4])
  geneName6 <- toString(pseudo_bulk.ORs$Gene[c+5])
  FeaturePlot(dataset,  reduction = 'tsne', features = c(geneName1, geneName2, geneName3, geneName4, geneName5, geneName6) , pt.size = .1, slot = "data", order=TRUE, ncol = 3)
  genes <- paste(toString(geneName1),toString(geneName2),toString(geneName3),toString(geneName4),toString(geneName5),toString(geneName6), sep="_")
  Vln.filename <- paste("FeaturePlot_",toString(genes),".pdf", sep="")
  ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 3000, height = 2000) #width 2000, height 1000
}
dev.off()

#Feature Plots for IRs
for (i in 1:20) {
  #i <- 1
  c <- (i*6)-5
  geneName1 <- toString(pseudo_bulk.IRs$Gene[c])
  geneName2 <- toString(pseudo_bulk.IRs$Gene[c+1])
  geneName3 <- toString(pseudo_bulk.IRs$Gene[c+2])
  geneName4 <- toString(pseudo_bulk.IRs$Gene[c+3])
  geneName5 <- toString(pseudo_bulk.IRs$Gene[c+4])
  geneName6 <- toString(pseudo_bulk.IRs$Gene[c+5])
  FeaturePlot(dataset,  reduction = 'tsne', features = c(geneName1, geneName2, geneName3, geneName4, geneName5, geneName6) , pt.size = .1, slot = "data", order=TRUE, ncol = 3)
  genes <- paste(toString(geneName1),toString(geneName2),toString(geneName3),toString(geneName4),toString(geneName5),toString(geneName6), sep="_")
  Vln.filename <- paste("FeaturePlot_",toString(genes),".pdf", sep="")
  ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 3000, height = 2000) #width 2000, height 1000
}
dev.off()

#Feature Plots for GRs
for (i in 1:20) {
  i <- 1
  c <- (i*6)-5
  geneName1 <- toString(pseudo_bulk.GRs$Gene[c])
  geneName2 <- toString(pseudo_bulk.GRs$Gene[c+1])
  geneName3 <- toString(pseudo_bulk.GRs$Gene[c+2])
  geneName4 <- toString(pseudo_bulk.GRs$Gene[c+3])
  geneName5 <- toString(pseudo_bulk.GRs$Gene[c+4])
  geneName6 <- toString(pseudo_bulk.GRs$Gene[c+5])
  FeaturePlot(dataset,  reduction = 'tsne', features = c(geneName1, geneName2, geneName3, geneName4, geneName5, geneName6) , pt.size = .1, slot = "data", order=TRUE, ncol = 3)
  genes <- paste(toString(geneName1),toString(geneName2),toString(geneName3),toString(geneName4),toString(geneName5),toString(geneName6), sep="_")
  Vln.filename <- paste("FeaturePlot_",toString(genes),".pdf", sep="")
  ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 3000, height = 2000) #width 2000, height 1000
}
dev.off()

################################################################### Neural Marker Dot Plot (Figure S9I, Figure S13H)

DotPlot(dataset, features = c(geneName5, geneName6, geneName3, geneName7))+ scale_y_discrete(limits = rev)
ggsave("neuron_markers.pdf", limitsize = FALSE, units = "px", width = 2200, height = 4000)


################################################################### Neuron Cluster Marker Dot Plot (Figure S14E)

receptorL <- c("Orco", "Ir25a", "Ir76b", "Gr2", "Gr3", "Or8", "Or49", "LOC5575210")
DotPlot(dataset, features = receptorL, dot.scale = 25) + scale_y_discrete(limits = rev) +
  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) + scale_y_discrete(limits = rev)
Vln.filename <- paste0("DotPlot_MaxPalp.pdf")
ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 4000, height = 2500)



