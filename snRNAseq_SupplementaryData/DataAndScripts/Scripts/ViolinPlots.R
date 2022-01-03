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

################################################################### QC Violin Plots (Figure S13C-E, Figure S14A-C)

#Combined Data for QC plots
dat_Vln <- dataset

#Violin plot separated by cluster
VlnPlot(dat_Vln, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Idents(object = dat_Vln) <- "all_cells"         #Or remove cluster assignments to visualize all cells on one violin plot
#Idents(object = dat_Vln) <- "neurons"

VlnPlot(dat_Vln, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Save plot
Vln.filename <- paste0("SupplementalViolinPlot_QC.pdf")
ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 4000, height = 2500)


rm(dat_Vln)

################################################################### Gene Expression Violin Plots (Figure S10D, S11A)
dat_Vln <- dataset

#Idents(object = dat_Vln) <- "all_cells"
Idents(object = dat_Vln) <- "neurons"

#Orco
geneName <- toString("Orco")
expression.cutoff <- 2
VlnPlot(dat_Vln, features = geneName , pt.size = .1, slot = "data", same.y.lims=TRUE, y.max=6.1) + NoLegend() + geom_hline(yintercept=expression.cutoff, linetype="dashed", color = "blue") & theme(axis.title.x = element_blank())
Vln.filename <- paste0("SupplementalViolinPlot_Orco_exp2.pdf")
ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 1000, height = 2000)

#Ir25a
geneName <- toString("Ir25a")
expression.cutoff <- 1
VlnPlot(dat_Vln, features = geneName , pt.size = .1, slot = "data", same.y.lims=TRUE, y.max=6.1) + NoLegend() + geom_hline(yintercept=expression.cutoff, linetype="dashed", color = "blue") & theme(axis.title.x = element_blank())
Vln.filename <- paste0("SupplementalViolinPlot_",toString(geneName),"_exp",toString(expression.cutoff),".pdf")
ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 1000, height = 2000) 

#Ir41k
geneName <- toString("Ir41k")
expression.cutoff <- 0.5
VlnPlot(dat_Vln, features = geneName , pt.size = .1, slot = "data", same.y.lims=TRUE, y.max=6.1) + NoLegend() + geom_hline(yintercept=expression.cutoff, linetype="dashed", color = "blue") & theme(axis.title.x = element_blank())
Vln.filename <- paste0("SupplementalViolinPlot_",toString(geneName),"_exp",toString(expression.cutoff),".pdf")
ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 1000, height = 2000)

#Or4
geneName <- toString("Or4")
expression.cutoff <- 0.5
VlnPlot(dat_Vln, features = geneName , pt.size = .1, slot = "data", same.y.lims=TRUE, y.max=6.1) + NoLegend() + geom_hline(yintercept=expression.cutoff, linetype="dashed", color = "blue") & theme(axis.title.x = element_blank())
Vln.filename <- paste0("SupplementalViolinPlot_",toString(geneName),"_exp",toString(expression.cutoff),".pdf")
ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 1000, height = 2000) 

#Or82
geneName <- toString("Or82")
expression.cutoff <- 0.5
VlnPlot(dat_Vln, features = geneName , pt.size = .1, slot = "data", same.y.lims=TRUE, y.max=6.1) + NoLegend() + geom_hline(yintercept=expression.cutoff, linetype="dashed", color = "blue") & theme(axis.title.x = element_blank())
Vln.filename <- paste0("SupplementalViolinPlot_",toString(geneName),"_exp",toString(expression.cutoff),".pdf")
ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 1000, height = 2000)

#Ir64a
geneName <- toString("Ir64a")
expression.cutoff <- 0.5
VlnPlot(dat_Vln, features = geneName , pt.size = .1, slot = "data", same.y.lims=TRUE, y.max=6.1) + NoLegend() + geom_hline(yintercept=expression.cutoff, linetype="dashed", color = "blue") & theme(axis.title.x = element_blank())
Vln.filename <- paste0("SupplementalViolinPlot_",toString(geneName),"_exp",toString(expression.cutoff),".pdf")
ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 1000, height = 2000)

#Ir41a
geneName <- toString("Ir41a")
expression.cutoff <- 0.5
VlnPlot(dat_Vln, features = geneName , pt.size = .1, slot = "data", same.y.lims=TRUE, y.max=6.1) + NoLegend() + geom_hline(yintercept=expression.cutoff, linetype="dashed", color = "blue") & theme(axis.title.x = element_blank())
Vln.filename <- paste0("SupplementalViolinPlot_",toString(geneName),"_exp",toString(expression.cutoff),".pdf")
ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 1000, height = 2000)

#Ir8a
geneName <- toString("Ir8a")
expression.cutoff <- 0.5
VlnPlot(dat_Vln, features = geneName , pt.size = .1, slot = "data", same.y.lims=TRUE, y.max=6.1) + NoLegend() + geom_hline(yintercept=expression.cutoff, linetype="dashed", color = "blue") & theme(axis.title.x = element_blank())
Vln.filename <- paste0("SupplementalViolinPlot_",toString(geneName),"_exp",toString(expression.cutoff),".pdf")
ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 1000, height = 2000)

#Gr77
geneName <- toString("Gr77")
expression.cutoff <- 0.5
VlnPlot(dat_Vln, features = geneName , pt.size = .1, slot = "data", same.y.lims=TRUE, y.max=6.1) + NoLegend() + geom_hline(yintercept=expression.cutoff, linetype="dashed", color = "blue") & theme(axis.title.x = element_blank())
Vln.filename <- paste0("SupplementalViolinPlot_",toString(geneName),"_exp",toString(expression.cutoff),".pdf")
ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 1000, height = 2000) 

#Or47
geneName <- toString("Or47")
expression.cutoff <- 0.5
VlnPlot(dat_Vln, features = geneName , pt.size = .1, slot = "data", same.y.lims=TRUE, y.max=6.1) + NoLegend() + geom_hline(yintercept=expression.cutoff, linetype="dashed", color = "blue") & theme(axis.title.x = element_blank())
Vln.filename <- paste0("SupplementalViolinPlot_",toString(geneName),"_exp",toString(expression.cutoff),".pdf")
ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 1000, height = 2000)


rm(dat_Vln)

################################################################### Violin Plot on all genes (Supplementary Data)
dat_Vln <- dataset
Idents(object = dat_Vln) <- "neurons"

nrow(pseudo_bulk.ORs) 
nrow(pseudo_bulk.IRs) 
nrow(pseudo_bulk.GRs) 

#Violin Plots for ORs
for (i in 1:20) {
  #i <- 1
  c <- (i*6)-5
  geneName1 <- toString(pseudo_bulk.ORs$Gene[c])
  geneName2 <- toString(pseudo_bulk.ORs$Gene[c+1])
  geneName3 <- toString(pseudo_bulk.ORs$Gene[c+2])
  geneName4 <- toString(pseudo_bulk.ORs$Gene[c+3])
  geneName5 <- toString(pseudo_bulk.ORs$Gene[c+4])
  geneName6 <- toString(pseudo_bulk.ORs$Gene[c+5])
  VlnPlot(dat_Vln, features = c(geneName1, geneName2, geneName3, geneName4, geneName5, geneName6) , pt.size = .1, slot = "data", same.y.lims=TRUE, y.max=6.1) + NoLegend()  & theme(axis.title.x = element_blank())
  genes <- paste(toString(geneName1),toString(geneName2),toString(geneName3),toString(geneName4),toString(geneName5),toString(geneName6), sep="_")
  Vln.filename <- paste("ViolinPlot_",toString(genes),".pdf", sep="")
  ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 2000, height = 3000)
}


#Violin Plots for IRs
for (i in 1:20) {
  #i <- 1
  c <- (i*6)-5
  geneName1 <- toString(pseudo_bulk.IRs$Gene[c])
  geneName2 <- toString(pseudo_bulk.IRs$Gene[c+1])
  geneName3 <- toString(pseudo_bulk.IRs$Gene[c+2])
  geneName4 <- toString(pseudo_bulk.IRs$Gene[c+3])
  geneName5 <- toString(pseudo_bulk.IRs$Gene[c+4])
  geneName6 <- toString(pseudo_bulk.IRs$Gene[c+5])
  VlnPlot(dat_Vln, features = c(geneName1, geneName2, geneName3, geneName4, geneName5, geneName6) , pt.size = .1, slot = "data", same.y.lims=TRUE, y.max=6.1) + NoLegend()  & theme(axis.title.x = element_blank())
  genes <- paste(toString(geneName1),toString(geneName2),toString(geneName3),toString(geneName4),toString(geneName5),toString(geneName6), sep="_")
  Vln.filename <- paste("ViolinPlot_",toString(genes),".pdf", sep="")
  ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 2000, height = 3000)
}

#Violin Plots for GRs
for (i in 1:20) {
  #i <- 1
  c <- (i*6)-5
  geneName1 <- toString(pseudo_bulk.GRs$Gene[c])
  geneName2 <- toString(pseudo_bulk.GRs$Gene[c+1])
  geneName3 <- toString(pseudo_bulk.GRs$Gene[c+2])
  geneName4 <- toString(pseudo_bulk.GRs$Gene[c+3])
  geneName5 <- toString(pseudo_bulk.GRs$Gene[c+4])
  geneName6 <- toString(pseudo_bulk.GRs$Gene[c+5])
  VlnPlot(dat_Vln, features = c(geneName1, geneName2, geneName3, geneName4, geneName5, geneName6) , pt.size = .1, slot = "data", same.y.lims=TRUE, y.max=6.1) + NoLegend()  & theme(axis.title.x = element_blank())
  genes <- paste(toString(geneName1),toString(geneName2),toString(geneName3),toString(geneName4),toString(geneName5),toString(geneName6), sep="_")
  Vln.filename <- paste("ViolinPlot_",toString(genes),".pdf", sep="")
  ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 2000, height = 3000)
}

rm(dat_Vln)


################################################################### Violin Plot by Cluster (Figure S14D, Supplementary Data)

nrow(pseudo_bulk.ORs) #36
nrow(pseudo_bulk.IRs) #40
nrow(pseudo_bulk.GRs) #40

#Violin Plots for ORs
for (i in 1:36) {
  #i <- 1
  c <- (i*6)-5
  geneName1 <- toString(pseudo_bulk.ORs$Gene[c])
  geneName2 <- toString(pseudo_bulk.ORs$Gene[c+1])
  geneName3 <- toString(pseudo_bulk.ORs$Gene[c+2])
  geneName4 <- toString(pseudo_bulk.ORs$Gene[c+3])
  geneName5 <- toString(pseudo_bulk.ORs$Gene[c+4])
  geneName6 <- toString(pseudo_bulk.ORs$Gene[c+5])
  VlnPlot(dataset, features = c(geneName1, geneName2, geneName3, geneName4, geneName5, geneName6) , pt.size = .1, slot = "data", same.y.lims=TRUE, y.max=6.1) + NoLegend()  & theme(axis.title.x = element_blank())
  genes <- paste(toString(geneName1),toString(geneName2),toString(geneName3),toString(geneName4),toString(geneName5),toString(geneName6), sep="_")
  Vln.filename <- paste("ViolinPlot_clusters_",toString(genes),".pdf", sep="")
  ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 2000, height = 3000)0
}


#Violin Plots for IRs
for (i in 1:20) {
  #i <- 1
  c <- (i*6)-5
  geneName1 <- toString(pseudo_bulk.IRs$Gene[c])
  geneName2 <- toString(pseudo_bulk.IRs$Gene[c+1])
  geneName3 <- toString(pseudo_bulk.IRs$Gene[c+2])
  geneName4 <- toString(pseudo_bulk.IRs$Gene[c+3])
  geneName5 <- toString(pseudo_bulk.IRs$Gene[c+4])
  geneName6 <- toString(pseudo_bulk.IRs$Gene[c+5])
  VlnPlot(dataset, features = c(geneName1, geneName2, geneName3, geneName4, geneName5, geneName6) , pt.size = .1, slot = "data", same.y.lims=TRUE, y.max=6.1) + NoLegend()  & theme(axis.title.x = element_blank())
  genes <- paste(toString(geneName1),toString(geneName2),toString(geneName3),toString(geneName4),toString(geneName5),toString(geneName6), sep="_")
  Vln.filename <- paste("ViolinPlot_clusters_",toString(genes),".pdf", sep="")
  ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 2000, height = 3000)
}

#Violin Plots for GRs
for (i in 1:20) {
  #i <- 1
  c <- (i*6)-5
  geneName1 <- toString(pseudo_bulk.GRs$Gene[c])
  geneName2 <- toString(pseudo_bulk.GRs$Gene[c+1])
  geneName3 <- toString(pseudo_bulk.GRs$Gene[c+2])
  geneName4 <- toString(pseudo_bulk.GRs$Gene[c+3])
  geneName5 <- toString(pseudo_bulk.GRs$Gene[c+4])
  geneName6 <- toString(pseudo_bulk.GRs$Gene[c+5])
  VlnPlot(dataset, features = c(geneName1, geneName2, geneName3, geneName4, geneName5, geneName6) , pt.size = .1, slot = "data", same.y.lims=TRUE, y.max=6.1) + NoLegend()  & theme(axis.title.x = element_blank())
  genes <- paste(toString(geneName1),toString(geneName2),toString(geneName3),toString(geneName4),toString(geneName5),toString(geneName6), sep="_")
  Vln.filename <- paste("ViolinPlot_clusters_",toString(genes),".pdf", sep="")
  ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 2000, height = 3000)
}
