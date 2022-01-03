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
load("AntennaNeurons.RData")

############################################################################################################### 
############################################################################################################### 
############################################################################################################### 
############################################################################################################### Figure 5F: Heatmap Ir64a

geneName <- toString("Ir64a")

#..............................................................................................Violin Plot
#VlnPlot(dat_Vln, features = geneName, pt.size = .1) + scale_y_continuous(breaks=seq(0,7,0.5)) + NoLegend() 
#Vln.filename <- paste("ViolinPlot_",toString(geneName),".pdf", sep="")
#ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 2000, height = 3200)

#..............................................................................................Large Heatmap & Clustering
expression.cutoff <- 0.5
expr <- FetchData(object = AntennaNeurons, vars = geneName)
dat <- AntennaNeurons[, which(x = expr > expression.cutoff)]
dim(dat@meta.data)

dat <- RunPCA(dat, npcs = 50, verbose = F)
dat <- FindNeighbors(dat, dims = 1:50, verbose = F)
dat <- FindClusters(dat, resolution = c(3))

ae <- AverageExpression(dat)                        #get avg expression of all genes in cluster       
ae <- ae$RNA                                        #pull out only relevant data                               
ae <- ae[order(ae[,1],decreasing=TRUE),]            #sort from largest to smallest                          
ae <- as.data.frame(ae)                             #convert to df      
ae <- tibble::rownames_to_column(ae, "Gene")        #make gene names a column                           
ae <- ae[!grepl("LOC*", ae$Gene),]                  #get rid of LOC & MC genes                     
ae <- ae[!grepl("MT*", ae$Gene),]
ae.ORs <- ae[grepl("Or*", ae$Gene),]
ae.IRs <- ae[grepl("Ir*", ae$Gene),]
ae.GRs <- ae[grepl("Gr*", ae$Gene),]
#ae <- rbind(ae.ORs, ae.IRs, ae.GRs)
ae <- rbind(ae.IRs, ae.ORs, ae.GRs)
#ae <- rbind(ae.GRs, ae.ORs, ae.IRs)
x <- ae[grepl("Orco", ae$Gene),]                    #Pull out co-receptors               
x[2,] <- ae[grepl("Ir25a", ae$Gene),]
x[3,] <- ae[grepl("Ir76b", ae$Gene),]
x[4,] <- ae[grepl("Ir8a", ae$Gene),]
ae <- ae[!grepl("Orco", ae$Gene),]                  #Delete co-receptors from main variable                     
ae <- ae[!grepl("Ir25a", ae$Gene),]
ae <- ae[!grepl("Ir76b", ae$Gene),]
ae <- ae[!grepl("Ir8a", ae$Gene),] 
ae <- rbind(x, ae)                                  #Add co-receptors to top     

DoHeatmap(object = dat, features = ae[,1], label=TRUE, group.bar.height =0.001, draw.lines=TRUE, size=4, angle=90) +
  theme(axis.text.y = element_text(size = 4)) #+ NoLegend()

#DoHeatmap(object = dat, features = ae[,1], label=TRUE, group.bar.height =0.001, draw.lines=TRUE, size=4, angle=90, slot="data") +            
#  theme(axis.text.y = element_text(size = 4)) #+ NoLegend()


#heatmap.filename <- paste("coexpression_heatmap_",toString(geneName),"_exp",toString(expression.cutoff),"_withLines.png", sep="")
#ggsave(heatmap.filename, limitsize = FALSE, units = "px", width = 4000, height = 2000)


#..............................................................................................Simplified heatmap

#Selected clusters of interest as identified in large coexpression heatmap
clusters.OI <- c(3,5,7,8,9,11,12)

#Selected genes of interest as identified in large coexpression heatmap
heatmapFeatures <- c('Ir64a',
                     'Or11',
                     'Or36',
                     'Or84',
                     'Or87',
                     'Or88',
                     'Or94',
                     'Or97',
                     'Or100',
                     'Or103',
                     'Or104',
                     'Or105',
                     'Or111',
                     'Or112',
                     'Or113')

#Subset cluster from filtered population of cells above threshold
cluster.subset <- subset(dat,idents = c(clusters.OI))

#Reorganizing clusters
my_levels <- c(3, 9, 12, 5, 11, 8, 7)
cluster.subset@active.ident <- factor(x = cluster.subset@active.ident, levels = my_levels)

DoHeatmap(object = cluster.subset, features = heatmapFeatures, label=TRUE, group.bar.height =0.01, draw.lines=TRUE, lines.width=4, size=4) +            #Make dotplot on sorted gene list
  theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colors = c("blue","black", "yellow"))

heatmap.filename <- paste("subset_heatmap_",toString(geneName),"_BlueBlackYellow_v3_ORs,ClustersReordered.png", sep="")
ggsave(heatmap.filename, limitsize = FALSE, width=9.2, height=6.3)
heatmap.filename <- paste("subset_heatmap_",toString(geneName),"_BlueBlackYellow_v3_ORs,ClustersReordered.pdf", sep="")
ggsave(heatmap.filename, limitsize = FALSE, width=9.2, height=6.3)


############################################################################################################### 
############################################################################################################### 
############################################################################################################### 
############################################################################################################### Figure 5G: Heatmap Or4

geneName <- toString("Or4")

#..............................................................................................Violin Plot
#VlnPlot(dat_Vln, features = geneName, pt.size = .1) + scale_y_continuous(breaks=seq(0,7,0.5)) + NoLegend() 
#Vln.filename <- paste("ViolinPlot_",toString(geneName),".pdf", sep="")
#ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 2000, height = 3200)


#..............................................................................................Large Heatmap & Clustering
expression.cutoff <- 0.5
expr <- FetchData(object = AntennaNeurons, vars = geneName)
dat <- AntennaNeurons[, which(x = expr > expression.cutoff)]
dim(dat@meta.data)

dat <- RunPCA(dat, npcs = 50, verbose = F)
dat <- FindNeighbors(dat, dims = 1:50, verbose = F)
dat <- FindClusters(dat, resolution = c(3))

ae <- AverageExpression(dat)                  
ae <- ae$RNA                                  
ae <- ae[order(ae[,1],decreasing=TRUE),]      
ae <- as.data.frame(ae)                       
ae <- tibble::rownames_to_column(ae, "Gene")  
ae <- ae[!grepl("LOC*", ae$Gene),]            
ae <- ae[!grepl("MT*", ae$Gene),]
ae.ORs <- ae[grepl("Or*", ae$Gene),]
ae.IRs <- ae[grepl("Ir*", ae$Gene),]
ae.GRs <- ae[grepl("Gr*", ae$Gene),]
ae <- rbind(ae.ORs, ae.IRs, ae.GRs)
#ae <- rbind(ae.IRs, ae.ORs, ae.GRs)
#ae <- rbind(ae.GRs, ae.ORs, ae.IRs)
x <- ae[grepl("Orco", ae$Gene),]              
x[2,] <- ae[grepl("Ir25a", ae$Gene),]
x[3,] <- ae[grepl("Ir76b", ae$Gene),]
x[4,] <- ae[grepl("Ir8a", ae$Gene),]
ae <- ae[!grepl("Orco", ae$Gene),]            
ae <- ae[!grepl("Ir25a", ae$Gene),]
ae <- ae[!grepl("Ir76b", ae$Gene),]
ae <- ae[!grepl("Ir8a", ae$Gene),] 
ae <- rbind(x, ae)                            

DoHeatmap(object = dat, features = ae[,1], label=TRUE, group.bar.height =0.001, draw.lines=TRUE, size=4, angle=90) + 
  theme(axis.text.y = element_text(size = 4)) #+ NoLegend()               #Make genes smaller font size

#DoHeatmap(object = cluster.subset, features = heatmapFeatures, label=TRUE, group.bar.height =0.01, draw.lines=TRUE, lines.width=4, size=4, slot = "data") + 
#  theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colors = c("blue","black", "yellow"))


#heatmap.filename <- paste("coexpression_heatmap_",toString(geneName),"_exp",toString(expression.cutoff),".png", sep="")
#ggsave(heatmap.filename, limitsize = FALSE, units = "px", width = 4000, height = 2000)


#..............................................................................................Small heatmap
clusters.OI <- c(0,1,4,5)

heatmapFeatures <- c('Or4',
                     'Or47',
                     'Or71',
                     'Or82',
                     'Or84',
                     'Or87',
                     'Or88',
                     'Or97',
                     'Or104',
                     'Or105',
                     'Or111',
                     'Or112',
                     'Or113')

cluster.subset <- subset(dat,idents = c(clusters.OI))

DoHeatmap(object = cluster.subset, features = heatmapFeatures, label=TRUE, group.bar.height =0.01, draw.lines=TRUE, lines.width=4, size=4) + 
  theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colors = c("blue","black", "yellow"))

heatmap.filename <- paste("subset_heatmap_",toString(geneName),"_BlueBlackYellow_v2.png", sep="")
ggsave(heatmap.filename, limitsize = FALSE, width=9.2, height=6.3)
heatmap.filename <- paste("subset_heatmap_",toString(geneName),"_BlueBlackYellow_v2.pdf", sep="")
ggsave(heatmap.filename, limitsize = FALSE, width=9.2, height=6.3)

############################################################################################################### 
############################################################################################################### 
############################################################################################################### 
############################################################################################################### Figure 5H: Heatmap Or82

geneName <- toString("Or82")

#..............................................................................................Violin Plot
#VlnPlot(dat_Vln, features = geneName, pt.size = .1) + scale_y_continuous(breaks=seq(0,7,0.5)) + NoLegend() 
#Vln.filename <- paste("ViolinPlot_",toString(geneName),".pdf", sep="")
#ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 2000, height = 3200) 0

#..............................................................................................Large Heatmap & Clustering
expression.cutoff <- 0.5
expr <- FetchData(object = AntennaNeurons, vars = geneName)
dat <- AntennaNeurons[, which(x = expr > expression.cutoff)]
dim(dat@meta.data)

dat <- RunPCA(dat, npcs = 50, verbose = F)
dat <- FindNeighbors(dat, dims = 1:50, verbose = F)
dat <- FindClusters(dat, resolution = c(0.5, 0.8, 1, 1.5, 2, 3)[6])

ae <- AverageExpression(dat)                  
ae <- ae$RNA                                  
ae <- ae[order(ae[,1],decreasing=TRUE),]      
ae <- as.data.frame(ae)                       
ae <- tibble::rownames_to_column(ae, "Gene")  
ae <- ae[!grepl("LOC*", ae$Gene),]            
ae <- ae[!grepl("MT*", ae$Gene),]
ae.ORs <- ae[grepl("Or*", ae$Gene),]
ae.IRs <- ae[grepl("Ir*", ae$Gene),]
ae.GRs <- ae[grepl("Gr*", ae$Gene),]
ae <- rbind(ae.ORs, ae.IRs, ae.GRs)
#ae <- rbind(ae.IRs, ae.ORs, ae.GRs)
#ae <- rbind(ae.GRs, ae.ORs, ae.IRs)
x <- ae[grepl("Orco", ae$Gene),]              
x[2,] <- ae[grepl("Ir25a", ae$Gene),]
x[3,] <- ae[grepl("Ir76b", ae$Gene),]
x[4,] <- ae[grepl("Ir8a", ae$Gene),]
ae <- ae[!grepl("Orco", ae$Gene),]            
ae <- ae[!grepl("Ir25a", ae$Gene),]
ae <- ae[!grepl("Ir76b", ae$Gene),]
ae <- ae[!grepl("Ir8a", ae$Gene),] 
ae <- rbind(x, ae)                            

DoHeatmap(object = dat, features = ae[,1], label=TRUE, group.bar.height =0.001, draw.lines=TRUE, size=4, angle=90) +      
  theme(axis.text.y = element_text(size = 4)) #+ NoLegend()

#heatmap.filename <- paste("coexpression_heatmap_",toString(geneName),"_exp",toString(expression.cutoff),"_withLines.png", sep="")
#ggsave(heatmap.filename, limitsize = FALSE, units = "px", width = 4000, height = 2000) 

#..............................................................................................Small heatmap
clusters.OI <- c(2,3,4,6,7,8,9)
heatmapFeatures <- c('Or82',
                     'Or3',
                     'Or4',
                     'Or6',
                     'Or28',
                     'Or47',
                     'Or71',
                     'Or80',
                     'Or81',
                     'Or104',
                     'Or105',
                     'Or111',
                     'Or112',
                     'Or113',
                     'Ir41c',
                     'Ir41l',
                     'Ir41m')

cluster.subset <- subset(dat,idents = c(clusters.OI))

DoHeatmap(object = cluster.subset, features = heatmapFeatures, label=TRUE, group.bar.height =0.01, draw.lines=TRUE, lines.width=4, size=4) +   
  theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colors = c("blue","black", "yellow"))

heatmap.filename <- paste("subset_heatmap_",toString(geneName),"_BlueBlackYellow.png", sep="")
ggsave(heatmap.filename, limitsize = FALSE, width=9.2, height=6.3)
heatmap.filename <- paste("subset_heatmap_",toString(geneName),"_BlueBlackYellow.pdf", sep="")
ggsave(heatmap.filename, limitsize = FALSE, width=9.2, height=6.3)




############################################################################################################### 
############################################################################################################### 
############################################################################################################### 
############################################################################################################### Figure 5I: Heatmap Ir41k

geneName <- toString("Ir41k")

#..............................................................................................Violin Plot
#VlnPlot(dat_Vln, features = geneName, pt.size = .1) + scale_y_continuous(breaks=seq(0,7,0.5)) + NoLegend() 
#Vln.filename <- paste("ViolinPlot_",toString(geneName),".pdf", sep="")
#ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 2000, height = 3200)

#..............................................................................................Large Heatmap & Clustering
expression.cutoff <- 0.5
expr <- FetchData(object = AntennaNeurons, vars = geneName)
dat <- AntennaNeurons[, which(x = expr > expression.cutoff)]
dim(dat@meta.data)
dim(AntennaNeurons@meta.data)

dat <- RunPCA(dat, npcs = 50, verbose = F)
dat <- FindNeighbors(dat, dims = 1:50, verbose = F)
dat <- FindClusters(dat, resolution = 3)

ae <- AverageExpression(dat)                   
ae <- ae$RNA                                   
ae <- ae[order(ae[,1],decreasing=TRUE),]       
ae <- as.data.frame(ae)                        
ae <- tibble::rownames_to_column(ae, "Gene")   
ae <- ae[!grepl("LOC*", ae$Gene),]             
ae <- ae[!grepl("MT*", ae$Gene),]
ae.ORs <- ae[grepl("Or*", ae$Gene),]
ae.IRs <- ae[grepl("Ir*", ae$Gene),]
ae.GRs <- ae[grepl("Gr*", ae$Gene),]
#ae <- rbind(ae.ORs, ae.IRs, ae.GRs)
ae <- rbind(ae.IRs, ae.ORs, ae.GRs)
#ae <- rbind(ae.GRs, ae.ORs, ae.IRs)
x <- ae[grepl("Orco", ae$Gene),]               
x[2,] <- ae[grepl("Ir25a", ae$Gene),]
x[3,] <- ae[grepl("Ir76b", ae$Gene),]
x[4,] <- ae[grepl("Ir8a", ae$Gene),]
ae <- ae[!grepl("Orco", ae$Gene),]             
ae <- ae[!grepl("Ir25a", ae$Gene),]
ae <- ae[!grepl("Ir76b", ae$Gene),]
ae <- ae[!grepl("Ir8a", ae$Gene),] 
ae <- rbind(x, ae)                             

DoHeatmap(object = dat, features = ae[,1], label=TRUE, group.bar.height =0.001, draw.lines=TRUE, size=4, angle=90) +    
  theme(axis.text.y = element_text(size = 4)) #+ NoLegend()

#DoHeatmap(object = dat, features = ae[,1], label=TRUE, group.bar.height =0.001, draw.lines=TRUE, size=4, angle=90, slot = "data") +    
#  theme(axis.text.y = element_text(size = 4)) #+ NoLegend() 

heatmap.filename <- paste("coexpression_heatmap_",toString(geneName),"_exp",toString(expression.cutoff),"_withLines.png", sep="")
ggsave(heatmap.filename, limitsize = FALSE, units = "px", width = 4000, height = 2000) 


#..............................................................................................Small heatmap

clusters.OI <- c(2,3,4,6,7,8)

heatmapFeatures <- c('Ir41k',
                     'Ir41j',
                     'Or11',
                     'Or38',
                     'Or94',
                     'Or97',
                     'Or100',
                     'Or104',
                     'Or105',
                     'Or111',
                     'Or112',
                     'Or113')


cluster.subset <- subset(dat,idents = c(clusters.OI))

DoHeatmap(object = cluster.subset, features = heatmapFeatures, label=TRUE, group.bar.height =0.01, draw.lines=TRUE, lines.width=4, size=4) +  
  theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colors = c("blue","black", "yellow"))

#DoHeatmap(object = cluster.subset, features = heatmapFeatures, label=TRUE, group.bar.height =0.01, draw.lines=TRUE, lines.width=4, size=4, slot = "data") +  
#  theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colors = c("black", "yellow"))

#..............................................................................................Reclustered.
#For reclustering, used 9 PCA and 1.9 resolution

dat.r <- cluster.subset

dat.r <- RunPCA(dat.r, features = VariableFeatures(object = dat.r))
#VizDimLoadings(dat.r, dims = 1:9, reduction = "pca") &theme(axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold"))
#DimHeatmap(dat.r, dims = 1:20, nfeatures = 20, cells = 500, balanced = T)
#DimPlot(dat.r, reduction = "pca")

#ElbowPlot(dat.r)
dat.r <- FindNeighbors(dat.r, dims = 1:9)
dat.r <- FindClusters(dat.r, resolution = 1.9)
table(dat.r@meta.data$seurat_clusters)
#DimPlot(dat,label.size = 4,repel = T,label = T)

clusters.OI <- c(0,1,3,4,5,6)
cluster.subset <- subset(dat.r,idents = c(clusters.OI))

DoHeatmap(object = cluster.subset, features = heatmapFeatures, label=TRUE, group.bar.height =0.01, draw.lines=TRUE, lines.width=4) + 
  theme(axis.text.y = element_text(size = 12)) + scale_fill_gradientn(colors = c("blue","black", "yellow")) #+ NoLegend()

heatmap.filename <- paste("subset_heatmap_",toString(geneName),"_BlueBlackYellow_v2.png", sep="")
ggsave(heatmap.filename, limitsize = FALSE, width=9.2, height=6.3)
heatmap.filename <- paste("subset_heatmap_",toString(geneName),"_BlueBlackYellow_v2.pdf", sep="")
ggsave(heatmap.filename, limitsize = FALSE, width=9.2, height=6.3)



############################################################################################################### 
############################################################################################################### 
############################################################################################################### 
############################################################################################################### Figure S11C: Heatmap Ir41a

geneName <- toString("Ir41a")

#dat.unscale <- dat
#dat.scale <- ScaleData(dat.unscale)
#dat <- dat.scale

#..............................................................................................Violin Plot
#VlnPlot(dat_Vln, features = geneName, pt.size = .1) + scale_y_continuous(breaks=seq(0,7,0.5)) + NoLegend() 
#Vln.filename <- paste("ViolinPlot_",toString(geneName),".pdf", sep="")
#ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 2000, height = 3200) 

#..............................................................................................Large Heatmap & Clustering
expression.cutoff <- 0.5
expr <- FetchData(object = AntennaNeurons, vars = geneName)
dat <- AntennaNeurons[, which(x = expr > expression.cutoff)]
dim(dat@meta.data)

dat <- RunPCA(dat, npcs = 50, verbose = F)
dat <- FindNeighbors(dat, dims = 1:50, verbose = F)
dat <- FindClusters(dat, resolution = c(0.5, 0.8, 1, 1.5, 2, 3)[6])

ae <- AverageExpression(dat)                  
ae <- ae$RNA                                  
ae <- ae[order(ae[,1],decreasing=TRUE),]      
ae <- as.data.frame(ae)                       
ae <- tibble::rownames_to_column(ae, "Gene")  
ae <- ae[!grepl("LOC*", ae$Gene),]            
ae <- ae[!grepl("MT*", ae$Gene),]
ae.ORs <- ae[grepl("Or*", ae$Gene),]
ae.IRs <- ae[grepl("Ir*", ae$Gene),]
ae.GRs <- ae[grepl("Gr*", ae$Gene),]
#ae <- rbind(ae.ORs, ae.IRs, ae.GRs)
ae <- rbind(ae.IRs, ae.ORs, ae.GRs)
#ae <- rbind(ae.GRs, ae.ORs, ae.IRs)
x <- ae[grepl("Orco", ae$Gene),]              
x[2,] <- ae[grepl("Ir25a", ae$Gene),]
x[3,] <- ae[grepl("Ir76b", ae$Gene),]
x[4,] <- ae[grepl("Ir8a", ae$Gene),]
ae <- ae[!grepl("Orco", ae$Gene),]            
ae <- ae[!grepl("Ir25a", ae$Gene),
ae <- ae[!grepl("Ir76b", ae$Gene),
ae <- ae[!grepl("Ir8a", ae$Gene),]
ae <- rbind(x, ae)      

DoHeatmap(object = dat, features = ae[,1], label=TRUE, group.bar.height =0.001, draw.lines=TRUE, size=4, angle=90) +
  theme(axis.text.y = element_text(size = 4)) #+ NoLegend() 

#DoHeatmap(object = dat, features = ae[,1], label=TRUE, group.bar.height =0.001, draw.lines=TRUE, size=4, angle=90, slot = "data") +
#  theme(axis.text.y = element_text(size = 4)) #+ NoLegend()

#heatmap.filename <- paste("coexpression_heatmap_",toString(geneName),"_exp",toString(expression.cutoff),"_withLines.png", sep="")
#ggsave(heatmap.filename, limitsize = FALSE, units = "px", width = 4000, height = 2000)


#..............................................................................................Small heatmap
clusters.OI <- c(0,3,5,8,12,13,14,15,17)

heatmapFeatures <- c('Ir41a',
                     'Ir41b',
                     'Ir41c',
                     'Ir41j',
                     'Ir41k',
                     'Or36',
                     'Or52',
                     'Or72',
                     'Or84',
                     'Or87',
                     'Or88',
                     'Or94',
                     'Or97',
                     'Or100',
                     'Or103',
                     'Or104',
                     'Or105',
                     'Or111',
                     'Or112',
                     'Or113',
                     'Or115',
                     'Or119',
                     'Or121',
                     'Or122',
                     'Or125',
                     'Gr77')

cluster.subset <- subset(dat,idents = c(clusters.OI))

DoHeatmap(object = cluster.subset, features = heatmapFeatures, label=TRUE, group.bar.height =0.01, draw.lines=TRUE, lines.width=4, size=4) +  
  theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colors = c("blue","black", "yellow"))

#DoHeatmap(object = cluster.subset, features = heatmapFeatures, label=TRUE, group.bar.height =0.01, draw.lines=TRUE, lines.width=4, size=4, slot = "data") +    
#  theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colors = c("blue","black", "yellow"))


heatmap.filename <- paste("subset_heatmap_",toString(geneName),"_BlueBlackYellow_v2.png", sep="")
ggsave(heatmap.filename, limitsize = FALSE, width=9.2, height=6.3)
heatmap.filename <- paste("subset_heatmap_",toString(geneName),"_BlueBlackYellow_v2.pdf", sep="")
ggsave(heatmap.filename, limitsize = FALSE, width=9.2, height=6.3)


############################################################################################################### 
############################################################################################################### 
############################################################################################################### 
############################################################################################################### Figure S11E: Heatmap Ir8a

geneName <- toString("Ir8a")

#..............................................................................................Violin Plot
#VlnPlot(dat_Vln, features = geneName, pt.size = .1) + scale_y_continuous(breaks=seq(0,7,0.5)) + NoLegend() 
#Vln.filename <- paste("ViolinPlot_",toString(geneName),".pdf", sep="")
#ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 2000, height = 3200)

#..............................................................................................Large Heatmap & Clustering
expression.cutoff <- 0.5
expr <- FetchData(object = AntennaNeurons, vars = geneName)
dat <- AntennaNeurons[, which(x = expr > expression.cutoff)]
dim(dat@meta.data)

dat <- RunPCA(dat, npcs = 50, verbose = F)
dat <- FindNeighbors(dat, dims = 1:50, verbose = F)
dat <- FindClusters(dat, resolution = c(0.5, 0.8, 1, 1.5, 2, 3)[6])

ae <- AverageExpression(dat)                  
ae <- ae$RNA                                  
ae <- ae[order(ae[,1],decreasing=TRUE),]      
ae <- as.data.frame(ae)                       
ae <- tibble::rownames_to_column(ae, "Gene")  
ae <- ae[!grepl("LOC*", ae$Gene),]            
ae <- ae[!grepl("MT*", ae$Gene),]
ae.ORs <- ae[grepl("Or*", ae$Gene),]
ae.IRs <- ae[grepl("Ir*", ae$Gene),]
ae.GRs <- ae[grepl("Gr*", ae$Gene),]
#ae <- rbind(ae.ORs, ae.IRs, ae.GRs)
ae <- rbind(ae.IRs, ae.ORs, ae.GRs)
#ae <- rbind(ae.GRs, ae.ORs, ae.IRs)
x <- ae[grepl("Orco", ae$Gene),]              
x[2,] <- ae[grepl("Ir25a", ae$Gene),]
x[3,] <- ae[grepl("Ir76b", ae$Gene),]
x[4,] <- ae[grepl("Ir8a", ae$Gene),]
ae <- ae[!grepl("Orco", ae$Gene),]            
ae <- ae[!grepl("Ir25a", ae$Gene),]
ae <- ae[!grepl("Ir76b", ae$Gene),]
ae <- ae[!grepl("Ir8a", ae$Gene),] 
ae <- rbind(x, ae)                            

DoHeatmap(object = dat, features = ae[,1], label=TRUE, group.bar.height =0.001, draw.lines=TRUE, size=4, angle=90) + 
  theme(axis.text.y = element_text(size = 4)) #+ NoLegend()

#heatmap.filename <- paste("coexpression_heatmap_",toString(geneName),"_exp",toString(expression.cutoff),"_withLines.png", sep="")
#ggsave(heatmap.filename, limitsize = FALSE, units = "px", width = 4000, height = 2000)


#..............................................................................................Small heatmap
clusters.OI <- c(1,2,5,7,8,9,10,11,13,16)
heatmapFeatures <- c('Ir8a',
                     'Ir31a1',
                     'Ir31a2',
                     'Ir64a',
                     'Ir75a',
                     'Ir75b',
                     'Ir75c',
                     'Ir75d',
                     'Ir75e',
                     'Ir75f',
                     'Ir75g',
                     'Ir75h',
                     'Ir75i',
                     'Ir75k',
                     'Ir75l',
                     'Or2',
                     'Or104',
                     'Or105',
                     'Or111',
                     'Or112',
                     'Or113')

cluster.subset <- subset(dat,idents = c(clusters.OI))

DoHeatmap(object = cluster.subset, features = heatmapFeatures, label=TRUE, group.bar.height =0.01, draw.lines=TRUE, lines.width=4, size=4) +    
  theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colors = c("blue","black", "yellow"))

heatmap.filename <- paste("subset_heatmap_",toString(geneName),"_BlueBlackYellow.png", sep="")
ggsave(heatmap.filename, limitsize = FALSE, width=9.2, height=6.3)
heatmap.filename <- paste("subset_heatmap_",toString(geneName),"_BlueBlackYellow.pdf", sep="")
ggsave(heatmap.filename, limitsize = FALSE, width=9.2, height=6.3)



############################################################################################################### 
############################################################################################################### 
############################################################################################################### 
############################################################################################################### Figure S11E: Heatmap Gr77

geneName <- toString("Gr77")

#..............................................................................................Violin Plot
#VlnPlot(dat_Vln, features = geneName, pt.size = .1) + scale_y_continuous(breaks=seq(0,7,0.5)) + NoLegend() 
#Vln.filename <- paste("ViolinPlot_",toString(geneName),".pdf", sep="")
#ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 2000, height = 3200)

#..............................................................................................Large Heatmap & Clustering
expression.cutoff <- 0.5
expr <- FetchData(object = AntennaNeurons, vars = geneName)
dat <- AntennaNeurons[, which(x = expr > expression.cutoff)]
dim(dat@meta.data)

dat <- RunPCA(dat, npcs = 50, verbose = F)
dat <- FindNeighbors(dat, dims = 1:50, verbose = F)
dat <- FindClusters(dat, resolution = c(0.5, 0.8, 1, 1.5, 2, 3)[6])

ae <- AverageExpression(dat)                 
ae <- ae$RNA                                 
ae <- ae[order(ae[,1],decreasing=TRUE),]     
ae <- as.data.frame(ae)                      
ae <- tibble::rownames_to_column(ae, "Gene") 
ae <- ae[!grepl("LOC*", ae$Gene),]           
ae <- ae[!grepl("MT*", ae$Gene),]
ae.ORs <- ae[grepl("Or*", ae$Gene),]
ae.IRs <- ae[grepl("Ir*", ae$Gene),]
ae.GRs <- ae[grepl("Gr*", ae$Gene),]
#ae <- rbind(ae.ORs, ae.IRs, ae.GRs)
#ae <- rbind(ae.IRs, ae.ORs, ae.GRs)
ae <- rbind(ae.GRs, ae.ORs, ae.IRs)
x <- ae[grepl("Orco", ae$Gene),]             
x[2,] <- ae[grepl("Ir25a", ae$Gene),]
x[3,] <- ae[grepl("Ir76b", ae$Gene),]
x[4,] <- ae[grepl("Ir8a", ae$Gene),]
ae <- ae[!grepl("Orco", ae$Gene),]           
ae <- ae[!grepl("Ir25a", ae$Gene),]
ae <- ae[!grepl("Ir76b", ae$Gene),]
ae <- ae[!grepl("Ir8a", ae$Gene),] 
ae <- rbind(x, ae)                           

DoHeatmap(object = dat, features = ae[,1], label=TRUE, group.bar.height =0.001, draw.lines=TRUE, size=4, angle=90) + 
  theme(axis.text.y = element_text(size = 4)) #+ NoLegend()

#heatmap.filename <- paste("coexpression_heatmap_",toString(geneName),"_exp",toString(expression.cutoff),"_withLines.png", sep="")
#ggsave(heatmap.filename, limitsize = FALSE, units = "px", width = 4000, height = 2000


#..............................................................................................Small heatmap
clusters.OI <- c(0,2,3,5,7)
heatmapFeatures <- c('Gr77',
                     'Or84',
                     'Or87',
                     'Or94',
                     'Or97',
                     'Or100',
                     'Or104',
                     'Or111',
                     'Or113',
                     'Or112')

cluster.subset <- subset(dat,idents = c(clusters.OI))

DoHeatmap(object = cluster.subset, features = heatmapFeatures, label=TRUE, group.bar.height =0.01, draw.lines=TRUE, lines.width=4, size=4) +   
  theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colors = c("blue","black", "yellow"))

heatmap.filename <- paste("subset_heatmap_",toString(geneName),"_BlueBlackYellow.png", sep="")
ggsave(heatmap.filename, limitsize = FALSE, width=9.2, height=6.3)
heatmap.filename <- paste("subset_heatmap_",toString(geneName),"_BlueBlackYellow.pdf", sep="")
ggsave(heatmap.filename, limitsize = FALSE, width=9.2, height=6.3)

############################################################################################################### 
############################################################################################################### 
############################################################################################################### 
############################################################################################################### Figure S11F: Heatmap Or47

geneName <- toString("Or47")

#..............................................................................................Violin Plot
#VlnPlot(dat_Vln, features = geneName, pt.size = .1) + scale_y_continuous(breaks=seq(0,7,0.5)) + NoLegend() 
#Vln.filename <- paste("ViolinPlot_",toString(geneName),".pdf", sep="")
#ggsave(Vln.filename, limitsize = FALSE, units = "px", width = 2000, height = 3200)

#..............................................................................................Large Heatmap & Clustering
expression.cutoff <- 0.5
expr <- FetchData(object = AntennaNeurons, vars = geneName)
dat <- AntennaNeurons[, which(x = expr > expression.cutoff)]
dim(dat@meta.data)

dat <- RunPCA(dat, npcs = 50, verbose = F)
dat <- RunUMAP(dat, dims = 1:50, verbose = F)
dat <- RunTSNE(dat, dims = 1:50, verbose = F)
dat <- FindNeighbors(dat, dims = 1:50, verbose = F)
dat <- FindClusters(dat, resolution = c(0.5, 0.8, 1, 1.5, 2, 3)[6])

ae <- AverageExpression(dat)                 
ae <- ae$RNA                                 
ae <- ae[order(ae[,1],decreasing=TRUE),]     
ae <- as.data.frame(ae)                      
ae <- tibble::rownames_to_column(ae, "Gene") 
ae <- ae[!grepl("LOC*", ae$Gene),]           
ae <- ae[!grepl("MT*", ae$Gene),]
ae.ORs <- ae[grepl("Or*", ae$Gene),]
ae.IRs <- ae[grepl("Ir*", ae$Gene),]
ae.GRs <- ae[grepl("Gr*", ae$Gene),]
ae <- rbind(ae.ORs, ae.IRs, ae.GRs)
#ae <- rbind(ae.IRs, ae.ORs, ae.GRs)
#ae <- rbind(ae.GRs, ae.ORs, ae.IRs)
x <- ae[grepl("Orco", ae$Gene),]             
x[2,] <- ae[grepl("Ir25a", ae$Gene),]
x[3,] <- ae[grepl("Ir76b", ae$Gene),]
x[4,] <- ae[grepl("Ir8a", ae$Gene),]
ae <- ae[!grepl("Orco", ae$Gene),]           
ae <- ae[!grepl("Ir25a", ae$Gene),]
ae <- ae[!grepl("Ir76b", ae$Gene),]
ae <- ae[!grepl("Ir8a", ae$Gene),] 
ae <- rbind(x, ae)

DoHeatmap(object = dat, features = ae[,1], label=TRUE, group.bar.height =0.001, draw.lines=TRUE, size=4, angle=90) + 
  theme(axis.text.y = element_text(size = 4)) #+ NoLegend() 

#heatmap.filename <- paste("coexpression_heatmap_",toString(geneName),"_exp",toString(expression.cutoff),"_withLines.png", sep="")
#ggsave(heatmap.filename, limitsize = FALSE, units = "px", width = 4000, height = 2000) 


#..............................................................................................Small heatmap
clusters.OI <- c(0,1,2,3,5,7,8,12,14)
heatmapFeatures <- c('Or47',
                     'Or3',
                     'Or4',
                     'Or26',
                     'Or27',
                     'Or30',
                     'Or32',
                     'Or51',
                     'Or50',
                     'Or54',
                     'Or55',
                     'Or57',
                     'Or71',
                     'Or72',
                     'Or77',
                     'Or78',
                     'Or79',
                     'Or82',
                     'Or94',
                     'Or97',
                     'Or100',
                     'Or103',
                     'Or104',
                     'Or105',
                     'Or111',
                     'Or112',
                     'Or113',
                     'Or114',
                     'Or115',
                     'Or123',
                     'Gr77')

cluster.subset <- subset(dat,idents = c(clusters.OI))

DoHeatmap(object = cluster.subset, features = heatmapFeatures, label=TRUE, group.bar.height =0.01, draw.lines=TRUE, lines.width=4, size=4) 
  theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colors = c("blue","black", "yellow"))

heatmap.filename <- paste("subset_heatmap_",toString(geneName),"_BlueBlackYellow.png", sep="")
ggsave(heatmap.filename, limitsize = FALSE, width=9.2, height=6.3)
heatmap.filename <- paste("subset_heatmap_",toString(geneName),"_BlueBlackYellow.pdf", sep="")
ggsave(heatmap.filename, limitsize = FALSE, width=9.2, height=6.3)


