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

############################# Pseudo bulk (may be required for code further on)

counts <- AntennaNeurons@assays$RNA@counts                            #Pull raw counts from Seurat object
pseudo_bulk <- rowSums(counts)                                             #Sum values across cells

pseudo_bulk <- as.data.frame(pseudo_bulk)                                  #Convert to dataframe
pseudo_bulk <- tibble::rownames_to_column(pseudo_bulk, "Gene")             #Pull out rownames as column

pseudo_bulk <- pseudo_bulk[!grepl("LOC*", pseudo_bulk$Gene),]              #get rid of LOC & MC genes (for just chemoreceptor analysis)                     #
pseudo_bulk <- pseudo_bulk[!grepl("MT*", pseudo_bulk$Gene),]                    
pseudo_bulk <- pseudo_bulk[order(pseudo_bulk[,2],decreasing=TRUE),]        #Sort genes in decreasing order 

pseudo_bulk.ORs <- pseudo_bulk[grepl("Or*", pseudo_bulk$Gene),]            #Separate ORs, IRs, GRs into three variables
pseudo_bulk.IRs <- pseudo_bulk[grepl("Ir*", pseudo_bulk$Gene),]
pseudo_bulk.GRs <- pseudo_bulk[grepl("Gr*", pseudo_bulk$Gene),]


###################### Heatmap BY CLUSTER

table(AntennaNeurons@meta.data$seurat_clusters)

#Turned PCA, FindNeighbors, FindClusters starting at cluster 20 (i=21)

for (i in 1:55) {
  #i <- 1
  c <- i-1
  print(c)
  dat <- subset(AntennaNeurons,idents = c(c))
  dim(dat@meta.data)[1]
  
  #dat <- RunPCA(dat, npcs = 50, verbose = F)
  #dat <- FindNeighbors(dat, dims = 1:50, verbose = F)
  #dat <- FindClusters(dat, resolution = 3)
  
  ae <- AverageExpression(dat)                            #get avg expression of all genes in cluster       
  ae <- ae$RNA                                        #pull out only relevant data                               
  ae <- ae[order(ae[,1],decreasing=TRUE),]        #sort from largest to smallest                          
  ae <- as.data.frame(ae)                             #convert to df      
  ae <- tibble::rownames_to_column(ae, "Gene")        #make gene names a column                           
  ae <- ae[!grepl("LOC*", ae$Gene),]              #get rid of LOC & MC genes                     
  ae <- ae[!grepl("MT*", ae$Gene),]
  ae.ORs <- ae[grepl("Or*", ae$Gene),]
  ae.IRs <- ae[grepl("Ir*", ae$Gene),]
  ae.GRs <- ae[grepl("Gr*", ae$Gene),]
  ae <- rbind(ae.ORs, ae.IRs, ae.GRs)
  #ae <- rbind(ae.IRs, ae.ORs, ae.GRs)
  #ae <- rbind(ae.GRs, ae.ORs, ae.IRs)
  
  x <- ae[grepl("Orco", ae$Gene),]                    #Pull out co-receptors               
  x[2,] <- ae[grepl("Ir25a", ae$Gene),]
  x[3,] <- ae[grepl("Ir76b", ae$Gene),]
  x[4,] <- ae[grepl("Ir8a", ae$Gene),]
  
  ae <- ae[!grepl("Orco", ae$Gene),]              #Delete co-receptors from main variable                     
  ae <- ae[!grepl("Ir25a", ae$Gene),]
  ae <- ae[!grepl("Ir76b", ae$Gene),]
  ae <- ae[!grepl("Ir8a", ae$Gene),] 
  
  ae <- rbind(x, ae)                                  #Add co-receptors to top     
  
  
  heatmap.title = paste("Antenna Neuron Cluster ",toString(c),", ",toString(dim(dat@meta.data)[1])," cells",sep="")
  DoHeatmap(object = dat, features = ae[,1], label=TRUE, group.bar.height =0.001, draw.lines=FALSE, size=2, angle=90) +            
    theme(axis.text.y = element_text(size = 4), legend.text = element_text(size = 8)) + ggtitle(heatmap.title) #+ NoLegend()       
  
  heatmap.filename <- paste("coexpression_heatmap_cluster",toString(c),".pdf", sep="")
  ggsave(heatmap.filename, limitsize = FALSE, units = "px", width = 4000, height = 2000) #width 2000, height 1000
}

dev.off()





############################# Heatmaps, ORs

nrow(pseudo_bulk.ORs)

for (c in 1:nrow(pseudo_bulk.ORs)) {
  #c <- 50
  geneName <- toString(pseudo_bulk.ORs$Gene[c])
  #geneName <- toString(rna_assay_genes[c])
  print(paste(toString(c), geneName))
  
  expression.cutoff <- 0.5
  expr <- FetchData(object = AntennaNeurons, vars = geneName)
  dat <- AntennaNeurons[, which(x = expr > expression.cutoff)]
  print(dim(dat@meta.data)[1])
  
  if (dim(dat@meta.data)[1] > 100) {
    
    dat <- RunPCA(dat, npcs = 50, verbose = F)
    dat <- FindNeighbors(dat, dims = 1:50, verbose = F)
    dat <- FindClusters(dat, resolution = 3)
    #Idents(object = dat) <- "0"
    
    ae <- AverageExpression(dat)                            #get avg expression of all genes in cluster       
    ae <- ae$RNA                                        #pull out only relevant data                               
    ae <- ae[order(ae[,1],decreasing=TRUE),]        #sort from largest to smallest                          
    ae <- as.data.frame(ae)                             #convert to df      
    ae <- tibble::rownames_to_column(ae, "Gene")        #make gene names a column                           
    ae <- ae[!grepl("LOC*", ae$Gene),]              #get rid of LOC & MC genes                     
    ae <- ae[!grepl("MT*", ae$Gene),]
    ae.ORs <- ae[grepl("Or*", ae$Gene),]
    ae.IRs <- ae[grepl("Ir*", ae$Gene),]
    ae.GRs <- ae[grepl("Gr*", ae$Gene),]
    ae <- rbind(ae.ORs, ae.IRs, ae.GRs)
    #ae <- rbind(ae.IRs, ae.ORs, ae.GRs)
    #ae <- rbind(ae.GRs, ae.ORs, ae.IRs)
    
    x <- ae[grepl("Orco", ae$Gene),]                    #Pull out co-receptors               
    x[2,] <- ae[grepl("Ir25a", ae$Gene),]
    x[3,] <- ae[grepl("Ir76b", ae$Gene),]
    x[4,] <- ae[grepl("Ir8a", ae$Gene),]
    
    ae <- ae[!grepl("Orco", ae$Gene),]              #Delete co-receptors from main variable                     
    ae <- ae[!grepl("Ir25a", ae$Gene),]
    ae <- ae[!grepl("Ir76b", ae$Gene),]
    ae <- ae[!grepl("Ir8a", ae$Gene),] 
    
    ae <- rbind(x, ae)                                  #Add co-receptors to top     
    
    heatmap.title = paste(toString(geneName), " cells (Normalized Expression>",toString(expression.cutoff),", reclustered)",
                          "\n",toString(dim(dat@meta.data)[1])," cells",sep="")
    #heatmap.title = paste(toString(geneName), " cells (Normalized Expression>",toString(expression.cutoff),", original clusters)",
    #                      "\n",toString(dim(dat@meta.data)[1])," cells",sep="")
    #heatmap.title = paste("rna_",toString(geneName), " cells (Normalized Expression>",toString(expression.cutoff),", original clusters)",
    #                      "\n",toString(dim(dat@meta.data)[1])," cells",sep="")
    DoHeatmap(object = dat, features = ae[,1], label=TRUE, group.bar.height =0.001, draw.lines=FALSE, size=2, angle=90) +            #Make dotplot on sorted gene list
      theme(axis.text.y = element_text(size = 4), legend.text = element_text(size = 8)) + ggtitle(heatmap.title) #+ NoLegend()               #Make genes smaller font size
    
    heatmap.filename <- paste("coexpression_heatmap_",toString(geneName),"_exp",toString(expression.cutoff),"_reclustered.pdf", sep="")
    ggsave(heatmap.filename, limitsize = FALSE, units = "px", width = 4000, height = 2000)
  
  }
  else { 
    #dat <- RunPCA(dat, npcs = 50, verbose = F)
    dat <- FindNeighbors(dat, dims = 1:50, verbose = F)
    #dat <- FindClusters(dat, resolution = 3)
    #Idents(object = dat) <- "0"
    
    ae <- AverageExpression(dat)                            #get avg expression of all genes in cluster       
    ae <- ae$RNA                                        #pull out only relevant data                               
    ae <- ae[order(ae[,1],decreasing=TRUE),]        #sort from largest to smallest                          
    ae <- as.data.frame(ae)                             #convert to df      
    ae <- tibble::rownames_to_column(ae, "Gene")        #make gene names a column                           
    ae <- ae[!grepl("LOC*", ae$Gene),]              #get rid of LOC & MC genes                     
    ae <- ae[!grepl("MT*", ae$Gene),]
    ae.ORs <- ae[grepl("Or*", ae$Gene),]
    ae.IRs <- ae[grepl("Ir*", ae$Gene),]
    ae.GRs <- ae[grepl("Gr*", ae$Gene),]
    ae <- rbind(ae.ORs, ae.IRs, ae.GRs)
    #ae <- rbind(ae.IRs, ae.ORs, ae.GRs)
    #ae <- rbind(ae.GRs, ae.ORs, ae.IRs)
    
    x <- ae[grepl("Orco", ae$Gene),]                    #Pull out co-receptors               
    x[2,] <- ae[grepl("Ir25a", ae$Gene),]
    x[3,] <- ae[grepl("Ir76b", ae$Gene),]
    x[4,] <- ae[grepl("Ir8a", ae$Gene),]
    
    ae <- ae[!grepl("Orco", ae$Gene),]              #Delete co-receptors from main variable                     
    ae <- ae[!grepl("Ir25a", ae$Gene),]
    ae <- ae[!grepl("Ir76b", ae$Gene),]
    ae <- ae[!grepl("Ir8a", ae$Gene),] 
    
    ae <- rbind(x, ae)                                  #Add co-receptors to top     
    
    #heatmap.title = paste(toString(geneName), " cells (Normalized Expression>",toString(expression.cutoff),", reclustered)",
    #                      "\n",toString(dim(dat@meta.data)[1])," cells",sep="")
    heatmap.title = paste(toString(geneName), " cells (Normalized Expression>",toString(expression.cutoff),", original clusters)",
                          "\n",toString(dim(dat@meta.data)[1])," cells",sep="")
    #heatmap.title = paste("rna_",toString(geneName), " cells (Normalized Expression>",toString(expression.cutoff),", original clusters)",
    #                      "\n",toString(dim(dat@meta.data)[1])," cells",sep="")
    DoHeatmap(object = dat, features = ae[,1], label=TRUE, group.bar.height =0.001, draw.lines=FALSE, size=2, angle=90) +            #Make dotplot on sorted gene list
      theme(axis.text.y = element_text(size = 4), legend.text = element_text(size = 8)) + ggtitle(heatmap.title) #+ NoLegend()               #Make genes smaller font size
    
    heatmap.filename <- paste("coexpression_heatmap_",toString(geneName),"_exp",toString(expression.cutoff),"_originalclusters.pdf", sep="")
    ggsave(heatmap.filename, limitsize = FALSE, units = "px", width = 4000, height = 2000)
    
    }
}

dev.off()


############################# Heatmaps, IRs

nrow(pseudo_bulk.IRs)

for (c in 1:nrow(pseudo_bulk.IRs)) {
  #c <- 50
  geneName <- toString(pseudo_bulk.IRs$Gene[c])
  #geneName <- toString(rna_assay_genes[c])
  print(paste(toString(c), geneName))
  
  expression.cutoff <- 0.5
  expr <- FetchData(object = AntennaNeurons, vars = geneName)
  dat <- AntennaNeurons[, which(x = expr > expression.cutoff)]
  print(dim(dat@meta.data)[1])
  
  if (dim(dat@meta.data)[1] > 100) {
    
    dat <- RunPCA(dat, npcs = 50, verbose = F)
    dat <- FindNeighbors(dat, dims = 1:50, verbose = F)
    dat <- FindClusters(dat, resolution = 3)
    #Idents(object = dat) <- "0"
    
    ae <- AverageExpression(dat)                            #get avg expression of all genes in cluster       
    ae <- ae$RNA                                        #pull out only relevant data                               
    ae <- ae[order(ae[,1],decreasing=TRUE),]        #sort from largest to smallest                          
    ae <- as.data.frame(ae)                             #convert to df      
    ae <- tibble::rownames_to_column(ae, "Gene")        #make gene names a column                           
    ae <- ae[!grepl("LOC*", ae$Gene),]              #get rid of LOC & MC genes                     
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
    
    ae <- ae[!grepl("Orco", ae$Gene),]              #Delete co-receptors from main variable                     
    ae <- ae[!grepl("Ir25a", ae$Gene),]
    ae <- ae[!grepl("Ir76b", ae$Gene),]
    ae <- ae[!grepl("Ir8a", ae$Gene),] 
    
    ae <- rbind(x, ae)                                  #Add co-receptors to top     
    
    heatmap.title = paste(toString(geneName), " cells (Normalized Expression>",toString(expression.cutoff),", reclustered)",
                          "\n",toString(dim(dat@meta.data)[1])," cells",sep="")
    #heatmap.title = paste(toString(geneName), " cells (Normalized Expression>",toString(expression.cutoff),", original clusters)",
    #                      "\n",toString(dim(dat@meta.data)[1])," cells",sep="")
    #heatmap.title = paste("rna_",toString(geneName), " cells (Normalized Expression>",toString(expression.cutoff),", original clusters)",
    #                      "\n",toString(dim(dat@meta.data)[1])," cells",sep="")
    DoHeatmap(object = dat, features = ae[,1], label=TRUE, group.bar.height =0.001, draw.lines=FALSE, size=2, angle=90) +            #Make dotplot on sorted gene list
      theme(axis.text.y = element_text(size = 4), legend.text = element_text(size = 8)) + ggtitle(heatmap.title) #+ NoLegend()               #Make genes smaller font size
    
    heatmap.filename <- paste("coexpression_heatmap_",toString(geneName),"_exp",toString(expression.cutoff),"_reclustered.pdf", sep="")
    ggsave(heatmap.filename, limitsize = FALSE, units = "px", width = 4000, height = 2000)
    
  }
  else { 
    #dat <- RunPCA(dat, npcs = 50, verbose = F)
    dat <- FindNeighbors(dat, dims = 1:50, verbose = F)
    #dat <- FindClusters(dat, resolution = 3)
    #Idents(object = dat) <- "0"
    
    ae <- AverageExpression(dat)                            #get avg expression of all genes in cluster       
    ae <- ae$RNA                                        #pull out only relevant data                               
    ae <- ae[order(ae[,1],decreasing=TRUE),]        #sort from largest to smallest                          
    ae <- as.data.frame(ae)                             #convert to df      
    ae <- tibble::rownames_to_column(ae, "Gene")        #make gene names a column                           
    ae <- ae[!grepl("LOC*", ae$Gene),]              #get rid of LOC & MC genes                     
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
    
    ae <- ae[!grepl("Orco", ae$Gene),]              #Delete co-receptors from main variable                     
    ae <- ae[!grepl("Ir25a", ae$Gene),]
    ae <- ae[!grepl("Ir76b", ae$Gene),]
    ae <- ae[!grepl("Ir8a", ae$Gene),] 
    
    ae <- rbind(x, ae)                                  #Add co-receptors to top     
    
    #heatmap.title = paste(toString(geneName), " cells (Normalized Expression>",toString(expression.cutoff),", reclustered)",
    #                      "\n",toString(dim(dat@meta.data)[1])," cells",sep="")
    heatmap.title = paste(toString(geneName), " cells (Normalized Expression>",toString(expression.cutoff),", original clusters)",
                          "\n",toString(dim(dat@meta.data)[1])," cells",sep="")
    #heatmap.title = paste("rna_",toString(geneName), " cells (Normalized Expression>",toString(expression.cutoff),", original clusters)",
    #                      "\n",toString(dim(dat@meta.data)[1])," cells",sep="")
    DoHeatmap(object = dat, features = ae[,1], label=TRUE, group.bar.height =0.001, draw.lines=FALSE, size=2, angle=90) +            #Make dotplot on sorted gene list
      theme(axis.text.y = element_text(size = 4), legend.text = element_text(size = 8)) + ggtitle(heatmap.title) #+ NoLegend()               #Make genes smaller font size
    
    heatmap.filename <- paste("coexpression_heatmap_",toString(geneName),"_exp",toString(expression.cutoff),"_originalclusters.pdf", sep="")
    ggsave(heatmap.filename, limitsize = FALSE, units = "px", width = 4000, height = 2000)
    
  }
}

dev.off()


############################# Heatmaps, GRs

nrow(pseudo_bulk.GRs)

for (c in 1:nrow(pseudo_bulk.GRs)) {
  #c <- 50
  geneName <- toString(pseudo_bulk.GRs$Gene[c])
  #geneName <- toString(rna_assay_genes[c])
  print(paste(toString(c), geneName))
  
  expression.cutoff <- 0.5
  expr <- FetchData(object = AntennaNeurons, vars = geneName)
  dat <- AntennaNeurons[, which(x = expr > expression.cutoff)]
  print(dim(dat@meta.data)[1])
  
  if (dim(dat@meta.data)[1] > 100) {
    
    dat <- RunPCA(dat, npcs = 50, verbose = F)
    dat <- FindNeighbors(dat, dims = 1:50, verbose = F)
    dat <- FindClusters(dat, resolution = 3)
    #Idents(object = dat) <- "0"
    
    ae <- AverageExpression(dat)                            #get avg expression of all genes in cluster       
    ae <- ae$RNA                                        #pull out only relevant data                               
    ae <- ae[order(ae[,1],decreasing=TRUE),]        #sort from largest to smallest                          
    ae <- as.data.frame(ae)                             #convert to df      
    ae <- tibble::rownames_to_column(ae, "Gene")        #make gene names a column                           
    ae <- ae[!grepl("LOC*", ae$Gene),]              #get rid of LOC & MC genes                     
    ae <- ae[!grepl("MT*", ae$Gene),]
    ae.ORs <- ae[grepl("Or*", ae$Gene),]
    ae.IRs <- ae[grepl("Ir*", ae$Gene),]
    ae.GRs <- ae[grepl("Gr*", ae$Gene),]
    #ae <- rbind(ae.ORs, ae.IRs, ae.GRs)
    #ae <- rbind(ae.IRs, ae.ORs, ae.GRs)
    ae <- rbind(ae.GRs, ae.ORs, ae.IRs)
    
    x <- ae[grepl("Orco", ae$Gene),]                    #Pull out co-receptors               
    x[2,] <- ae[grepl("Ir25a", ae$Gene),]
    x[3,] <- ae[grepl("Ir76b", ae$Gene),]
    x[4,] <- ae[grepl("Ir8a", ae$Gene),]
    
    ae <- ae[!grepl("Orco", ae$Gene),]              #Delete co-receptors from main variable                     
    ae <- ae[!grepl("Ir25a", ae$Gene),]
    ae <- ae[!grepl("Ir76b", ae$Gene),]
    ae <- ae[!grepl("Ir8a", ae$Gene),] 
    
    ae <- rbind(x, ae)                                  #Add co-receptors to top     
    
    heatmap.title = paste(toString(geneName), " cells (Normalized Expression>",toString(expression.cutoff),", reclustered)",
                          "\n",toString(dim(dat@meta.data)[1])," cells",sep="")
    #heatmap.title = paste(toString(geneName), " cells (Normalized Expression>",toString(expression.cutoff),", original clusters)",
    #                      "\n",toString(dim(dat@meta.data)[1])," cells",sep="")
    #heatmap.title = paste("rna_",toString(geneName), " cells (Normalized Expression>",toString(expression.cutoff),", original clusters)",
    #                      "\n",toString(dim(dat@meta.data)[1])," cells",sep="")
    DoHeatmap(object = dat, features = ae[,1], label=TRUE, group.bar.height =0.001, draw.lines=FALSE, size=2, angle=90) +            #Make dotplot on sorted gene list
      theme(axis.text.y = element_text(size = 4), legend.text = element_text(size = 8)) + ggtitle(heatmap.title) #+ NoLegend()               #Make genes smaller font size
    
    heatmap.filename <- paste("coexpression_heatmap_",toString(geneName),"_exp",toString(expression.cutoff),"_reclustered.pdf", sep="")
    ggsave(heatmap.filename, limitsize = FALSE, units = "px", width = 4000, height = 2000)
    
  }
  else { 
    #dat <- RunPCA(dat, npcs = 50, verbose = F)
    dat <- FindNeighbors(dat, dims = 1:50, verbose = F)
    #dat <- FindClusters(dat, resolution = 3)
    #Idents(object = dat) <- "0"
    
    ae <- AverageExpression(dat)                            #get avg expression of all genes in cluster       
    ae <- ae$RNA                                        #pull out only relevant data                               
    ae <- ae[order(ae[,1],decreasing=TRUE),]        #sort from largest to smallest                          
    ae <- as.data.frame(ae)                             #convert to df      
    ae <- tibble::rownames_to_column(ae, "Gene")        #make gene names a column                           
    ae <- ae[!grepl("LOC*", ae$Gene),]              #get rid of LOC & MC genes                     
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
    
    ae <- ae[!grepl("Orco", ae$Gene),]              #Delete co-receptors from main variable                     
    ae <- ae[!grepl("Ir25a", ae$Gene),]
    ae <- ae[!grepl("Ir76b", ae$Gene),]
    ae <- ae[!grepl("Ir8a", ae$Gene),] 
    
    ae <- rbind(x, ae)                                  #Add co-receptors to top     
    
    #heatmap.title = paste(toString(geneName), " cells (Normalized Expression>",toString(expression.cutoff),", reclustered)",
    #                      "\n",toString(dim(dat@meta.data)[1])," cells",sep="")
    heatmap.title = paste(toString(geneName), " cells (Normalized Expression>",toString(expression.cutoff),", original clusters)",
                          "\n",toString(dim(dat@meta.data)[1])," cells",sep="")
    #heatmap.title = paste("rna_",toString(geneName), " cells (Normalized Expression>",toString(expression.cutoff),", original clusters)",
    #                      "\n",toString(dim(dat@meta.data)[1])," cells",sep="")
    DoHeatmap(object = dat, features = ae[,1], label=TRUE, group.bar.height =0.001, draw.lines=FALSE, size=2, angle=90) +            #Make dotplot on sorted gene list
      theme(axis.text.y = element_text(size = 4), legend.text = element_text(size = 8)) + ggtitle(heatmap.title) #+ NoLegend()               #Make genes smaller font size
    
    heatmap.filename <- paste("coexpression_heatmap_",toString(geneName),"_exp",toString(expression.cutoff),"_originalclusters.pdf", sep="")
    ggsave(heatmap.filename, limitsize = FALSE, units = "px", width = 4000, height = 2000)
    
  }
}

dev.off()





