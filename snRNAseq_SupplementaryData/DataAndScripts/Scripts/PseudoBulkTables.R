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

counts <- dataset@assays$RNA@counts                            #Pull raw counts from Seurat object

pseudo_bulk <- rowSums(counts)                                             #Sum values across cells

pseudo_bulk <- as.data.frame(pseudo_bulk)                                  #Convert to dataframe
pseudo_bulk <- tibble::rownames_to_column(pseudo_bulk, "Gene")             #Pull out rownames as column

pseudo_bulk <- pseudo_bulk[!grepl("LOC*", pseudo_bulk$Gene),]              #get rid of LOC & MC genes (for just chemoreceptor analysis)
pseudo_bulk <- pseudo_bulk[!grepl("MT*", pseudo_bulk$Gene),]                    
pseudo_bulk <- pseudo_bulk[order(pseudo_bulk[,2],decreasing=TRUE),]        #Sort genes in decreasing order 

#pseudo_bulk.ORs <- pseudo_bulk[grepl("Or*", pseudo_bulk$Gene),]            #Separate ORs, IRs, GRs into three variables
#pseudo_bulk.IRs <- pseudo_bulk[grepl("Ir*", pseudo_bulk$Gene),]
#pseudo_bulk.GRs <- pseudo_bulk[grepl("Gr*", pseudo_bulk$Gene),]

#head(pseudo_bulk.ORs)

#view(x)

pseudo_bulk_normalized <- pseudo_bulk
pseudo_bulk_raw <- pseudo_bulk

normalizedData <- dataset@assays$RNA@data
rawData <- dataset@assays$RNA@counts

pseudo_bulk_normalized$NormalizedData_CellCount_0.5 <- NA 
pseudo_bulk_normalized$NormalizedData_CellCount_1.0 <- NA 
pseudo_bulk_normalized$NormalizedData_CellCount_1.5 <- NA 
pseudo_bulk_normalized$NormalizedData_CellCount_2.0 <- NA 

pseudo_bulk_raw$RawData_CellCount_0.5 <- NA 
pseudo_bulk_raw$RawData_CellCount_1.0 <- NA 
pseudo_bulk_raw$RawData_CellCount_1.5 <- NA 
pseudo_bulk_raw$RawData_CellCount_2.0 <- NA 


######### Normalized

for (c in 1:nrow(pseudo_bulk_normalized)) {                                        #Loop for calculating number of cells expressing chemoreceptor genes
  #c <- 1
  geneName <- toString(pseudo_bulk_normalized$Gene[c])                             #Pull gene name 
  print(paste(toString(c), geneName))
  
  geneRow <- subset(normalizedData, rownames(normalizedData) %in% geneName)                                                                #Set expression cut off
  
  pseudo_bulk_normalized$NormalizedData_CellCount_0.5[c] <- sum(geneRow > 0.5) 
  pseudo_bulk_normalized$NormalizedData_CellCount_1.0[c] <- sum(geneRow > 1.0) 
  pseudo_bulk_normalized$NormalizedData_CellCount_1.5[c] <- sum(geneRow > 1.5) 
  pseudo_bulk_normalized$NormalizedData_CellCount_2.0[c] <- sum(geneRow > 2.0) 
  
}

for (c in 1:nrow(pseudo_bulk_raw)) {                                        #Loop for calculating number of cells expressing chemoreceptor genes
  #c <- 1
  geneName <- toString(pseudo_bulk_raw$Gene[c])                             #Pull gene name 
  print(paste(toString(c), geneName))
  
  geneRow <- subset(rawData, rownames(rawData) %in% geneName)                                                                #Set expression cut off
  
  pseudo_bulk_raw$RawData_CellCount_0.5[c] <- sum(geneRow > 0.5) 
  pseudo_bulk_raw$RawData_CellCount_1.0[c] <- sum(geneRow > 1.0) 
  pseudo_bulk_raw$RawData_CellCount_1.5[c] <- sum(geneRow > 1.5) 
  pseudo_bulk_raw$RawData_CellCount_2.0[c] <- sum(geneRow > 2.0) 
  
}

write.csv(pseudo_bulk_normalized, "Antenna_PseudoBulk_NormalizedData.csv", row.names=FALSE)
#write.csv(pseudo_bulk_normalized, "MaxPalp_PseudoBulk_NormalizedData.csv", row.names=FALSE)

write.csv(pseudo_bulk_raw, "Antenna_PseudoBulk_RawData.csv", row.names=FALSE)
#write.csv(pseudo_bulk_raw, "MaxPalp_PseudoBulk_RawData.csv", row.names=FALSE)


