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
  library(stats)
  library(gridExtra)
  library(scater)
  library(celda)
  library(DropletUtils)
  library(harmony)
  library(DoubletFinder)
  library(ggrepel)
  library(ggplot2)
  library(reshape2)
  library(biomaRt)
})


#.....................................................................................
# SET UP
# ....................................................................................
rm(list = ls())   #remove all objects
setwd("E:/Github/Zenodo-Antenna")
antenna.neuron <- readRDS("9_neurons_harmony_rm0_res4.rds")  #Uploaded on Zenodo as "SeuratObject2_Antenna_mergedBatches_Neurons.rds"
pipeline <- "harmony_rmclusts_res4"
dat.assay <- "SCT"
DefaultAssay(antenna.neuron) <- dat.assay

outputFeaturePlot <- 'Antenna_Figures4_NeuronPseudobulk_Output'
dir.create(outputFeaturePlot)
setwd(outputFeaturePlot)

######################################################################################
# Create table of raw counts (RNA assay)
######################################################################################
dat.assay <- "RNA"

#Separate seurat objects
seurat.0819 <- subset(antenna.neuron, subset = batch == "batch1")
dat <- seurat.0819

mtx <- GetAssayData(dat, slot ="counts", assay=dat.assay)
pseudo_bulk <- rowSums(mtx)
#write.csv(pseudo_bulk, "PseudoBulk_RNA_counts_batch1_0819.csv", row.names=TRUE)

df <- as.data.frame(pseudo_bulk)
names(df)[1] <- "Counts"   #Rename column by index
df$Genes <- rownames(df)

df_nompC <- df[grep("LOC5575210", df$Genes),]
df_chemoR <- df[- grep("LOC", df$Genes),]
df_chemoR <- df_chemoR[- grep("MT", df_chemoR$Genes),]
df_chemoR[(nrow(df_chemoR)+1),] <- df_nompC[1,]

df_chemoR.0819 <- as.data.frame(df_chemoR)
p <- ggplot(df_chemoR.0819, aes(x= reorder(Genes, -Counts), y=Counts)) + geom_point()
p <- p + scale_y_continuous(trans='log10')
p <- p + ggtitle("Chemoreceptor expression (PseudoBulk)") + xlab("Receptor Gene") + ylab("Counts (log10)")
p1 <- p + theme(axis.text.x = element_text(angle = 90, size = 12))
p1

rm(mtx, dat, df, df_chemoR, df_nompC, p, p1, seurat.0819)

# ....................................................................................
seurat.1006 <- subset(antenna.neuron, subset = batch == "batch2")
dat <- seurat.1006

mtx <- GetAssayData(dat, slot = "counts", assay=dat.assay)
pseudo_bulk <- rowSums(mtx)
#write.csv(pseudo_bulk, "PseudoBulk_RNA_counts_batch2_1006.csv", row.names=TRUE)

df <- as.data.frame(pseudo_bulk)
names(df)[1] <- "Counts"   #Rename column by index
df$Genes <- rownames(df)

df_nompC <- df[grep("LOC5575210", df$Genes),]
df_chemoR <- df[- grep("LOC", df$Genes),]
df_chemoR <- df_chemoR[- grep("MT", df_chemoR$Genes),]
df_chemoR[(nrow(df_chemoR)+1),] <- df_nompC[1,]

df_chemoR.1006 <- as.data.frame(df_chemoR)
p <- ggplot(df_chemoR.1006, aes(x= reorder(Genes, -Counts), y=Counts)) + geom_point()
p <- p + scale_y_continuous(trans='log10')
p <- p + ggtitle("Chemoreceptor expression (PseudoBulk)") + xlab("Receptor Gene") + ylab("Counts (log10)")
p2 <- p + theme(axis.text.x = element_text(angle = 90, size = 12))
p2

p1+p2

b1_b2 <- merge(df_chemoR.0819,df_chemoR.1006,by="Genes")
colnames(b1_b2)
names(b1_b2)[names(b1_b2) == 'Counts.x'] <- paste0('batch1.0819.',dat.assay)
names(b1_b2)[names(b1_b2) == 'Counts.y'] <- paste0('batch2.1006.',dat.assay)
#write.csv(b1_b2, "PseudoBulk_RNA_counts.csv", row.names=TRUE)

# ....................................................................................

b1_b2.RNA <- b1_b2

rm(mtx, dat, df, df_chemoR, df_nompC, p, p1, p2, seurat.1006, df_chemoR.0819, df_chemoR.1006)

######################################################################################
# Create table of raw counts (SCT assay)
######################################################################################
dat.assay <- "SCT"

#Separate seurat objects
seurat.0819 <- subset(antenna.neuron, subset = batch == "batch1")
dat <- seurat.0819

mtx <- GetAssayData(dat, slot ="counts", assay=dat.assay)
pseudo_bulk <- rowSums(mtx)
#write.csv(pseudo_bulk, "PseudoBulk_RNA_counts_batch1_0819.csv", row.names=TRUE)

df <- as.data.frame(pseudo_bulk)
names(df)[1] <- "Counts"   #Rename column by index
df$Genes <- rownames(df)

df_nompC <- df[grep("LOC5575210", df$Genes),]
df_chemoR <- df[- grep("LOC", df$Genes),]
df_chemoR <- df_chemoR[- grep("MT", df_chemoR$Genes),]
df_chemoR[(nrow(df_chemoR)+1),] <- df_nompC[1,]

df_chemoR.0819 <- as.data.frame(df_chemoR)
p <- ggplot(df_chemoR.0819, aes(x= reorder(Genes, -Counts), y=Counts)) + geom_point()
p <- p + scale_y_continuous(trans='log10')
p <- p + ggtitle("Chemoreceptor expression (PseudoBulk)") + xlab("Receptor Gene") + ylab("Counts (log10)")
p1 <- p + theme(axis.text.x = element_text(angle = 90, size = 12))
p1

rm(mtx, dat, df, df_chemoR, df_nompC, p, p1, seurat.0819)

# ....................................................................................
seurat.1006 <- subset(antenna.neuron, subset = batch == "batch2")
dat <- seurat.1006

mtx <- GetAssayData(dat, slot = "counts", assay=dat.assay)
pseudo_bulk <- rowSums(mtx)
#write.csv(pseudo_bulk, "PseudoBulk_RNA_counts_batch2_1006.csv", row.names=TRUE)

df <- as.data.frame(pseudo_bulk)
names(df)[1] <- "Counts"   #Rename column by index
df$Genes <- rownames(df)

df_nompC <- df[grep("LOC5575210", df$Genes),]
df_chemoR <- df[- grep("LOC", df$Genes),]
df_chemoR <- df_chemoR[- grep("MT", df_chemoR$Genes),]
df_chemoR[(nrow(df_chemoR)+1),] <- df_nompC[1,]

df_chemoR.1006 <- as.data.frame(df_chemoR)
p <- ggplot(df_chemoR.1006, aes(x= reorder(Genes, -Counts), y=Counts)) + geom_point()
p <- p + scale_y_continuous(trans='log10')
p <- p + ggtitle("Chemoreceptor expression (PseudoBulk)") + xlab("Receptor Gene") + ylab("Counts (log10)")
p2 <- p + theme(axis.text.x = element_text(angle = 90, size = 12))
p2

p1+p2

b1_b2 <- merge(df_chemoR.0819,df_chemoR.1006,by="Genes")
colnames(b1_b2)
names(b1_b2)[names(b1_b2) == 'Counts.x'] <- paste0('batch1.0819.',dat.assay)
names(b1_b2)[names(b1_b2) == 'Counts.y'] <- paste0('batch2.1006.',dat.assay)
#write.csv(b1_b2, "PseudoBulk_RNA_counts.csv", row.names=TRUE)
b1_b2.SCT <- b1_b2

rm(mtx, dat, df, df_chemoR, df_nompC, p, p1, p2, seurat.1006, df_chemoR.0819, df_chemoR.1006)

# ....................................................................................
b1_b2.SCT <- b1_b2

b1_b2 <- merge(b1_b2.RNA,b1_b2.SCT,by="Genes")


######################################################################################
# Combined counts from both batches for SCT and RNA assays
######################################################################################

dat <- antenna.neuron
mtx <- GetAssayData(dat, slot = "counts", assay="SCT")
pseudo_bulk <- rowSums(mtx)

df <- as.data.frame(pseudo_bulk)
names(df)[1] <- "Counts"   #Rename column by index
df$Genes <- rownames(df)

df_nompC <- df[grep("LOC5575210", df$Genes),]
df_chemoR <- df[- grep("LOC", df$Genes),]
df_chemoR2 <- df_chemoR[- grep("MT", df_chemoR$Genes),]
df_chemoR2[(nrow(df_chemoR2)+1),] <- df_nompC[1,]

df_chemoR_sort <- df_chemoR2[order(-df_chemoR2$Counts),]

df_chemoR2.antenna <- df_chemoR2

all_counts <- merge(b1_b2,df_chemoR2.antenna,by="Genes")
names(all_counts)[names(all_counts) == 'Counts'] <- 'batchsum_SCT'


# ....................................................................................
dat <- antenna.neuron
mtx <- GetAssayData(dat, slot = "counts", assay="RNA")
pseudo_bulk <- rowSums(mtx)

df <- as.data.frame(pseudo_bulk)
names(df)[1] <- "Counts"   #Rename column by index
df$Genes <- rownames(df)

df_nompC <- df[grep("LOC5575210", df$Genes),]
df_chemoR <- df[- grep("LOC", df$Genes),]
df_chemoR2 <- df_chemoR[- grep("MT", df_chemoR$Genes),]
df_chemoR2[(nrow(df_chemoR2)+1),] <- df_nompC[1,]

df_chemoR_sort <- df_chemoR2[order(-df_chemoR2$Counts),]

df_chemoR2.antenna <- df_chemoR2

all_counts <- merge(all_counts,df_chemoR2.antenna,by="Genes")
names(all_counts)[names(all_counts) == 'Counts'] <- 'batchsum_RNA'

rm(b1_b2, mtx, dat, df, df_chemoR, df_chemoR2, 
 df_nompC,df_chemoR_sort, df_chemoR2.antenna)

write.csv(all_counts, "Antenna_pseudoBulk_chemoreceptor_counts.csv", row.names=TRUE)

######################################################################################
# FIGURE S5F
######################################################################################

theme_set(theme_cowplot())

dat <- all_counts[, c("Genes","batch1.0819.RNA", "batch2.1006.RNA", "batchsum_SCT")]
dat_sort <- dat[order(-dat$batchsum_SCT),]

dat_sort$Genes <- as.character(dat_sort$Genes)
dat_sort$Genes <- factor(dat_sort$Genes, levels=c(dat_sort$Genes))

df <- subset(dat_sort, select = -c(batchsum_SCT))
d <- melt(df, id.vars="Genes")

# Everything on the same plot
p <- ggplot(d, aes(Genes,value, col=variable)) + 
  geom_point(size=1) + scale_y_continuous(trans='log10')
p <- p + ggtitle("Chemoreceptor expression (PseudoBulk)") + xlab("Receptor Gene") + ylab("Counts (log10)")
p <- p + theme(axis.text.x = element_text(angle =90, size = 6))
p

filename <- paste0("15_pseudobulk_allgenes.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 2000, p)


#For only top 93 genes
df.abbrev<-df[!(df$Gene=="LOC5575210"),]

df.abbrev <- df.abbrev[1:93,]
d <- melt(df.abbrev, id.vars="Genes")

# Everything on the same plot
p <- ggplot(d, aes(Genes,value, col=variable)) + 
  geom_point(size=1) + scale_y_continuous(trans='log10')
p <- p + ggtitle("Chemoreceptor expression (PseudoBulk)") + xlab("Receptor Gene") + ylab("Counts (log10)")
p <- p + theme(axis.text.x = element_text(angle =90, size = 6))
p

filename <- paste0("15_pseudobulk_top93.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 2000, p)




