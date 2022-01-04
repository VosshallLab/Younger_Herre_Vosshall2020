# Table of Contents (And Where to Find Large Data Files)
Updated January 3, 2022 - for questions, contact ogoldman@rockefeller.edu  

## Repositories
All supplementary analysis is in the following repository:
https://github.com/VosshallLab/Younger_Herre_Vosshall2020

Antenna snRNAseq data files, scripts, and expression matricies are on Zenodo:
https://zenodo.org/record/5818543

Maxillary Palp snRNAseq data files, scripts, and expression matricies are on Zenodo:
https://zenodo.org/record/5818952

For raw sequencing files for both antenna and maxillary palp see NCBI/GEO number: GSE192978

## *Antenna: List of Data Items*
| TISSUE 	| LOCATION 	| ITEM DESCRIPTION 	| DATA FILE NAME (all on [Zenodo](https://zenodo.org/record/5818543)) 	| SCRIPT FILE NAME (all on [Zenodo](https://zenodo.org/record/5818543)	|
|---	|---	|---	|---	|---	|
| Antenna 	| Figure 5B 	| All Cells Marker Heatmap 	| antenna_6.neuralClusters.Rdata 	| antenna_7.expressedReceptors.R 	|
| Antenna 	| Figure 5C 	| Neurons tSNE 	| antenna_6.neuralClusters.Rdata 	| antenna_7.expressedReceptors.R 	|
| Antenna 	| Figure 5D 	| Dot plot (manually edited, for full version see Figure S10A) 	| N/A 	| N/A 	|
| Antenna 	| Figure 5E 	| Chord Plot 	| antenna_6.neuralClusters.Rdata 	| antenna_7.expressedReceptors.R 	|
| Antenna 	| Figure 5F-I 	| Simplified Coexpression Heatmap (Ir64a, Or4, Or82, I41k) 	| SeuratObject_AntennaNeurons.RData 	| SimplifiedCoexpressionHeatmaps.R 	|
| Antenna 	| Figure S9B 	| Doublet Removal Scatter Plot 	| antenna_2.DoubletFinder.08192021.rds antenna_2.DoubletFinder.10062021.rds 	| antenna_2.doubletfinder.R 	|
| Antenna 	| Figure S9C-E 	| Genes, UMIs, Mitochondrial Genes Violin Plots 	| antenna_3.merge.filter.rds 	| antenna_3.integrate.R 	|
| Antenna 	| Figure S9F 	| Batch Correction Feature Plot 	| antenna_3.merge.filter.rds antenna_4.seurat.combined.rds 	| antenna_3.integrate.R antenna_4.anchor.R 	|
| Antenna 	| Figure S9G 	| All Cells tSNE 	| antenna_4.seurat.combined.rds 	| antenna_4.anchor.R 	|
| Antenna 	| Figure S9H 	| Cell Marker Feature Plot 	| SeuratObject_AntennaNeurons.RData 	| FeaturePlots.R 	|
| Antenna 	| Figure S9I 	| Neural Marker Dot Plot 	| SeuratObject_AntennaNeurons.RData 	| FeaturePlots.R 	|
| Antenna 	| Figure S10A 	| Chemosensor Dot Plot 	| antenna_6.neuralClusters.Rdata 	| antenna_7.expressedReceptors.R 	|
| Antenna 	| Figure S10B 	| Dot Plot Criteria Scatter Plot 	| antenna_6.neuralClusters.Rdata 	| antenna_7.expressedReceptors.R 	|
| Antenna 	| Figure S10C 	| Cluster Coexpression Heatmap 	| SeuratObject_AntennaNeurons.RData 	| CoexpressionHeatmaps.R 	|
| Antenna 	| Figure S10D 	| Gene Expression Violin Plot (Orco, Ir25a) 	| SeuratObject_AntennaNeurons.RData 	| ViolinPlots.R 	|
| Antenna 	| Figure S10E 	| Cluster Gene Expression Scatter Plot 	| antenna_6.neuralClusters.Rdata 	| antenna_7.expressedReceptors.R 	|
| Antenna 	| Figure S10F 	| Cell Gene Expression Scatter Plot 	| antenna_6.neuralClusters.Rdata 	| antenna_7.expressedReceptors.R 	|
| Antenna 	| Figure S10G 	| Gene Expression Venn Diagram 	| antenna_6.neuralClusters.Rdata 	| antenna_7.expressedReceptors.R 	|
| Antenna 	| Figure S11A 	| Gene Expression Violin Plot (Ir41a, Or82, Or4, Ir64a, Gr77, Ir8a, Or47,   Ir41k) 	| SeuratObject_AntennaNeurons.RData 	| ViolinPlots.R 	|
| Antenna 	| Figure S11B 	| Coexpression Heatmap Ir41k 	| SeuratObject_AntennaNeurons.RData 	| CoexpressionHeatmaps.R 	|
| Antenna 	| Figure S11C-F 	| Simplified Coexpression Heatmap (Ir41a, Ir8a, Gr77, Or47) 	| SeuratObject_AntennaNeurons.RData 	| SimplifiedCoexpressionHeatmaps.R 	|
| Antenna 	| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/PseudobulkTables) 	| PseudoBulk Tables 	| SeuratObject_AntennaNeurons.RData 	| PseudoBulkTables.R 	|
| Antenna 	| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/Antenna_FeaturePlots) 	| Gene Expression Feature Plots 	| SeuratObject_AntennaNeurons.RData 	| FeaturePlots.R 	|
| Antenna 	| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/Antenna_ViolinPlots) 	| Gene Expression Violin Plots 	| SeuratObject_AntennaNeurons.RData 	| ViolinPlots.R 	|
| Antenna 	| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/Antennae_BroadCoexpressionHeatmaps_byOriginalClusters) 	| Broad Coexpression Heatmaps by Neuron Cluster 	| SeuratObject_AntennaNeurons.RData 	| CoexpressionHeatmaps.R 	|
| Antenna 	| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/Antenna_BroadCoexpressionHeatmaps_byChemoreceptor) 	| Broad Coexpression Heatmaps by Chemoreceptor 	| SeuratObject_AntennaNeurons.RData 	| CoexpressionHeatmaps.R 	|
| Antenna 	| [Zenodo](https://zenodo.org/record/5818543) 	| Raw Expression Matrix (Batch 1    & 2) 	| antenna1_raw_feature_bc_matrix.h5 antenna2_raw_feature_bc_matrix.h5 	| N/A 	|
| Antenna 	| [Zenodo](https://zenodo.org/record/5818543) 	| Filtered Expression Matrix (Batch 1    & 2) 	| antenna1_filtered_feature_bc_matrix.h5 antenna2_filtered_feature_bc_matrix.h5 	| N/A 	|
| Antenna 	| [Zenodo](https://zenodo.org/record/5818543) 	| Seurat Object: All Antenna Cells (Batch 1 & 2 combined) 	| SeuratObject_AntennaAllCells.rds 	| N/A 	|
| Antenna 	| [Zenodo](https://zenodo.org/record/5818543) 	| Seurat Object: All Antenna Neurons (Batch 1 & 2 combined) 	| SeuratObject_AntennaNeurons.rds 	| N/A 	|


## *Maxillary Palp: List of Data Items*
| TISSUE 	| LOCATION 	| ITEM DESCRIPTION 	| DATA FILE NAME (all on [Zenodo](https://zenodo.org/record/5818952)) 	| SCRIPT FILE NAME (all on [Zenodo](https://zenodo.org/record/5818952)	|
|---	|---	|---	|---	|---	|
| Maxillary Palp 	| Figure 7B 	| Cell Marker Heatmap  	| maxpalp_4.neuralClusters.RData 	| maxpalp_5.expressedReceptors.R 	|
| Maxillary Palp 	| Figure 7C 	| tSNE 	| maxpalp_3.neuron.Rdata 	| maxpalp_3.check_receptor.R 	|
| Maxillary Palp 	| Figure 7D 	| Neuron Gene Heatmap 	| maxpalp_4.neuralClusters.RData 	| maxpalp_5.expressedReceptors.R 	|
| Maxillary Palp 	| Figure 7E 	| Chord Plot 	| maxpalp_4.neuralClusters.RData 	| maxpalp_5.expressedReceptors.R 	|
| Maxillary Palp 	| Figure 7F-M 	| Feature Plot 	| SeuratObject_MaxPalpAllCells.RData 	| FeaturePlots.R 	|
| Maxillary Palp 	| Figure S13A 	| Orco Feature Plot (DeconX illustration) 	| maxpalp_1.Decontx.08192021.Rdata 	| maxpalp_1.decontX.R 	|
| Maxillary Palp 	| Figure S13B 	| Doublet Removal Scatter Plot 	| maxpalp_2.DoubletFinder.rds 	| maxpalp_2.doubletfinder.R 	|
| Maxillary Palp 	| Figure S13C-E 	| Genes, UMIs, Mitochondrial Genes Violin Plots 	| SeuratObject_MaxPalpAllCells.RData 	| ViolinPlots.R 	|
| Maxillary Palp 	| Figure S13F 	| Cell Marker Feature Plot 	| SeuratObject_MaxPalpAllCells.RData 	| FeaturePlots.R 	|
| Maxillary Palp 	| Figure S13G 	| Dot Plot Criteria Scatter Plot 	| maxpalp_4.neuralClusters.RData 	| maxpalp_5.expressedReceptors.R 	|
| Maxillary Palp 	| Figure S13H 	| Neural Marker Dot Plot 	| SeuratObject_MaxPalpAllCells.RData 	| FeaturePlots.R 	|
| Maxillary Palp 	| Figure S14A-C 	| Genes, UMIs, Mitochondrial Genes Violin Plots (Neuron Clusters) 	| SeuratObject_MaxPalpNeurons.RData 	| ViolinPlots.R 	|
| Maxillary Palp 	| Figure S14D 	| Gene Expression Violin Plot by Cluster (Gr1, Gr2, Gr3, Ir25a, Orco, Or8,   Or49, Ir76b) 	| SeuratObject_MaxPalpNeurons.RData 	| ViolinPlots.R 	|
| Maxillary Palp 	| Figure S14E 	| Neuron Cluster Marker Dot Plot 	| SeuratObject_MaxPalpNeurons.RData 	| FeaturePlots.R 	|
| Maxillary Palp 	| Figure S14F-G 	| Cell Gene Expression Scatter Plot 	| maxpalp_4.neuralClusters.RData 	| maxpalp_5.expressedReceptors.R 	|
| Maxillary Palp 	| Figure S14H 	| Gene Expression Feature Plot 	| SeuratObject_MaxPalpAllCells.RData 	| FeaturePlots.R 	|
| Maxillary Palp 	| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/PseudobulkTables) 	| PseudoBulk Tables 	| SeuratObject_MaxPalpAllCells.RData 	| PseudoBulkTables.R 	|
| Maxillary Palp 	| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/MaxPalp_FeaturePlots) 	| Gene Expression Feature Plots 	| SeuratObject_MaxPalpAllCells.RData 	| FeaturePlots.R 	|
| Maxillary Palp 	| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/MaxPalp_ViolinPlots) 	| Gene Expression Violin Plots 	| SeuratObject_MaxPalpNeurons.RData 	| ViolinPlots.R 	|
| Maxillary Palp 	| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/MaxPalp_ViolinPlots_byCluster) 	| Gene Expression Violin Plots by Cluster 	| SeuratObject_MaxPalpNeurons.RData 	| ViolinPlots.R 	|
| Maxillary Palp 	| [Zenodo](https://zenodo.org/record/581895) 	| Raw Expression Matrix 	| maxpalp_raw_feature_bc_matrix.h5 	| N/A 	|
| Maxillary Palp 	| [Zenodo](https://zenodo.org/record/5818952) 	| Filtered Expression Matrix 	| maxpalp_filtered_feature_bc_matrix.h5 	| N/A 	|
| Maxillary Palp 	| [Zenodo](https://zenodo.org/record/5818952) 	| Seurat   Object: All Maxillary Palp Cells 	| SeuratObject_MaxPalpAllCells.RData 	| N/A 	|
| Maxillary Palp 	| [Zenodo](https://zenodo.org/record/5818952) 	| Seurat   Object: All Maxillary Palp Neurons 	| SeuratObject_MaxPalpNeurons.RData 	| N/A 	|
