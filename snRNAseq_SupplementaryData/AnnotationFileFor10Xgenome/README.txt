Annotation file used with AaegL5.0 genome for 10X reference genome generation: 20190215_BRCTSC_Sup18_ref_AaegL5.0_top_level_forPipe_rmtRNA_CDS_notxid.gtf

Modified from https://github.com/VosshallLab/Jove_Vosshall_2020/tree/master/RNAseq_merged_annotation
Contains manually annotated chemoreceptors (from NCBI_annotation_UPDATED_chemoreceptors_CORRECTED.gtf)

Contains ncRNA (example: LOC23687648)

Additional formatted for scRNAseq:
- mitochondrial genes have been appended with "MT-" (for Seurat handling)
- mitochondrial genes have been editted to be "Exon" rather than "CDS"
- tRNAs are removed (to avoid issues in cellranger)

