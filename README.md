# burn-skin-cell-type-timeline
Provided by Xueling Li, email: xlli@cmpt.ac.cn
Codes and CIBERSORTx group-mode deconvolution results (GSE8056_group) on GSE8056 dataset for burn-skin-cell-type-timeline reprository. 
Codes:
dataAugm.R is for data augmentation when the number of bulk samples is smaller than the number of cell types;
genename_countoTPM.R transform raw counts data to TPM format.
getGEO.R provides demo for microarray downloading and gene ID mapping for both microarray and RNA-Seq data.
scRNA_GSE169147.R is the preprocessing scRNA-Seq dataset of GSE169147, sampling and preparing data for uploading to CIBERSORTx for cell type signature construction and group-mode deconvolution.
scRNA_GSE174072.R is the preprocessing scRNA-Seq dataset of GSE174072, sampling and preparing, and combinding with GSE169147 before uploading to CIBERSORTx for cell type signature construction and group-mode deconvolution.
simuSkin_ori.R provides the pipeline for differential expressed genes identification for each deconvoluted cell type.
