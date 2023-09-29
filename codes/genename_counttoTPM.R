
#setwd("H:\\Validation_dataSet")

prepData<-function(){
library(readxl)
library(openxlsx)
counts_data = read.xlsx("GSE139028_Counts_Table_for_RNA-Seq_Data.xlsx",colNames = T) 
# "GSE139028_Counts_Table_for_RNA-Seq_Data.xlsx" was downloaded from GEO.
genelength = read.table("All_hg19gene_len.txt",header = T,sep = '\t')  # path to "hg19gene_len.txt" is different.
counts_data <- entrz2symbol(counts_data)
#counts_data <- data.frame(Gene = rownames(counts_data),counts_data)
colnames(counts_data)[1] <- "Gene"
counts_data <- merge(genelength,counts_data,by = "Gene")
}
# boxplot(datatpms)

##############fpkmToTpm计算

# fpkmToTpm <- function(fpkm){
#   exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
# }
# datatpms <- apply(data,2,fpkmToTpm)


############TPM计算
counts2tpm <- function(counts_data){
mycounts <- counts_data

mycounts$kb <- mycounts$Length / 1000
rownames(mycounts) <- mycounts[,1]
countdata <- mycounts[,c(-1,-2,-3,-length(mycounts[1,]))]
rpk <- countdata / mycounts$kb
#rpk
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
#head(tpm)
tpm <- data.frame(tpm)
tpm <- data.frame(GeneSymbol = rownames(tpm),tpm)
tpm <- tpm[,-1]
write.table(tpm, file="GS3139028_TPM.txt",sep="\t",quote=F)
return(tpm)
}


######gene name mapping
entrz2symbol <- function(counts_data_input){

# loaded dependent packages
library(clusterProfiler)
library(org.Hs.eg.db)  #org.Mm.eg.db mouse
# check gene names provided by org.Hs.eg.db, which are able to map between themselves
  keytypes(org.Hs.eg.db)

# Transform gene from ENTREZID to SYMBOL
ENTREZ_ID <- counts_data$GeneID
# Transform gene id with bitr()function
genes <- bitr(ENTREZ_ID,fromType="ENTREZID",toType=c("ENTREZID","SYMBOL"),OrgDb="org.Hs.eg.db")
# Check the gene mapping results
genes <- data.frame(Geneid = genes$ENTREZID,SYMBOL=genes$SYMBOL )
counts_data$Geneid <- counts_data$GeneID
counts_data <- merge(genes,counts_data,by = "Geneid")
counts_data=aggregate(.~SYMBOL,mean,data=counts_data[,-1])  ##average of the genes with same symbols.
rownames(counts_data) <- counts_data[,1];
#counts_data <- counts_data[,-1]
return(counts_data)
}



