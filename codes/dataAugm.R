############# Data augmentation #################
# Prepred by Xiaoyi Fei, Xueling Li and Min Zhu
# Contact Xueling Li via email: xlli@cmpt.ac.cn or xuelingli16@foxmail.com
# Generate 200 random values with a fraction of one unit standard deviation of gene expression profiles of x of each group.

rnorm_2 <- function(x) {
  rnorm(200,0,runif(1,min = 0,max = 1)*sd(x))
}

## input of dataAugm:
  #  num_samples, number of samples in the group for the group-wise data augmentation
  #  num_genes, number of genes in each sample.
# HAN, the augmented data for each group
dataAugm <- function(GSE8056_MAS5_groupi){
  
  # num_genes <- dim(exprData)[1]
  # num_samples <- dim(exprData)[2]
  # HAN <- data.frame(matrix(NA,num_genes,num_samples))
  HAN <- GSE8056_MAS5_groupi ## for dataset GSE139028, AE1,4,5 is normal skin, the remaining are burned
  HAN <- log2(HAN+1)
  ### generate 200 samples with the same mean and less than one but more than 0 of one unit of original standard deviation.
  norm1 <- data.frame(t(data.frame(apply(HAN,1,rnorm_2))))   
  

  for(i in 4:200){
    HAN[,i]  <-data.frame(rowMeans(HAN[,1:3]) + norm1[,i-3])
    colnames(HAN)[i] <- paste("C",i,sep="_") ###### for each group's column names, where the first three columns are original data, and the remaining is the augmented data.
  }
  boxplot(as.matrix(HAN))
  HAN <- 2^HAN-1
  HAN <- data.frame(Gene = rownames(HAN),HAN)
  barplot(as.matrix(t(colSums(HAN[,-1]))) )
  HAN_s <- HAN[,-1:-3]
  for (i in 1:15) {
    HAN_s1 <-  HAN_s[,sample(ncol(HAN_s),67,replace = FALSE)]
    HAN_s1 <- data.frame(HAN[,1:3],HAN_s1) 
    write.table(HAN_s1,paste("Path_to_GSE8056\\GSE8056\\Simulated\\D_8_17_",i,".txt",sep = ""),sep ="\t", row.names =F, col.names =TRUE, quote =TRUE)
  }
}

function(){
augmented_bulk_data_control <- dataAugm(GSE8056_MAS5_control)
augmented_bulk_data_early_stage <- dataAugm(GSE8056_MAS5_early_stage)
augmented_bulk_data_middle_stage <- dataAugm(GSE8056_MAS5_middle_stage)
augmented_bulk_data_late_stage <- dataAugm(GSE8056_MAS5_late_stage)
}