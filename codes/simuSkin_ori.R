################
# Prepred by Xueling Li and Min Zhu
# Contact Xueling Li via email: xlli@cmpt.ac.cn or xuelingli16@foxmail.com

############ Get the average means of 15 simulations for one cell type
# e.g. csGEPs_B03_ct<- csGEPs_B03[[i]]$Neutrophil
# e.g. stderr_B03_ct <- stderr_B03[[i]]$Neutrophi
# ix, the element index of a R list transformed from group deconvolution GEP data (JobX_GEPs_Filtered.TXT files)...
# ,,,, downloaded from CIBERSORTx webserver. 
# csGEPs_B03_ct, the cell type specific gene expression profiles from JobX_GEPs_Filtered.TXT files. 
getAverageMean <- function(csGEPs_B03_ct,ix){
    b03 <- data.frame(matrix(data=NA,nrow=21654,ncol=15))
    for(i in 1:15){
    b03[,i] <- csGEPs_B03_ct[[i]][ix]
    }
  gname <- row.names(csGEPs_B03_ct[[i]][ix])
  row.names(b03)<-gname
  b03[b03==1]<-NA
  mean_b03 <- rowMeans(b03, na.rm = TRUE)
  return(mean_b03) 
}

## Averaged standard errors of 15 simulations for one cell type
# ix, the element index of a list of group deconvolution data  (_JobX_GEPs_Filtered.TXT files)...
# ,,,, downloaded from CIBERSORTx webserver. 
# csGEPs_B03_ct, the cell type specific gene expression profiles from _JobX_GEPs_StdErrs.TXT files. 
getAverageErr <- function(stderr_B03_ct,ix){
  b03_err <- data.frame(matrix(data=NA,nrow=21654,ncol=15))
  for(i in 1:15){
     b03_err[,i] <- stderr_B03_ct[[i]][ix]
  }
  gname <- row.names(stderr_B03_ct[[i]][ix])
  row.names(b03_err)<-gname
  
  mean_bvar <- rowMeans(b03_err,na.rm=TRUE)
  return(mean_bvar) 
}


############ Load  one group csGEPs of four groups, i.e., control, 0_3,4_7, and 8_17 post-burn skin 
# fpath: the absolute or relateve path to the _JobX_GEPs_StdErrs.TXT files
loadcsGEPs <- function(fpath,fctrl){
  
  len <- length(fctrl)
  csGEPs_c <- list()
  
  for (i in 1:len){
    cuPath <- paste0(fpath,fctrl[i])
    cuFile <- list.files(cuPath)
    # fgepC <- list.files(cuPath,pattern = "_GEPs")[1]
    fgepC <- list.files(cuPath,pattern = "_GEPs_Filtered")[1]
    csGEPs_c[[i]] <- read.table(paste0(cuPath,"\\",fgepC),header=T,sep="\t",check.names=F,row.names=1)
  }
 
  avemean_b03_allCT <- sapply(1:ncol(csGEPs_c[[i]]), function(ix) (getAverageMean(csGEPs_c,ix)))
  return(avemean_b03_allCT)
}

##### the same as loadcsGEPs function except modifying control by ignoring CD4 and CD8 extra cell types
############load  one group of csGEPs of four groups, i.e., control, 0_3,4_7, and 8_17 post burns of skin 
loadcsGEPs2 <- function(fpath,fctrl){
  len <- length(fctrl)
  csGEPs_c <- list()
  for (i in c(1:3,5:15)){
    cuPath <- paste0(fpath,fctrl[i])
    cuFile <- list.files(cuPath)
    # fgepC <- list.files(cuPath,pattern = "_GEPs")[1]
    fgepC <- list.files(cuPath,pattern = "_GEPs_Filtered")[1]
    csGEPs_c[[i]] <- read.table(paste0(cuPath,"\\",fgepC),header=T,sep="\t",check.names=F,row.names=1)
  }
  avemean_b03_allCT <- sapply(1:ncol(csGEPs_c[[i]]), function(ix) (getAverageMean(csGEPs_c,ix)))
  
  return(avemean_b03_allCT)
}



###########load StdErr of each group from _JobX_GEPs_StdErrs.TXT files
loadStdErr <- function(fpath,fctrl){
  len <- length(fctrl)
  stderr_c <-list()
  for (i in 1:len){
    cuPath <- paste0(fpath,fctrl[i])
    cuFile <- list.files(cuPath)
    ferrC <- list.files(cuPath,pattern = "_StdErrs")
    stderr_c[[i]] <- read.table(paste0(cuPath,"\\",ferrC),header=T,sep="\t",check.names=F,row.names=1)
  }
  aveerr_b03_allCT <- sapply(1:ncol(stderr_c[[i]]), function(ix) (getAverageErr(stderr_c,ix)))
  return(aveerr_b03_allCT)
  
}

####### same as loadStdErr function except modifying by ignoring the 4th file with CD4+ and CD8+ T cells
###########load StdErr of one group
loadStdErr2 <- function(fpath,fctrl){
  len <- length(fctrl)
  stderr_c <-list()
  for (i in c(1:3,5:15)){
    cuPath <- paste0(fpath,fctrl[i])
    cuFile <- list.files(cuPath)
    ferrC <- list.files(cuPath,pattern = "_StdErrs")
    stderr_c[[i]] <- read.table(paste0(cuPath,"\\",ferrC),header=T,sep="\t",check.names=F,row.names=1)
  }
  aveerr_b03_allCT <- sapply(1:ncol(stderr_c[[i]]), function(ix) (getAverageErr(stderr_c,ix)))
  return(aveerr_b03_allCT)
  
}

# ############## load Q value files
# loadQ <- function(fpath,fctrl){
#   len <- length(fctrl)
#   Q_c <-list()
#   for (i in 1:len){
#     cuPath <- paste0(fpath,fctrl[i])
#     cuFile <- list.files(cuPath)
#     ferrC <- list.files(cuPath,pattern = "_Qvals")
#     Q_c[[i]] <- read.table(paste0(cuPath,"\\",ferrC),header=T,sep="\t",check.names=F,row.names=1)
#   }
#   return(Q_c)
#   
# }

####### modify by ignoring the 4th file with CD4+ and CD8+ T cells
###########load Qval of one group
# loadQ2 <- function(fpath,fctrl){
#   len <- length(fctrl)
#   stderr_c <-list()
#   for (i in c(1:3,5:15)){
#     cuPath <- paste0(fpath,fctrl[i])
#     cuFile <- list.files(cuPath)
#     ferrC <- list.files(cuPath,pattern = "_Qvals")
#     Q_c[[i]] <- read.table(paste0(cuPath,"\\",ferrC),header=T,sep="\t",check.names=F,row.names=1)
#   }
#   return( Q_c)
# }
# 


#### apply the averaging cell type mean and std functions to all list
##############load four groups,i.e.,control, 0_3,4_7, and 8_17 post burns of skin  
load_allcsGEPs <-function(){
  # setwd("D:\\工作\\导师\\安医\\FXY\\新建文件夹\\GSE8056_group")
  # fpath <- "D:\\工作\\导师\\安医\\FXY\\新建文件夹\\GSE8056_group\\"
  setwd("Path_To_GSE8056_files\\GSE8056_group")
  fpath <- "Path_To_GSE8056_files\\GSE8056_group\\"
  
  fctrl <- list.files(fpath, pattern = "_C_") 
  fB03 <- list.files(fpath, pattern = "_B_0_3")
  fB47 <- list.files(fpath, pattern = "_B_4_7")
  fB817 <- list.files(fpath, pattern = "_B_8_17")
  
  avemeanAllCT_c<-loadcsGEPs2(fpath,fctrl)
  avemeanAllCT_b03<-loadcsGEPs(fpath,fB03)
  avemeanAllCT_b47<-loadcsGEPs(fpath,fB47)
  avemeanAllCT_b817<-loadcsGEPs(fpath,fB817)

  aveMean_4Grp <-list(c=avemeanAllCT_c,b03=avemeanAllCT_b03,b47=avemeanAllCT_b47,b817=avemeanAllCT_b817)
  return(aveMean_4Grp)
}

################################ separate the error from mean as standalone function
load_allcsSTD <-function(){
  # setwd("D:\\工作\\导师\\安医\\FXY\\新建文件夹\\GSE8056_group")
  # fpath <- "D:\\工作\\导师\\安医\\FXY\\新建文件夹\\GSE8056_group\\", the path to GSE8056 group
  
  setwd("Path_To_GSE8056_files\\GSE8056_group")
  fpath <- "Path_To_GSE8056_files\\GSE8056_group\\"
  
  fctrl <- list.files(fpath, pattern = "_C_") 
  fB03 <- list.files(fpath, pattern = "_B_0_3")
  fB47 <- list.files(fpath, pattern = "_B_4_7")
  fB817 <- list.files(fpath, pattern = "_B_8_17")

  avestderr_c<-loadStdErr2(fpath,fctrl)
  avestderr_b03<-loadStdErr(fpath,fB03)
  avestderr_b47<-loadStdErr(fpath,fB47)
  avestderr_b817<-loadStdErr(fpath,fB817)
  
  aveStderr_4Grp <-list(c=avestderr_c,b03=avestderr_b03,b47=avestderr_b47,b817<-avestderr_b817)
  return(aveStderr_4Grp)
}


# load_allQ <-function(){
#   setwd("D:\\工作\\导师\\安医\\FXY\\新建文件夹\\GSE8056_group")
#   fpath <- "D:\\工作\\导师\\安医\\FXY\\新建文件夹\\GSE8056_group\\"
#   
#   fctrl <- list.files(fpath, pattern = "_C_") 
#   fB03 <- list.files(fpath, pattern = "_B_0_3")
#   fB47 <- list.files(fpath, pattern = "_B_4_7")
#   fB817 <- list.files(fpath, pattern = "_B_8_17")
#   
#   Q_c<-loadQ2(fpath,fctrl)
#   Q_b03<-loadQ(fpath,fB03)
#   Q_b47<-loadQ(fpath,fB47)
#   Q_b817<-loadSQ(fpath,fB817)
#   
#   Q_4Grp <-list(c=Q_c,b03=Q_b03,b47=Q_b47,b817<-Q_b817)
#   return(aveStderr_4Grp)
# }




############ single contrast#################
## m1, the mean expression profiles of the control of one of the seven cell types;
## m2, the mean expression profiles of early, middle or late stage of one of the seven cell types;
## sd1, the standard deviations of the control of one of the seven cell types;
## sd2, the standard deviations of early, middle or late stage of one of the seven cell types;
### thres, threshold for q-values, i.e., the corrected p-values by BH methods

grpTest2 <- function(m1,m2,sd1,sd2,thres){
  
  q_z <-list()
   for (i in 1:ncol(m1)){  
    
    vBetaZ <- sapply(1:nrow(m1), function(j) (m1[j,i]-m2[j,i])/sqrt(sd1[j,i]^2+sd2[j,i]^2))
    ZPS <- 2*pnorm(-abs(vBetaZ))
    qvals <- p.adjust(ZPS, method = "BH")
   
    zval =vBetaZ[which(qvals < thres)]
     len1 <- length(zval)
     df <-data.frame(matrix(data=NA,nrow = len1))
    
    df$geneSymbol =row.names(m1)[which( qvals < thres)]
    df$Zval <-zval
    df$qval = qvals[which(qvals < thres)]
    df$mean_b = m1[,i][which(qvals < thres)]
    df$mean_c =m2[,i][which(qvals < thres)]
   
    df<-df[,-1]
    q_z[[i]]=df
    
  }
  
  names(q_z) <- c("Neutrophils","Fibroblasts","Melanocyte","Macrophage", "Keratinocytes", "T_cell","Endothelial")
  return(q_z)
  
}


##############      export result
# df, output of grpTest2
# file1, specified excel file names to save the csDEGs
exportGroupResult <-function(df,file1){
  library("openxlsx")
  library(dplyr)
 
  data_names <- c("Neutrophils","Fibroblasts","Melanocyte","Macrophage", "Keratinocytes", "T_cell","Endothelial")
  work_book1 <- createWorkbook()
                
  for(i in 1:7) {
    addWorksheet(work_book1,sheetName = data_names[i])
    # df <- bind_rows(q_z[[i]])
    # df <- df[order(df[[i]]$Zval),]
    writeData(work_book1, data_names[i],df[[i]], rowNames = F, colNames = T )
    saveWorkbook(work_book1, file=paste0(file1,".xlsx"),overwrite = TRUE)
  }
}


##########final group contrast #############
## aveMean_4Grp[[j]], means of seven cell types in each post-burn stage
#aveMean_4Grp[[1]], means of seven cell types in control
#aveStderr_4Grp[[j]], standard deviations of seven cell types in each post-burn stage
#aveStderr_4Grp[[1]], standard deviations of seven cell types in  Control
 groupContrast <-function(){
   aveMean_4Grp <-load_allcsGEPs()
   aveStderr_4Grp<-load_allcsSTD()
   
   
   # qz_all <- lapply(2:4,function(j) grpTest2(aveMean_4Grp[[j]],aveMean_4Grp[[1]],aveStderr_4Grp[[j]],aveStderr_4Grp[[1]],0.05))
   #fileA <- c("b03_vs_c","b47_vs_c","b817_vs_c")
   # lapply(1:3,function(j) exportGroupResult(qz_all[[j]],fileA[j]))
   qz_all_filtered <- lapply(2:4,function(j) grpTest2(aveMean_4Grp[[j]],aveMean_4Grp[[1]],aveStderr_4Grp[[j]],aveStderr_4Grp[[1]],0.05))
   fileA <- c("b03_vs_c_filtered","b47_vs_c_filtered","b817_vs_c_filtered")
   setwd("D:\\工作\\导师\\安医\\FXY\\新建文件夹")
   lapply(1:3,function(j) exportGroupResult(qz_all_filtered[[j]],fileA[j]))
   save(qz_all,"qz_all.RData")
   return(qz_all_filtered)
   
  }
 
 ##############
 enrichAnal <- function(csDEGs){
   library(openxlsx)# read .xlsx file
   library(ggplot2)# bar and dot plots
   library(stringr)# gene ID mapping
   library(enrichplot)#GO,KEGG,GSEA
   library(clusterProfiler)#GO,KEGG,GSEA
   library(GOplot)# for chord diagram，chord table plot，clustering plot
   library(DOSE)
   library(ggnewscale)
   library(topGO)# plot the pathway network diagram
   library(circlize)# Plot the enrichment analysis circle
   library(ComplexHeatmap)# Draw a legend
   library("openxlsx")
   # ################################ (1)GO enrichment analysis
  # csDEGs <- DEGs$geneSymbol[DEGs$Zval<0] ## down- or up=regulated genes
   csgene <- bitr(csDEGs,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = "org.Hs.eg.db")
   #csgene_EN <- bitr(csDEGs,fromType = 'SYMBOL',toType = 'ENSEMBLE',OrgDb = "org.Hs.eg.db")
   go <- enrichGO(csgene$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2, qvalueCutoff = 1,keyType = 'ENTREZID')
   close
   # }
   # write.csv(go,file=paste0("go_neutrophil_csDEGs_",outName,".csv")) ##
  return(go)
 }  
 
################# gene mapping from ENTRZID to SYMBOL
gmap <- function(go){
    len <- dim(go)[1]
    if(len!=0){
  entrz <- lapply(1:len, function(i)  str_split(go$geneID[[i]],"/"))
  entrz2sym <- lapply(1:len, function(i) bitr(entrz[[i]][[1]],fromType ='ENTREZID',toType = 'SYMBOL' ,OrgDb = "org.Hs.eg.db"))
  gSym <- vector(mode = "character")
  for(j in 1:len){
    len2 <- length(as.vector(entrz2sym[[j]]$SYMBOL))
    str1 <- entrz2sym[[j]]$SYMBOL[1]
    for(i in 2:len2) {str1 <-paste(str1,entrz2sym[[j]]$SYMBOL[i],sep="/")}
    gSym[j]<-str1
  }
  go$symbol <- gSym
  go <-go[,-c(2,9)]
    }else{
  go <-go[,-2] 
    }
    
  return(go) 
}

############# visualization of enriched GO terms
visalGO <- function(go){
  ######  (2) visualize the GO enrichment
  p1_go <- ggplot(head(go@result,20),aes(x=GeneRatio,y=Description))+
    geom_point(aes(size = `Count`, color = `pvalue`), shape = 16)+
    scale_size(range = c(1, 10))+scale_color_gradient(low = "#cf0000",high = "#f5c6c6")
  plot(p1_go)
  p1_go + theme(text=element_text(size=20), #change font size of all text
                axis.text=element_text(size=20), #change font size of axis text
                # axis.text.x=element_text(size=20,angle = 45, vjust = 1, hjust=1)),
                # axis.text.y=element_text(size=20,angle = 45, vjust = 1, hjust=1)),
                axis.title=element_text(size=20), #change font size of axis titles
                plot.title=element_text(size=20), #change font size of plot title
                legend.text=element_text(size=20), #change font size of legend text
                legend.title=element_text(size=20)) #change font size of legend title
}
 
################# batch processing for non-filtered GO terms
  batchGOenr1 <-function(qz_all){
   outName <- c("Neutrophils","Fibroblasts","Melanocyte","Macrophage", "Keratinocytes", "T_cell","Endothelial")
   tpb <- c("03","47","817")
   go_all_ori <-list()
   for (i in 1:3){
     go_all_ori[[i]]<-lapply(1:7,function(j) enrichAnal(qz_all[[i]][[j]]))
     names(go_all_ori[[i]])<-outName
    }
   names(go_all_ori)<-tpb
   return(go_all_ori) 
 }
  
################ save GO data frame for filtered GO terms of all cell-type specific differentially expressed genes (csDEGs)
  batchGOenr2 <-function(go_all_ori){
   outName <- c("Neutrophils","Fibroblasts","Melanocyte","Macrophage", "Keratinocytes", "T_cell","Endothelial")
   tpb <- c("03_f","47_f","817_f")
   tpb_g <- c("03_g0_filteredg_all","47_go_filteredg_all","817_go_filteredg_all")
   go_df <- list()
   for (i in 1:3){
     go_df[[i]]<-lapply(1:7,function(j) gmap(go_all_ori[[i]][[j]]))
     names(go_df[[i]])<-outName
   }
   names(go_df)<-tpb
   lapply(1:3,function(i) exportGroupResult(go_df[[i]],tpb_g[i]))
   return(go_df) 
  }
 
  ########  for down regulated csDEGs 
 batchGOenr_down <-function(qz_all){
    outName <- c("Neutrophils","Fibroblasts","Melanocyte","Macrophage", "Keratinocytes", "T_cell","Endothelial")
    tpb <- c("03","47","817")
   
    tpb_dn <- c("03_dn","47_dn","817_dn")
    go_df_dn <- list()
    for (i in 1:3){
      go_all_ori <-lapply(1:7,function(j) enrichAnal(qz_all[[i]][[j]]$geneSymbol[qz_all[[i]][[j]]$Zval<0]))
      go_df_dn[[i]]<-lapply(1:7,function(j) gmap(go_all_ori[[j]]@result))
    }
    names(go_df_dn)<-tpb
   lapply(1:3,function(i) exportGroupResult(go_df_dn[[i]],tpb_dn[i]))
    return(go_df_dn) 
  }
 ########  for up regulated csDEGs  
 batchGOenr_up <-function(qz_all){
   outName <- c("Neutrophils","Fibroblasts","Melanocyte","Macrophage", "Keratinocytes", "T_cell","Endothelial")
   tpb <- c("03","47","817")
   tpb_up <- c("03_up","47_up","817_up")
   go_df_up <- list()
   for (i in 1:3){
     go_all_ori <-lapply(1:7,function(j) enrichAnal(qz_all[[i]][[j]]$geneSymbol[qz_all[[i]][[j]]$Zval>0]))
     go_df_up[[i]]<-lapply(1:7,function(j) gmap(go_all_ori[[j]]@result))
   }
   names(go_df_up)<-tpb
   lapply(1:3,function(i) exportGroupResult(go_df_up[[i]],tpb_up[i]))
   return(go_df_up) 
 }
 
 
