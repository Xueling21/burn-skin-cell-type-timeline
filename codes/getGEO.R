
######## For GEO microarray downloading and preprocessing (a demo)
setwd("H:\\burnSCRNA\\test_data")
#affy  .CEL format  GSE59732 ##################################################
library(affy)
library(stringr)
library(readxl)
library(openxlsx)
affy_data1=ReadAffy(celfile.path = "GSE103744_RAW")
## standardization
library(dplyr)
#probes_expr1=rma(affy_data1) %>% exprs() %>% as.data.frame() # with log2 format
probes_expr1=mas5(affy_data1) %>% exprs() %>% as.data.frame() # without log2 format for CIBERSORT, need mas5 preprocessing
gpl570 = read.delim("GSE103744_family.xml/GPL10558-tbl-1.txt", header=FALSE, stringsAsFactors=FALSE)
gpl570=gpl570[!grepl("///",gpl570$V11),c(1,11)]
colnames(gpl570)=c("ID_REF","SYMBOL")

colnames(probes_expr1) <- str_split(colnames(probes_expr1),"_", simplify = TRUE)[,1]

head(gpl570)
library(XML)
library(plyr)
MyContact=ldply(xmlToList("GSE103744_family.xml\\GSE103744_family.xml"), function(x)data.frame(x,stringsAsFactors = F))
gse_pheno=MyContact[1:(nrow(MyContact)-2),c(7,seq(14,32,2))] #phenodata

colnames(gse_pheno)=c("sample",MyContact[1,seq(15,33,2)])

gse_pheno=as.data.frame(apply(gse_pheno,2,function(x)gsub(" |\n","",x)))
gse_pheno=unique(gse_pheno)
rownames(gse_pheno)=gse_pheno$sample

probes_expr1$ID_REF=rownames(probes_expr1)
probes_expr1=merge(gpl570,probes_expr1,by="ID_REF")
rowMan=apply(probes_expr1[,3:ncol(probes_expr1)],1,function(x)mean(as.numeric(x),na.rm=T))
probes_expr1 =probes_expr1[order(rowMan,decreasing =T),]
probes_expr1 =probes_expr1[!duplicated(probes_expr1$SYMBOL),]
rownames(probes_expr1)=probes_expr1$SYMBOL
probes_expr1=probes_expr1[-1,3:ncol(probes_expr1)]
probes_expr1[1:3,1:3]
colnames(probes_expr1)<-rownames(gse_pheno)
library(readxl)
library(openxlsx)
write.table(probes_expr1,"GSE103744_mas5_neut.txt",sep ="\t", row.names =F, col.names =TRUE, quote =FALSE)


#####geoChina('GSE97368') downloadable####
#options(stringsAsFactors = F)#When using as.data.frame，set stringsAsFactors to be FALSE to avoid character automatically changing into factor type.
# load wide usable packages
library(GEOmirror)
library(idmap1)
library(idmap2) 
library(idmap3)
# remotes::install_github("jmzeng1314/idmap3")
# check the content
ls('package:GEOmirror')
ls('package:idmap1')
ls('package:idmap2')
ls('package:idmap3')

# download data
geoChina('GSE103744')
load("GSE103744_eSet.Rdata")
gset
# extract data information 
pdata=pData(gset[[1]]) 


# extract expression matrix 
a=exprs(gset[[1]])
dim(a) 
# annotation and probe information
a[1:4,1:4]
gset[[1]]@annotation
#发现GPL17077在idmap2包中
ids <- getIDs("GPL17077")#失败
ids <- idmap2::get_soft_IDs('GPL17077')#成功
head(ids)
genes_expr <- filterEM(probes_expr,ids )
write.table(genes_expr,"GSE97368_genes_expr.txt",sep ="\t", row.names =F, col.names =TRUE, quote =FALSE)


#####geoChina('GSE182616')不可以下载数据的####
#options(stringsAsFactors = F)#在调用as.data.frame的时，将stringsAsFactors设置为FALSE可以避免character类型自动转化为factor类型
#加载万能包
library(AnnoProbe)
library(GEOmirror)
library(GEOquery)
library(idmap1)
library(idmap2) 
library(idmap3)
setwd("H:\\GSE182616")
# remotes::install_github("jmzeng1314/idmap3")

#了解一下每个包
ls('package:GEOmirror')
ls('package:idmap1')  #第一个万能芯片探针ID注释平台R包
ls('package:idmap2')  #第二个万能芯片探针ID注释平台R包
ls('package:idmap3')  #第三个万能芯片探针ID注释平台R包
#下载数据
gset<-getGEO(filename ="GSE182616_series_matrix.txt.gz", destdir=".", AnnotGPL = T)#加载下载的文件
# 提取表型数据信息 
pdata=pData(gset)
write.xlsx(pdata,"GSE182616_pdata.xlsx",rowNames = T)

# 提取表达矩阵 0
eSet=gset 
probes_expr <- exprs(eSet)
dim(probes_expr) 
boxplot(genes_expr[1:30,])

# 芯片平台的设计注释信息
GPL <- gset@annotation
ids <- idmap2::get_soft_IDs(GPL)#一列是探针信息二列为基因名
#将表达矩阵中的标签更换成相应的基因名
genes_expr <- filterEM(probes_expr,ids )
write.table(genes_expr,"GSE182616_genes_expr.txt",sep ="\t", row.names =F, col.names =TRUE, quote =FALSE)


#############Agilent  数据处理####################

##设置临床路径，这里我已经下载好了GSE182616的原始数据GSE182616_RAW
setwd("F:\\GSE182616")
##加载R包
library(AnnoProbe)
library(GEOmirror)
library(GEOquery)
library(idmap1)
library(idmap2) 
library(idmap3)
library(Biobase)
library(limma)
library(readxl)
library(openxlsx)

#提供一个GSE编号就可以下载啦。因为表达矩阵是处理过的，
#我们不要，所以只提取临床信息表格，从中获得分组信息。
gset<-getGEO(filename ="GSE182616_series_matrix.txt.gz", destdir=".", AnnotGPL = T)#加载下载的文件
# 提取表型数据信息 
pdata=pData(gset)
write.xlsx(pdata,"GSE182616_pdata.xlsx",rowNames = T)


##获取原始数据

file1 = list.files("GSE182616_RAW3", pattern = "txt.gz",full.names = T)
dat_aver=read.maimages(file=file,source = "agilent",green.only=FALSE) ##单色芯片

dat_aver1<-dat_aver
##前期处理
####背景校正和标准化处理####
dat_aver=backgroundCorrect(dat_aver,"normexp",normexp.method = "rma",offset = 50)
dat_aver=normalizeBetweenArrays(dat_aver,method = "quantile")
dat_aver=avereps(dat_aver,ID=dat_aver$genes$ProbeName)
exprset=as.data.frame(dat_aver$E)
probes_expr <- exprset

# 芯片平台的设计注释信息
GPL <- gset@annotation
ids <- idmap2::get_soft_IDs(GPL)#一列是探针信息二列为基因名
#将表达矩阵中的标签更换成相应的基因名
genes_expr <- filterEM(probes_expr,ids )
write.table(genes_expr,"GSE182616_genes_expr.txt",sep ="\t", row.names =F, col.names =TRUE, quote =FALSE)


