
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
#options(stringsAsFactors = F)#When using as.data.frame��set stringsAsFactors to be FALSE to avoid character automatically changing into factor type.
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
#����GPL17077��idmap2����
ids <- getIDs("GPL17077")#ʧ��
ids <- idmap2::get_soft_IDs('GPL17077')#�ɹ�
head(ids)
genes_expr <- filterEM(probes_expr,ids )
write.table(genes_expr,"GSE97368_genes_expr.txt",sep ="\t", row.names =F, col.names =TRUE, quote =FALSE)


#####geoChina('GSE182616')�������������ݵ�####
#options(stringsAsFactors = F)#�ڵ���as.data.frame��ʱ����stringsAsFactors����ΪFALSE���Ա���character�����Զ�ת��Ϊfactor����
#�������ܰ�
library(AnnoProbe)
library(GEOmirror)
library(GEOquery)
library(idmap1)
library(idmap2) 
library(idmap3)
setwd("H:\\GSE182616")
# remotes::install_github("jmzeng1314/idmap3")

#�˽�һ��ÿ����
ls('package:GEOmirror')
ls('package:idmap1')  #��һ������оƬ̽��IDע��ƽ̨R��
ls('package:idmap2')  #�ڶ�������оƬ̽��IDע��ƽ̨R��
ls('package:idmap3')  #����������оƬ̽��IDע��ƽ̨R��
#��������
gset<-getGEO(filename ="GSE182616_series_matrix.txt.gz", destdir=".", AnnotGPL = T)#�������ص��ļ�
# ��ȡ����������Ϣ 
pdata=pData(gset)
write.xlsx(pdata,"GSE182616_pdata.xlsx",rowNames = T)

# ��ȡ������� 0
eSet=gset 
probes_expr <- exprs(eSet)
dim(probes_expr) 
boxplot(genes_expr[1:30,])

# оƬƽ̨�����ע����Ϣ
GPL <- gset@annotation
ids <- idmap2::get_soft_IDs(GPL)#һ����̽����Ϣ����Ϊ������
#����������еı�ǩ��������Ӧ�Ļ�����
genes_expr <- filterEM(probes_expr,ids )
write.table(genes_expr,"GSE182616_genes_expr.txt",sep ="\t", row.names =F, col.names =TRUE, quote =FALSE)


#############Agilent  ���ݴ���####################

##�����ٴ�·�����������Ѿ����غ���GSE182616��ԭʼ����GSE182616_RAW
setwd("F:\\GSE182616")
##����R��
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

#�ṩһ��GSE��žͿ�������������Ϊ��������Ǵ������ģ�
#���ǲ�Ҫ������ֻ��ȡ�ٴ���Ϣ���񣬴��л�÷�����Ϣ��
gset<-getGEO(filename ="GSE182616_series_matrix.txt.gz", destdir=".", AnnotGPL = T)#�������ص��ļ�
# ��ȡ����������Ϣ 
pdata=pData(gset)
write.xlsx(pdata,"GSE182616_pdata.xlsx",rowNames = T)


##��ȡԭʼ����

file1 = list.files("GSE182616_RAW3", pattern = "txt.gz",full.names = T)
dat_aver=read.maimages(file=file,source = "agilent",green.only=FALSE) ##��ɫоƬ

dat_aver1<-dat_aver
##ǰ�ڴ���
####����У���ͱ�׼������####
dat_aver=backgroundCorrect(dat_aver,"normexp",normexp.method = "rma",offset = 50)
dat_aver=normalizeBetweenArrays(dat_aver,method = "quantile")
dat_aver=avereps(dat_aver,ID=dat_aver$genes$ProbeName)
exprset=as.data.frame(dat_aver$E)
probes_expr <- exprset

# оƬƽ̨�����ע����Ϣ
GPL <- gset@annotation
ids <- idmap2::get_soft_IDs(GPL)#һ����̽����Ϣ����Ϊ������
#����������еı�ǩ��������Ӧ�Ļ�����
genes_expr <- filterEM(probes_expr,ids )
write.table(genes_expr,"GSE182616_genes_expr.txt",sep ="\t", row.names =F, col.names =TRUE, quote =FALSE)

