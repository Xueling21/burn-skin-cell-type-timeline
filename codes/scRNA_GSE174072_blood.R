

library(Matrix)
setwd("H:\\burnSCRNA\\GSE174072")
load("GSE124880.RData")
library(dplyr)
library(Seurat)
library(patchwork)
library(Rtsne) 
library(ggplot2)
library(tidyr)
library(cowplot)
library("FactoMineR")
library("factoextra")
library("ROCR")
##########Read the data1############
## =============1.Load the dataset

#data <- readRDS("H:\\burnSCRNA\\GSM5285694_ACK_cell.counts.matrices.rds")
#data_exon <- data.matrix(data$exon) 
#data_HIP015 <- readRDS("H:\\burnSCRNA\\GSM5285708_HIP015_ACK_cell.counts.matrices.rds")
#data_HIP015_exon <- data.matrix(data_HIP015$exon) 
data_HIP002 <- readRDS("H:\\burnSCRNA\\GSM5285706_HIP002_ACK_cell.counts.matrices.rds") # the path to data is different
data_HIP002_exon <- data.matrix(data_HIP002$exon)

blood <- CreateSeuratObject(counts = data_HIP002_exon, 
                           min.cells = 3,
                           min.features = 200)
## =============2.QC quality control analysis
#Two parameters for QC： 
#1.Number of unique features detected per cell (unique features represent the number of genes detected by a cell and can be adjusted according to the quality of the data)
#2.The proportion of mitochondrial genes detected in each cell, theoretically accounts for only a small part of the mitochondrial genome compared to the nuclear genome.
#When damaged or dead cells often exhibit a large amount of mitochondrial contamination, cells with too high a proportion of mitochondrial gene expression are filtered.
grep(pattern = "^MT-", rownames(blood), value = T) #Search rownames (blood) for strings containing "MT-"
#Use the PercentageFeatureSet function to calculate mitochondrial QC index, percent.mt the mitochondrial gene ratio
blood[["percent.mt"]] <- PercentageFeatureSet(blood, pattern = "^MT-")

#  nFeature_RNA  represents the number of genes measured per cell，
#  nCount       represents the sum of the expression of all genes measured in each cell，
#  percent.mt  represents the proportion of mitochondrial genes measured.
head(blood@meta.data, 5)
summary(blood@meta.data$nCount_RNA)

VlnPlot(blood, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

data <- blood@meta.data
library(ggplot2)
p <- ggplot(data = data, aes(x=nFeature_RNA)) + geom_density()
plot1 <- FeatureScatter(blood, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(blood, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p +plot1 + plot2

#Remove cells with too high a proportion of mitochondrial gene expression, and some extreme cells.
blood <- subset(blood, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & nCount_RNA < 5000 & percent.mt < 5)
blood


## =============3.Standardization
# LogNormalize
blood <- NormalizeData(blood, normalization.method = "LogNormalize", scale.factor = 10000)

normalized.data <- blood[["RNA"]]@data
normalized.data[1:20,1:3]
dim(normalized.data)

## =============4.Identification of highly variable genes
#Identification of genes with high cell-to-cell expression variation (feature selection)（feature selection）
#The purpose of this step is to identify genes with very different amounts of expression from cell to cell for subsequent identification of cell types.
#We use the default parameter, the "VST" method, to select 2000 highly variable genes.
blood <- FindVariableFeatures(blood, selection.method = "vst", nfeatures = 2000)
dim(blood)

# Extract the top 10 highly variable genes
top10 <- head(VariableFeatures(blood), 10)
top10

# Demonstrate highly variable genes
plot1 <- VariableFeaturePlot(blood)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


#4.Cell classification
#Before classification, the data set must first be dimensionally reduced
## =============5.Scaling the data
blood <- ScaleData(blood)
scale.data <- blood[["RNA"]]@scale.data
dim(scale.data)
scale.data[1:10,1:4]

all.genes <- rownames(blood)
blood <- ScaleData(blood, features = all.genes)


## =============6.降维
blood <- RunPCA(blood, features = VariableFeatures(object = blood))

# 可视化
VizDimLoadings(blood, dims = 1:2, reduction = "pca")
DimPlot(blood, reduction = "pca")
DimHeatmap(blood, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(blood, dims = 1:15, cells = 500, balanced = TRUE)
#2)Define the "dimensions" of the dataset
#Here we need to select the number of principal components for subsequent cell classification.
#The "dimension" defined here does not represent the number of cell types.
#It is a parameter that is used when classifying cells.

## =============7.确定使用PC个数
#JackStraw和Elbow都可以决定数据的“维度”。但是Elbow比较直观，我们选择Elbow结果进行解读。
#可以看到，主成分（PC）7到10之间，数据的标准差基本不在下降.
#所以我们需要在7到10之间进行选择，为了尊重官网的建议，我们选取10，即前10个主成分用于细胞的分类。
blood <- JackStraw(blood, num.replicate = 100)
blood <- ScoreJackStraw(blood, dims = 1:20)
JackStrawPlot(blood, dims = 1:15)
ElbowPlot(blood)
#blood2 <- blood

## =============8.对细胞聚类
# 首先基于PCA空间构建一个基于欧氏距离的KNN图
blood <- FindNeighbors(blood, dims = 1:10)
blood <- FindClusters(blood, resolution = 1.6)
#这里我们设置了dims = 1:10 即选取前10个主成分来分类细胞。
#分类的结果如下,可以看到，细胞被分为9个类别。
head(Idents(blood), 5)


## =============9.将细胞在低维空间可视化UMAP/tSNE
blood <- RunUMAP(blood, dims = 1:10)
blood <- RunTSNE(blood, dims = 1:10)

# 可视化
DimPlot(blood, reduction = "umap", label = T, label.size = 5)
DimPlot(blood, reduction = "tsne", label = T, label.size = 5)
#saveRDS(blood, file = "data/blood_tutorial.rds") 

#提取各个细胞类型的marker gene
#提取各个细胞类型的marker gene
#利用 FindMarkers 命令，可以找到找到各个细胞类型中与其他类别的差异表达基因，作为该细胞类型的生物学标记基因。
#其中ident.1参数设置待分析的细胞类别，min.pct表示该基因表达数目占该类细胞总数的比例

#find all markers of cluster 1
cluster1.markers <- FindMarkers(blood, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
#利用 DoHeatmap 命令可以可视化marker基因的表达
blood.markers <- FindAllMarkers(blood, only.pos = TRUE, min.pct = 0.25)
#?FindMarkers
top3 <- blood.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(blood, features = top3$gene) + NoLegend()

#6.探索感兴趣的基因
#Seurat提供了许多方法使我们能够方便的探索感兴趣的基因在各个细胞类型中的表达情况
#Seurat提供了许多方法使我们能够方便的探索感兴趣的基因在各个细胞类型中的表达情况
VlnPlot(blood, features = c( "BANK1", "IGHD","CD19","CD79A","MS4A1","Igha"))
#monocytes
VlnPlot(blood, features = c("TET2", "ASGR2", "APOBEC3A","CFP","ASGR1","ATP6V1B2"))
#Neutrophil
VlnPlot(blood, features = c("CREB5", "CDA", "TREM1","S100A12","APOBEC3A","CASP5"))
#Neutrophil#################
VlnPlot(blood, features = c("CSF3R", "PILRA", "CD10","MMP9","MDSC","Ly6G","CD11b"))
#CD8
VlnPlot(blood, features = c("CD8","CD3","CD8A","CD8B","CD39","CX3CR1","CXCR5","CXCR6" ,"GNLY","GZMB","GZMK","ITGA1","ITGAE","NKG2D","NKG7","PRF1"))
#CD4
VlnPlot(blood, features = c("CD4","CD3","ABCG2","AFF3","ANXA7","BCL2","CCDC66","CD11C","CD3E","CD69","CCR7","CD28","SELL"))
# NK
VlnPlot(blood, features = c("CD56","GZMB","CD4","NCAM1","NCR1","NKG7","CCL5","CD3D"))

#这种展示方法把基因表达量映射到tSNE结果中，同样可以直观的看到基因表达的特异性。


####################7.利用先验知识定义细胞类型
#通过对比我们鉴定的marker gene与已发表的细胞类型特意的基因表达marker，
#可以定义我们划分出来的细胞类群。最后，给我们定义好的细胞类群加上名称
new.cluster.ids <- c("Neutrophils", "Neutrophils", "Neutrophils", "CD4_T", "Neutrophils","Neutrophils",
                     "Neutrophils", "Neutrophils","CD8_T", "Neutrophils","Neutrophils", "Neutrophils",
                     "CD8_T","NK_cells", "NK_cells", " B_cells","Neutrophils", "CD4_T", "Monocytes"," B_cells")


names(new.cluster.ids) <- levels(blood2)
blood2 <- RenameIdents(blood2, new.cluster.ids)
DimPlot(blood, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

blood2 <- subset(blood, subset = seurat_clusters != 24 & seurat_clusters != 19 & seurat_clusters != 23&blood$seurat_clusters != 22&blood$seurat_clusters != 20)
DimPlot(blood, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


## =============10.使用SingleR进行细胞注释
library(Seurat)
library(SingleR)
refdata <- ref_Monaco
head(colnames(refdata))
head(rownames(refdata))

unique(refdata@colData@listData[["label.main"]])


blood <- readRDS("data/mydata_tutorial.rds")

testdata <- GetAssayData(blood, slot="data")

dim(testdata)
testdata[1:30,1:4]

clusters <- blood@meta.data$seurat_clusters
table(clusters)


cellpred <- SingleR(test = testdata,
                    ref = refdata,
                    labels = refdata$label.main,
                    clusters = clusters,
                    assay.type.test = "logcounts",
                    assay.type.ref = "logcounts")
#mca_result2 <- scMCA(scdata = dat, numbers_plot = 3)


str(cellpred,max.level = 3)
metadata <- cellpred@metadata
head(metadata)

celltype = data.frame(ClusterID = rownames(cellpred), 
                      celltype = cellpred$labels, 
                      stringsAsFactors = F)
celltype

p = plotScoreHeatmap(cellpred, clusters = rownames(cellpred), order.by = "cluster")

#绘制聚类图
new.cluster.ids3 <- celltype[,2]
names(new.cluster.ids3) <- levels(blood)
blood <- RenameIdents(blood, new.cluster.ids3)

DimPlot(blood, reduction = "umap", label = TRUE, pt.size = 1.2)
DimPlot(blood, reduction = "tsne", label = TRUE, pt.size = 1.2)


levels(blood2) 
my_cols2 <- c("Neutrophils"= "#CD3333", "CD4_T"="#FED976","CD8_T"="#A65628","NK_cells"="#ADD8E6"," B_cells"="#E6E6FA","Monocytes"="#C1FFC1" )
DimPlot(blood2, reduction = "tsne", label = F,label.size = 3, pt.size =1.5,cols = my_cols2)+
   theme(axis.text = element_text(size = rel(1.5),color="black"), ## 设置标签
         strip.text  = element_text(size = rel(3)),
         # plot.title  = element_text(size = rel(1.7),hjust=0.5,vjust=0.5),
         #panel.grid.minor = element_line(size = rel(3)),              ##
         axis.title.x=element_text(size=rel(1.5),hjust=0.5,vjust=0.5), ##对x轴的坐标轴名称进行改动
         axis.title.y=element_text(size=rel(1.5),hjust=0.5,vjust=0.5), ##对y轴的坐标轴名称进行改动
         axis.line =element_line(size=0.01,lineend='butt'), ##对x轴的刻度线进行改动
         axis.ticks.length =  unit(0.5,"lines"),
         panel.border = element_rect(color="black",fill = NA,size=2.5),
         legend.title=element_blank() ,
         legend.key.size = unit(1, "line"),  #legend.key.width与legend.key.height 功能的整合
         legend.text = element_text(size = rel(1)), #标签字体大小
         legend.key.width = unit(0.5,"cm"), #
         legend.key.height = unit(1,"cm") 
         #legend.position = 'top'#控制图例放置在图中上下左右（'top','bottom','right','left'）
         
   )+
   coord_fixed() + guides(color=guide_legend(override.aes = list(size=3)))#,
# shape = guide_legend(override.aes = list(size=3)))
#dev.off()


ident <- data.frame(blood2@active.ident) 
ident_cluster <- data.frame(blood@active.ident) 
ident_2 <- subset(ident_cluster,ident_cluster$blood2.active.ident != 2)
ident_2 <- subset(ident_2,ident_cluster$blood2.active.ident != 3)
ident_2 <- subset(ident_2,ident_cluster$blood2.active.ident != 7)


cell_HIP002 <- data.frame(blood2@assays$RNA@data)
ACK_HIP002_cell_2[1:4,1:4]
colnames(cell_HIP002) <- ident[match(colnames(cell_HIP002),rownames(ident)),1]
cell_nember <- match(rownames(ident_2),colnames(ACK_HIP002_cell_2))
cell_nember <- na.omit(cell_nember)
ACK_HIP002_cell_2 <- ACK_HIP002_cell_2[,cell_nember]
colnames(ACK_HIP002_cell_2) <- ident[match(colnames(ACK_HIP002_cell_2),rownames(ident)),1]

colnames(ACK_HIP002_cell) <- ident[match(rownames(ident_2),rownames(ident)),1]

table(colnames(cell_HIP002))
new.cluster.ids <- c("Neutrophils", "Neutrophils", "Neutrophils", "CD4_T", "Neutrophils","Neutrophils",
                     "Neutrophils", "Neutrophils","CD8_T", "Neutrophils","Neutrophils", "Neutrophils",
                     "CD8_T","NK_cells", "NK_cells", " B_cells","Neutrophils", "CD4_T", "Monocytes"," B_cells")


ACK_HIP002 <- cell_HIP002_tpm
B_cells <- sample(ACK_HIP002[,grep("B_cells",colnames(ACK_HIP002))], 400 ) 
CD4_T <- sample(ACK_HIP002[,grep("CD4_T",colnames(ACK_HIP002))], 600 ) 
Monocytes <- sample(ACK_HIP002[,grep("Monocytes",colnames(ACK_HIP002))], 300 ) 
Neutrophils <- sample(ACK_HIP002[,grep("Neutrophils",colnames(ACK_HIP002))], 1000 ) 
NK_cells <- sample(ACK_HIP002[,grep("NK_cells",colnames(ACK_HIP002))], 600 ) 
CD8_T <- sample(ACK_HIP002[,grep("CD8_T",colnames(ACK_HIP002))], 600 ) 
CELL <- c(rep("B_cells",400), rep("CD4_T",600),rep("CD8_T",600),rep("NK",600),rep("Monocytes",300),rep("Neutrophils",1000))

HIP002cell_200 <- data.frame(B_cells,CD4_T,CD8_T,NK_cells,Monocytes,Neutrophils)
HIP002cell_200 <- data.frame(GeneSymbol = rownames(HIP002cell_200),HIP002cell_200)
colnames(HIP002cell_200)[-1] <- CELL
write.table(HIP002cell_200, file ="HIP002cell_200.txt", sep ="\t", row.names =F, col.names =TRUE, quote =FALSE)

Neutrophils <- data.frame(gene = rownames(Neutrophils),Neutrophils)
GSE169147_cell_200 <- GSE169147_cell_200[,-1]
colnames(GSE169147_cell_200)[1] <- "gene"
GSE169147_cell_200 <- data.frame(gene = rownames(GSE169147_cell_200),GSE169147_cell_200)
GSE169147_Neu <- merge(Neutrophils,GSE169147_cell_200,by = "gene")
colnames(GSE169147_Neu)[1] <- "GeneSymbol"
write.table(GSE169147_Neu,"GSE169147_Neu_2.txt",sep ="\t", row.names =F, col.names =TRUE, quote =FALSE)


