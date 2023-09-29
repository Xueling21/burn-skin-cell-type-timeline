

library(Matrix)
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
library("stringr")
library(cowplot)
library(harmony)
##########load the normal skin scRNA-Seq data 1############
GSM5176929 <- read.table("GSM5176929_counts_scHealthy1.txt",header = T,sep = '\t')
GSM5176930 <- read.table("GSM5176930_counts_scHealthy2.txt",header = T,sep = '\t')
GSM5176931 <- read.table("GSM5176931_counts_scHealthy3.txt",header = T,sep = '\t')
skin.data <- data.frame(GSM5176929,GSM5176930,GSM5176931)
#pancreas.data <- GSM5176930
platform1 <- data.frame(platform = rep("GSM5176929",819)) 
platform2 <- data.frame(platform = rep("GSM5176930",4813))
platform3 <- data.frame(platform = rep("GSM5176931",2590))
platform <- data.frame(rbind(platform1,platform2,platform3))

#create the Seurat object of the three normal skin samples
metadata <- platform
skin <- CreateSeuratObject(skin.data, meta.data = metadata)
skin@meta.data$platform <- metadata
mydata <- CreateSeuratObject(counts = mat, 
                             min.cells = 3,
                             min.features = 200)

## =============2.QC(quality control)
# the QC considers two parameters 
#1.the number of the unique features per cell. Here the unique feature is gene, which can be adjusted according to data quality.
#2. The gene number ratio of mitochondrion. Theoretically, the mitochondrial genome is only a small part compared to the nuclear genome.
#Damaged or dead cells often exhibit a large amount of mitochondrial contamination, so cells with a high proportion of mitochondrial gene expression are filtered.
grep(pattern = "^MT-", rownames(skin), value = T) #Search rownames (skin) for strings containing "MT-".
#Use the PercentageFeatureSet function to calculate the mitochondrial QC finger, percent.mt the mitochondrial gene ratio
skin[["percent.mt"]] <- PercentageFeatureSet(skin, pattern = "^MT-")

#  nFeature_RNA  Represents the number of genes measured per cell，
#  nCount       Represents the sum of the expression of all genes measured per cell，
#  percent.mt  Represents the proportion of mitochondrial genes measured
head(skin@meta.data, 5)
summary(skin@meta.data$nCount_RNA)
VlnPlot(skin, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

data <- skin@meta.data
library(ggplot2)
p <- ggplot(data = data, aes(x=nFeature_RNA)) + geom_density()
plot1 <- FeatureScatter(skin, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(skin, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p +plot1 + plot2

#Remove cells with too high a proportion of mitochondrial gene expression, and some extreme cells.
skin <- subset(skin, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 4000 & percent.mt < 5)

## =============3.standardization
# LogNormalize
skin <- NormalizeData(skin, normalization.method = "LogNormalize", scale.factor = 10000)
normalized.data <- skin[["RNA"]]@data
normalized.data[1:20,1:3]
dim(normalized.data)

## =============4.Identification of highly variable genes
#Identification of genes with high cell-to-cell expression (feature selection)
#The purpose of this step is to identify genes with very different amounts of expression from cell to cell for subsequent identification of cell types.
#We use the default parameter, the "VST" method, to select 2000 highly variable genes.
skin <- FindVariableFeatures(skin, selection.method = "vst", nfeatures = 2000)
dim(skin)

# Extract the top 10 highly variable genes
top10 <- head(VariableFeatures(skin), 10)
top10

# Demonstrate highly variable genes
plot1 <- VariableFeaturePlot(skin)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


#4.Cell classification
#Before classification, the data set must first be dimensionally reduced
## =============5.Scaling the data
skin <- ScaleData(skin)
scale.data <- skin[["RNA"]]@scale.data
dim(scale.data)
scale.data[1:10,1:4]

all.genes <- rownames(skin)
skin <- ScaleData(skin, features = all.genes)


## =============6 Dimension reduction
skin <- RunPCA(skin, features = VariableFeatures(object = skin))

# Visualization
VizDimLoadings(skin, dims = 1:2, reduction = "pca")
DimPlot(skin, reduction = "pca")
DimHeatmap(skin, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(skin, dims = 1:15, cells = 500, balanced = TRUE)
#2)Define the "dimensions" of the dataset
#Here we need to select the number of principal components for subsequent cell classification.
#The "dimension" defined here does not represent the number of cell types，
#It is a parameter that is used when classifying cells.

## =============7.Determine the number of PCs used
#JackStraw and Elbowcan determine the "dimension" of the data. But Elbow is more intuitive, and we chose Elbow results for interpretation.
#It can be seen that between 7 and 10 principal component (PC), the standard deviation of the data is basically not decreasing.
#So we need to choose between 7 and 10, and in order to respect the recommendations of the official website, we pick 10, that is, the first 10 principal components for the classification of cells.
skin <- JackStraw(skin, num.replicate = 100)
skin <- ScoreJackStraw(skin, dims = 1:20)
JackStrawPlot(skin, dims = 1:20)
ElbowPlot(skin)
#skin2 <- skin

## =============8.Cluster cells
# Firstly, a KNN plot based on Euclidean distance is constructed based on PCA space
skin <- FindNeighbors(skin, dims = 1:30)
skin <- FindClusters(skin, resolution = 1.5)
#Here we set dims = 1:30, that is, the first 10 principal components are selected to classify cells.
#The results of the classification are as follows, and it can be seen that the cells are divided into 9 categories.
head(Idents(skin), 5)
saveRDS(object = skin,file = "data/skin_cluster.rds") # for later singleR annotation usage.

## =============9.Visualize cells in low-dimensional space (UMAP/tSNE)
skin_umap <- RunUMAP(skin, dims = 1:24,n.neighbors = 5)
skin_tsne <- RunTSNE(skin, dims = 1:24,seed.use = 1)

# 可视化
DimPlot(skin, reduction = "umap", label = T, label.size = 5)
DimPlot(skin, reduction = "tsne", label = T, label.size = 5)
DimPlot(object = skin, reduction = "tsne", pt.size = .1, group.by ="platform" , label = TRUE, repel = TRUE)

#saveRDS(skin, file = "data/skin_tutorial.rds") 
#Visualization 2


   #Extract marker genes for each cell type
   #Using the FindMarkers command, differentially expressed genes of each cell type against other cell types can be found as biomarkers for that cell type.
   #The ident.1 parameter sets the cell types to be analyzed, and the min.pct indicates the proportion of the number of genes to the total number of cells of this cell type
   
   #find all markers of cluster 1
   cluster1.markers <- FindMarkers(skin, ident.1 = 1, min.pct = 0.25)
   head(cluster1.markers, n = 5)
   #The DoHeatmap command allows visualization of marker gene expression
   skin.markers <- FindAllMarkers(skin, only.pos = TRUE, min.pct = 0.25)
   #?FindMarkers
   top3 <- skin.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
   DoHeatmap(skin, features = top3$gene) + NoLegend()
   
   #6.Explore genes of interest
   #Seurat provieded many methods to easily explore the expression of genes of interest in various cell types
   #Seurat offers a number of methods that allow us to easily explore the expression of genes of interest in various cell types
   VlnPlot(skin, features = c( "BANK1", "IGHD","CD19","CD79A","MS4A1","Igha"))
   #monocytes
   VlnPlot(skin, features = c("TET2", "ASGR2", "APOBEC3A","CFP","ASGR1","ATP6V1B2"))
   #Neutrophil
   VlnPlot(skin, features = c("CREB5", "CDA", "TREM1","S100A12","APOBEC3A","CASP5"))
   #Neutrophil#################
   VlnPlot(skin, features = c("CSF3R", "PILRA", "CD10","MMP9","MDSC","Ly6G","CD11b"))
   
   VlnPlot(skin, features = c("CCL5", "OSM", "GZMB","FGFBP2","CD79A","MS4A1","Igha"))
   
   VlnPlot(skin, features = c("CCL5", "CLIC3", "GZMB","FGFBP2","CD79A","MS4A1","Igha"))
   
   VlnPlot(skin, features = c("ARL4C","C15orf48","CD83","LMNA","MMP9","PIM3","PXDC1","STK17B"))
   #CD8
   VlnPlot(skin, features = c("CD8","CD3","CD8A","CD8B","CD39","CX3CR1","CXCR5","CXCR6" ,"GNLY","GZMB","GZMK","ITGA1","ITGAE","NKG2D","NKG7","PRF1"))
   #CD4
   VlnPlot(skin, features = c("CD4","CD3","ABCG2","AFF3","ANXA7","BCL2","CCDC66","CD11C","CD3E","CD69","CCR7","CD28","SELL"))
   #角质
   VlnPlot(skin, features = c("KRT1","KRT5","KRT14","CK14","KRT15","ALDH","ATP1B1","CD44"))
   # marc    
   VlnPlot(skin, features = c("CD68", "C1QA", "CSF1R", "FCGR3A", "CD163","CD206","ANXA1","CD204","DAB2","MRC1","TYMP"))
   #mono
   VlnPlot(skin, features = c("CD14", "S100A9", "FCN1"))
   #mast
   VlnPlot(skin, features = c("TPSAB1", "CPA3"))
   VlnPlot(skin, features = c("CD3","TCRgd", "NCAM1", "NCR1"))
   VlnPlot(skin, features = c("CD4", "CD45RA"))
   #DC
   VlnPlot(skin, features = c("CD1C","CD2C","CD11C","CD1A","CD14","CD209","CD86","ITGAX","MHCII"))
   #cytotoxic
   VlnPlot(skin, features = c("CD355","CD8","granzymeB","GZMA","GZMB","IFN-γ","MAPKAPK5"))
   # NK
   VlnPlot(skin, features = c("CD56","GZMB","CD4","NCAM1","NCR1","NKG7","CCL5","CD3D"))
   #melanocyte
   VlnPlot(skin, features = c("MLANA","PMEL","DCT","MITF","HMB45","APOD"))
   #Epithelial_cells
   VlnPlot(skin, features = c( "CD31","VWF","PECAM1","CLDN5","CDH5","SERPINE1"))
    # stem
   VlnPlot(skin, features = c(  "NANOG","OCT4",
                                "MCSP","CD200","CD29","CD34","CD46","CD49f",
                                "CD133","SOX2",
                                "ACTA2","COL1A1","COL1A2","MYL9","TAGLN",
                                "CD13","SSEA-4","TRA-1-60","Cytokeratin-15","LGR5","LGR6","LRIG1","ABCB5","Beta-catenin"))
   # fibroblates
   VlnPlot(skin, features = c("COL1A1","DCN","APOD","COL1A2","FN1","LUM","CD90","CFD"))
   #neutrophil
   VlnPlot(skin, features = c("FCGR3B","CD10","CD16","CD16b","CD66B","CMTM2","VNN2"))
   #neurons
   VlnPlot(skin, features = c("PGP9.5","TUJ1"))
   VlnPlot(skin, features = c("NANOG","OCT4","SOX2"))
   
 
      #you can plot raw counts as well
   VlnPlot(skin, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
   
   FeaturePlot(skin_umap, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
   #This display method maps the amount of gene expression to the UMAP results, and the specificity of gene expression can also be intuitively seen.
   FeaturePlot(skin,reduction = "tsne", features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
   #This display method maps the amount of gene expression to the tSNE results, and the specificity of gene expression can also be intuitively seen.
     
   
   ####################7.Define cell types using prior knowledge
   #By comparing the marker gene we identified with the purposeful gene expression marker of the published cell type，
   #We can define the cell groups we divide. Finally, add names to our defined cell groups


# ref <- celldex::NovershternHematopoieticData()
# class(ref)
# table(ref$label.main)
# 
# hpca.se <- celldex::HumanPrimaryCellAtlasData()
# class(hpca.se)
# table(hpca.se$label.main)
# table(refdata$label.main)
# unique(hpca.se@colData@listData[["label.main"]])


   ## =============10.Cell annotation using SingleR
   library(Seurat)
   library(SingleR)
   refdata <- ref_Human_all
   #refdata <- ref_Monaco
   
   head(colnames(refdata))
   head(rownames(refdata))
   
   unique(refdata@colData@listData[["label.main"]])
   
   skin <- readRDS("data/skin_cluster.rds") # the input of singleR is the output of Seurat FindClusters() function.
   
   testdata <- GetAssayData(skin, slot="data")
   
   dim(testdata)
   testdata[1:30,1:4]
   
   clusters <- skin@meta.data$seurat_clusters
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
   
   #Plot the clusters
   new.cluster.ids3 <- celltype[,2]
   names(new.cluster.ids3) <- levels(skin)
   skin <- RenameIdents(skin, new.cluster.ids3)
   
   DimPlot(skin, reduction = "umap", label = TRUE, pt.size = 1.2)

   
   skin2 <- subset(skin, subset = seurat_clusters != 11 & seurat_clusters != 8 & seurat_clusters != 5&skin$seurat_clusters != 15&skin$seurat_clusters != 21&seurat_clusters != 22&seurat_clusters != 23&seurat_clusters != 18 )
   new.cluster.ids <- c("CD4_T", "Melanocyte","Endothelial", "Fibroblasts", "Macrophage", "CD8_T", "Keratinocytes")
   # skin4 <- skin2   
   names(new.cluster.ids) <- levels(skin2)
   skin4 <- RenameIdents(skin4, new.cluster.ids)
   DimPlot(skin, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
   DimPlot(skin4, reduction = "tsne", label = TRUE, pt.size =2) + NoLegend()
   
   my_cols <- c("CD4_T"="#FED976","Melanocyte"="#4DAF4A","Endothelial"= "#F781BF","Fibroblasts"="#377EB8","Macrophage"="#984EA3", "CD8_T"="#A65628", "Keratinocytes"="#FF7F00")
   DimPlot(skin4, reduction = "tsne", label = F,label.size = 3, pt.size =1.5,cols = my_cols)+
      theme(axis.text = element_text(size = rel(1.5),color="black"), ## 设置标签
            strip.text  = element_text(size = rel(3)),
            # plot.title  = element_text(size = rel(1.7),hjust=0.5,vjust=0.5),
            #panel.grid.minor = element_line(size = rel(3)),              ##
            axis.title.x=element_text(size=rel(1.5),hjust=0.5,vjust=0.5), ##Change the axis name of the x-axis
            axis.title.y=element_text(size=rel(1.5),hjust=0.5,vjust=0.5), ##Change the axis name of the y-axis
            axis.line =element_line(size=0.01,lineend='butt'), ##Make changes to the tick marks on the x-axis
            axis.ticks.length =  unit(0.5,"lines"),
            panel.border = element_rect(color="black",fill = NA,size=2.5),
            legend.title=element_blank() ,
            legend.key.size = unit(1, "line"),  #legend.key.width与legend.key.height Integration of functions
            legend.text = element_text(size = rel(1)), #Label font size
            legend.key.width = unit(0.5,"cm"), #
            legend.key.height = unit(1,"cm") 
            #legend.position = 'top'#control legend position in the plot up, down, left, right（'top','bottom','right','left'）
      )+coord_fixed() + guides(color=guide_legend(override.aes = list(size=3)))#,
   
   
   
ident <- data.frame(skin2@active.ident) 
ident_cluster <- data.frame(skin22@active.ident) 
ident_2 <- subset(ident_cluster,ident_cluster$skin22.active.ident != 2)
ident_2 <- subset(ident_2,ident_cluster$skin22.active.ident != 3)
ident_2 <- subset(ident_2,ident_cluster$skin22.active.ident != 7)

GSE169147_cell <- data.frame(skin2@assays$RNA@data)
GSE169147_cell_2[1:4,1:4]
colnames(GSE169147_cell) <- ident[match(colnames(GSE169147_cell),rownames(ident)),1]
colnames(GSE169147_cell) <- ident[match(rownames(ident_2),rownames(ident)),1]

### transform GSE169147 counts to tpm 
GSE169147 <- GSE169147_tpm 
Fibroblasts <- sample(GSE169147[,grep("Fibroblasts",colnames(GSE169147))], 200 ) 
Melanoma <- sample(GSE169147[,grep("Melanoma",colnames(GSE169147))], 200 ) 
Macrophage <- sample(GSE169147[,grep("Macrophage",colnames(GSE169147))], 200 ) 
Keratinocytes <- sample(GSE169147[,grep("Keratinocytes",colnames(GSE169147))], 150 ) 
CD4_T <- sample(GSE169147[,grep("CD4_T",colnames(GSE169147))], 200 ) 
CD8_T <- sample(GSE169147[,grep("CD8_T",colnames(GSE169147))], 200 ) 
Endothelial <- sample(GSE169147[,grep("Endothelial",colnames(GSE169147))], 200 ) 
Adipocyte <- sample(tpm,200)

CELL <- c(rep("Fibroblasts",200), rep("Melanoma",200),rep("Macrophage",200),rep("Keratinocytes",150),rep("CD4_T",200),rep("CD8_T",200),rep("Endothelial",200))

GSE169147_cell_200 <- data.frame(GeneSymbol = rownames(GSE169147_cell_300),Fibroblasts,Melanoma,Macrophage,Keratinocytes,CD4_T,CD8_T,Endothelial)
colnames(GSE169147_cell_300)[-1] <- CELL

write.table(GSE169147_cell_200,"GSE169147_cell_200.txt",sep ="\t", row.names =F, col.names =TRUE, quote =FALSE)



