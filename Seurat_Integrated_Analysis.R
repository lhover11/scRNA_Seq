setwd("path/to/wd")
library(Seurat)
library(cowplot)
library(dplyr)
library(viridis)
####################SEURAT INTEGRATION#################################
#Because 1/3 of my samples were processed using 10X 3' Single Cell v3 chemistry and 2/3 were processed using 10x v2 chemistry
#I used Seurat integration to find common cell types.  Without integration, cells cluster by 10x chemistry


#Performed QC, filtered all samples and created object in script: CreateSeuratObject.R

load("AllSamples.Rdata")  #Seurat object of all 15 samples merged (tumor and control samples)

#split object by 10x version, here termed chemistry
in.list <- SplitObject(AllSamples, split.by = "chemistry")

#log normalize and find variable features:
in.list <- lapply(X = in.list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

#Identify anchors between our datasets:
anchors <- FindIntegrationAnchors(object.list = in.list, dims = 1:20)

#Use anchors to integrate datasets
combined <- IntegrateData(anchorset = anchors, dims = 1:20)

DefaultAssay(combined) <- "integrated"

#Standard Seurat workflow for visualization and clustering
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)

DimPlot(combined, label=T)  #25 clusters
DimPlot(combined, reduction = "umap", group.by = "Age")
DimPlot(combined, reduction = "umap", group.by = "Cre")

#Extract the number of cells per cluster by our Cre driver:
table<-table(Idents(combined), combined$Cre)
#write.csv(table, "cellspercluster.csv")

#To find markers and look at differences in gene expression, switch to the RNA assay
DefaultAssay(combined) <- "RNA"
combined<- NormalizeData(combined, verbose = TRUE)
combined <- ScaleData(combined, verbose=TRUE)
AllClusters<-FindAllMarkers(combined, only.pos = TRUE, logfc.threshold = 0.5)
write.csv(AllClusters, file="AllClusterMarkers.csv")

#Store averages of all clusters:
Combined.Avg.Exp<-AverageExpression(combined, return.seurat = FALSE,  slot = "data", use.scale = FALSE, verbose = TRUE)

#After looking through our cluster markers and our proportions of cells/cluster we can identify which clusters are 
#comprised of normal cells and which clusters are comprised of tumor cells

##Select tumor cells and recluster those:
Tumor<-subset(combined, idents=c("0", "1", "2", "3", "4", "7", "10", "12", "13", "17"))

#Remove any remaining control cells
Tumor<-subset(Tumor, subset = orig.ident =="Ctrl2991", invert=TRUE)
Tumor<-subset(Tumor, subset = orig.ident =="Ctrl2940", invert=TRUE)
Tumor<-subset(Tumor, subset = orig.ident =="Ctrl3148", invert=TRUE)

Tumor
#21775 features across 37559 samples within 2 assays #Active assay: RNA (21557 features)
#1 other assay present: integrated
#2 dimensional reductions calculated: pca, umap
DefaultAssay(Tumor) <- "integrated"

# Run the standard workflow for visualization and clustering
Tumor <- ScaleData(Tumor, verbose = TRUE)
Tumor <- RunPCA(Tumor, npcs = 30, verbose = TRUE)
# t-SNE and Clustering
Tumor <- FindNeighbors(Tumor, reduction = "pca", dims = 1:30, verbose=TRUE)
Tumor <- FindClusters(Tumor, resolution = 0.2, verbose=TRUE)  #9 clusters
#Tumor <- FindClusters(Tumor, resolution = 0.4) #Look at clusters at other resolutions
#Tumor <- FindClusters(Tumor, resolution = 0.3)

Tumor <- RunUMAP(Tumor, reduction = "pca", dims = 1:20)
Tumor <- RunTSNE(Tumor, reduction = "pca", check_duplicates = FALSE, dims = 1:20)

# Visualization
DimPlot(Tumor, label=T)
DimPlot(Tumor, split.by = "Age", label=T)

table(Idents(Tumor), Tumor$Age)

#Change colors of clusters:
library(RColorBrewer)

#Colors from Rcolor brewer paired palette:
DimPlot(Tumor, cols=c("#a6cee3", "#fb9a99", "#B2DF8A", "#33A02C", "#1f78b4", "#ff7f00", "#e31a1c", "#6a3d9a", "#ff7f00"), label=F)

#Find cluster markers:
DefaultAssay(Tumor) <- "RNA"
Tumor<- NormalizeData(Tumor, verbose = TRUE)
Tumor <- ScaleData(Tumor, verbose=TRUE)
AllClusters<-FindAllMarkers(Tumor, only.pos = TRUE, logfc.threshold = 0.5)
write.csv(AllClusters, file="TumorMarkers.csv")

#Cells per cluster for each individual tumor:
samples<-table(Idents(Tumor), Tumor$orig.ident)

#Save tumor cells Seurat object for continued analysis
