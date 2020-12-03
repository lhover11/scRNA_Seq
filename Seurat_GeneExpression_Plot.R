setwd("path/to/wd")
library(Seurat)
library(dplyr)
load("SeuratObject.Rdata")  #Here my object is named "Tumor"

DefaultAssay(Tumor)<-"RNA"
# Normalize and scale RNA data 
Tumor<-NormalizeData(Tumor, verbose = TRUE)
Tumor<-ScaleData(Tumor, verbose=TRUE)

Idents(Tumor)<-"integrated_snn_res.0.2"  #Setting identities as the cluster labels

#Arrange meta data in the order you wnat for the figures
Tumor$Age<-factor(Tumor$Age, levels=c("Pup", "Adolescent", "Adult"))

#Get out a list of all cells that have expression of the gene of interest (here expression>0)
GfapExp<-WhichCells(Tumor, expression = Gfap>0)

#Show on umap all cells that have expression of Gfap, show these cells in green
DimPlot(Tumor, cells.highlight=GfapExp, split.by="Age", cols.highlight = c("forestgreen"), sizes.highlight =1)

#Now split by age of tumors
DimPlot(Tumor, cells.highlight=GfapExp, split.by="Sample", cols.highlight = c("forestgreen"), sizes.highlight = .5)

#Display expression in violin plot:
VlnPlot(Tumor, features=c("Gfap"), y.max=1)
