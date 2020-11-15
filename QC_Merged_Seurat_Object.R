#setwd("path/to/wd")
library(Seurat)
library(cowplot)
library(dplyr)


#Script to create Seurat objects and perform QC.  After performing QC on
#each sample separately, I'll merge the objects and create 1 Seurat object.

#To identify low-quality cells I will assess the total counts, genes expressed
#and percent mitochondrial reads per cell in each sample
#The filtering will be different for each individual sample, each sample will have a different 
#distribution of counts, genes expressed and percent mitochondrial reads.

#Example control object
CTL_1<-Read10X(data.dir="Sample1_mm10_v3/outs/filtered_feature_bc_matrix")
Ctrl_1<-CreateSeuratObject(counts = CTL_1, project = "Ctrl_1", min.cells = 5)
Ctrl_1$Cre<-"CTRL"
Ctrl_1$Sample<-"Sample_Name"
Ctrl_1$Age<-"Adult"
Ctrl_1$Chemistry<-"V3"
Ctrl_1[["percent.mt"]]<-PercentageFeatureSet(Ctrl_1, pattern="^mt-")
#Before filtering: 5643 cells with 16882 genes
VlnPlot(Ctrl_1, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
VlnPlot(Ctrl_1, features="percent.mt", y.max=20) #get better look at pt. mito reads
Ctrl_1<-subset(Ctrl_1, subset=nFeature_RNA > 500  & nFeature_RNA <6000 & percent.mt<10)
Ctrl_1
#After filtering: 5180 cells with 16882 genes

#Example p60-70 induced tumor cells
GFAP_A1<-Read10X(data.dir="Sample2_mm120_v3/outs/filtered_feature_bc_matrix")
Gfap_A1<-CreateSeuratObject(counts=GFAP_A1, project="Gfap_A1", min.cells = 5)
Gfap_A1$Cre<-"GFAP"
Gfap_A1$Sample<-"Sample_Name"
Gfap_A1$Age<-"Adult"
Gfap_A1$Chemistry<-"V3"
Gfap_A1[["percent.mt"]]<-PercentageFeatureSet(Gfap_A1, pattern="^mt-")
Gfap_A1
#Before filtering: 6608 cells with 18177 genes
VlnPlot(Gfap_A1, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
VlnPlot(Gfap_A1, features=c("percent.mt"), y.max=15)
Gfap_A1<-subset(Gfap_A1, subset=nFeature_RNA > 500 & nFeature_RNA <7500 & percent.mt<8)
Gfap_A1
#5230 cells with 18177 genes
 
#Example p18 sample
NES_1<-Read10X(data.dir="Sample1/filtered_gene_bc_matrices/mm10")
Nes_1<-CreateSeuratObject(counts=NES_1, project ="Nes_1", min.cells = 5)
Nes_1$Cre<-"NES"
Nes_1$Sample<-"Sample_Name"
Nes_1$Age<-"Pup"
Nes_1$Chemistry<-"V2"
Nes_1[["percent.mt"]]<-PercentageFeatureSet(object=Nes_1, pattern="^mt-")
VlnPlot(Nes_1, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
Nes_1<-subset(Nes_1, subset=nFeature_RNA > 500 & nFeature_RNA <7500 & percent.mt<5)
Nes_1
#2211 cells, 17048 genes

#Repeat for each each sample subsetting cells based on number of features and percent 
#mitochondrial reads.   Filters will be different for each sample.


#Merge all Seurat objects:
All.combined<-merge(Ctrl_1, y = c(Ctrl_2, Ctrl_3, Gfap_A1, Gfap_A2, Gfap_A3, Gfap_A4,
                                    Nes_1, Nes_2, Nes_3, Nes_4,
                                    Gfap_Ad1, Gfap_Ad2, Gfap_Ad3, Gfap_Ad4), #Seurat objects
                                    add.cell.ids = c("Sample1", "Sample2", "Sample3",
                                                     "Sample4", "Sample5", "sample6",
                                                     "Sample7", "Sample8", "Sample9",
                                                     "Sample10", "Sample11", "sample12"), #Sample names
                                   project = "Combined")

#check to make sure all samples are present in the object
table(All.combined$orig.ident)

#Save merged Seurat object
save(All.combined, file="AllSamples.Rdata") 
