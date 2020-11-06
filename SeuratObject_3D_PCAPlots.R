setwd("path/to/your/wd")
library(Seurat)
library(cowplot)
library(dplyr)
library(viridis)
library(datapasta)
library(plotly)


##Here we create a 3D PCA to get a better feel for the separation of the single cell clusters using plotly

#Load in your seurat object, at this point cells have already been assigned to clusters
#load("SeuratObject.Rdata")

#First look at 2D PCA plots to see which populations/gene signatures are separated on PC1-PC4 

#Here my Seurat object is named tumor and I want to set the cell identities to the cluster assignments under "integrated_snn_res.0.2"
Idents(Tumor)<-"integrated_snn_res.0.2"  
DimPlot(Tumor, label=T)
DimPlot(Tumor, dims=1:2, reduction="pca")
DimPlot(Tumor, dims=2:3, reduction="pca")
DimPlot(Tumor, dims=3:4, reduction="pca")

#PC1 separates cycling cells from the rest of the tumor cells

#Now I will plot 3D PCAs without cycling cells to look at how the remaining cells clusters.  I want to plot each age group independently 
#to see how the tumors differ at each time point (Pup, Adolescent and Adult)

#1st subset by age:
Idents(Tumor)<-"Age"
Pup<-subset(Tumor, idents=c("Pup"))

Idents(Pup)<-"integrated_snn_res.0.2"
Pup<-subset(Pup, idents=c("0", "1", "4", "6", "7"))  #non-cycling cells
Pup  #6182 cells
DimPlot(Pup, dims=1:2, reduction="pca")
DimPlot(Pup, dims=2:3, reduction="pca", cols=c("#a6cee3", "#fb9a99", "#1f78b4", "#e31a1c", "#6a3d9a"), pt.size=1)

#3D plot
# Get out PC1,2 and 3 values with cluster ID
PCs_Pup <- FetchData(Pup, vars = c("PC_1", "PC_2", "PC_3", "integrated_snn_res.0.2"))
table(PCs_Pup$integrated_snn_res.0.2)

#Set colors to match those of the UMAP plot, clusters not being shown assigned to "white" (2,3,8)
set1<-c("#a6cee3", "#fb9a99", "white", "white", "#1f78b4", "white", "#e31a1c", "#6a3d9a", "white")

Pup_Plot<-plot_ly(data = PCs_Pup, 
        x = ~PC_1, y = ~PC_2, z = ~PC_3, 
        color = PCs_Pup$integrated_snn_res.0.2, colors = set1,
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 3, width=1))

Pup_Plot <- Pup_Plot %>% layout(scene = list(
  xaxis = list(range = c(-20, 30)),
  yaxis = list(range = c(-30, 30)),
  zaxis=list(range=c(-10,35))))

Pup_Plot 


Idents(Tumor)<-"Age"
Adult<-subset(Tumor, idents=c("Adult"))
Idents(Adult)<-"orig.ident"
#Adult<-subset(Adult, idents=c("Sample4"), invert=TRUE)  #Remove only sample with large scale CNAs

Idents(Adult)<-"integrated_snn_res.0.2"
Adult<-subset(Adult, idents=c("0", "1", "2", "3", "4", "6", "7"))
Adult #10332 cells
DimPlot(Adult, dims=2:3, reduction="pca") 

#3D Plot
# Get out PC1,2 and 3 with cluster ID
PCs_Adult <- FetchData(Adult, vars = c("PC_1", "PC_2", "PC_3", "integrated_snn_res.0.2"))

#Same color scheme
set1<-c("#a6cee3", "#fb9a99", "blue", "blue", "#1f78b4", "white", "#e31a1c", "#6a3d9a", "white")

Adult_plot<-plot_ly(data = PCs_Adult, 
           x = ~PC_1, y = ~PC_2, z = ~PC_3, 
           color = PCs_Adult$integrated_snn_res.0.2, colors = set1,
           type = "scatter3d",
           mode = "markers",
           marker = list(size = 3, width=1))

Adult_plot <- Adult_plot %>% layout(scene = list(
  xaxis = list(range = c(-20, 30)),
  yaxis = list(range = c(-30, 30)),
  zaxis=list(range=c(-10,30))))


Adult_plot


Idents(Tumor)<-"Age"
Adoles<-subset(Tumor, idents=c("Adolescent"))

Idents(Adoles)<-"integrated_snn_res.0.2"
Adoles<-subset(Adoles, idents=c("0", "1", "4", "6", "7"))
Adoles #7718
DimPlot(Adoles, dims=2:3, reduction="pca", cols=c("#a6cee3", "#fb9a99", "#1f78b4", "#e31a1c", "#6a3d9a"), pt.size=1)

#3D plot:
# Get out PC1,2 and 3 with cluster ID
PCs_Adol <- FetchData(Adoles, vars = c("PC_1", "PC_2", "PC_3", "integrated_snn_res.0.2"))

Adol_plot<-plot_ly(data = PCs_Adol, 
           x = ~PC_1, y = ~PC_2, z = ~PC_3, 
           color = PCs_Adol$integrated_snn_res.0.2, colors = set1,
           type = "scatter3d",
           mode = "markers",
           marker = list(size = 3, width=1))

Adol_plot <- Adol_plot %>% layout(scene = list(
  xaxis = list(range = c(-20, 30)),
  yaxis = list(range = c(-30, 30)),
  zaxis=list(range=c(-10,30))))


Adol_plot

