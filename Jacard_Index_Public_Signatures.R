setwd("path/to/wd")
library(ggplot2)
library(tidyverse)

#Script to calculate the Jaccard Index to determine the similarity between
#my tumor signatures and 435 published single cell and bulk RNAseq signatures.

#Jaccard Index/Similarity: the intersect of A+B / ((A+B)-length of intersect)

#Table with all published cell signatures (435), each row is a signature, 
#column 1 signature name
#column 2 signature description
#column >=3 genes in the signature
Table<-read.table("Master_SingleCell_List_for_Jaccard.txt",
                  sep="\t", header=FALSE, stringsAsFactors = FALSE)  
Table[1:5,1:5]
row.names(Table)<-Table$V1
Table<-Table[,-c(1,2)]
#Each sign name as rownames, only genes in the signatures in the table

#Read in tumor signatures: logFC>0.5 for each tumor cluster versus all other cells
#(normal and tumor combined)

Sign<-read.delim("Tumor_Signatures.txt", sep="\t", stringsAsFactors = FALSE)
Sign<-Sign[,1:8]


#Jaccard Index:
dat<-data.frame(matrix(NA, ncol = 8, nrow = 435))  #start with "empty" dataframe,
#you will add the Jaccard index to each cell,
#start with each cell containing NA

colnames(dat)<-colnames(Sign)  #set column names of new df to my tumor signatures
row.names(dat)<-row.names(Table) #set rownames to the published signatures

#For each column assignment below run the for loop, run each column through the for loop
#separately (run col<-1 through loop then repeat with col<-2)
col<-1
col<-2
col<-3
col<-4
col<-5
col<-6
col<-7
col<-8

#For loop to run through each row of Table (each row = a published signature (435))
#and calculate the jaccard similarity\with each column of Sign (my Tumor signatures)
for(row in seq_len(nrow(Table))) {
  genes <- Table[row,] #genes for published signature
  genes <- unique(genes[genes != ""]) #remove empty cells
  My_list<-Sign[,col]  #genes for each tumor sign
  My_list <-unique(My_list[My_list != ""]) #remove empty cells
  x<- 0+row #row index (1-435)
  #Calculate Jaccard index:
  I<-length(intersect(My_list,genes))  
  S <- I/((length(My_list)+length(genes))-I)
  dat[x, col]<-S  #in my new df, for each specified row & column enter the jaccard similarity index
}

#Repeat for loop until you've calculated similarity for each tumor signature

write.csv(dat, "JaccardFinalList.csv")

########################################
#Visualization of results:
Jaccard<-dat
Jaccard$X<-row.names(Jaccard)

#Select top 10 signatures with highest jaccard score for each of my clusters

#First gather the data into 3 columns: pub. set, my tumor sign and the jaccard score
Gather<-gather(Jaccard, key = "Cluster", value = "jaccard",-X)

#Select top 10 highest jaccard scores for each of my tumor signatures
Gather2<-Gather %>% group_by(Cluster) %>% slice_max(jaccard, n=10, with_ties=FALSE)

#Set order of signatures so we can order the published datasets for plotting purposes
Gather2$Cluster<-factor(Gather2$Cluster, levels=c("S.Phase", "G2M.Phase", "Progenitor", "OPC", "Oligo",
                                                     "Mes", "Astro", "Immune"))

#Arrange so jaccard indices are in descending order for each tumor sign.
Gather2<-Gather2 %>% group_by(Cluster) %>% arrange(desc(jaccard), .by_group = TRUE)

#Get the order of the published signatures:
Order<-unique(Gather2$X)

#Now go back to df that includes the scores for each signature and order it according to the
#the order of published sets we just specified 
Data<-Jaccard[match(Order, Jaccard$X),]
row.names(Data)<-Data$X

##Dot plot:

#Keep order info
Data$Order<-1:nrow(Data)
Data2<-gather(Data, key = "Cluster", value = "jaccard", -X, -Order)

#Order of my tumor clusters for the x axis
Data2$Cluster<-factor(Data2$Cluster, levels=c("S.Phase", "G2M.Phase", "Progenitor", "OPC", "Oligo",
                                              "Mes", "Astro", "Immune"))

ggplot(Data2, aes(x=Cluster, y=reorder(X, -Order))) + 
  geom_point(aes(size=jaccard, color=jaccard)) +  
  scale_size(range=c(-.5, 4))+
  theme_classic(base_size=14) +
  scale_color_gradient(low ="#0dad8d", high="#1164b4")+  #blue-red
  ylab("Single Cell Data Sets") +
  ggtitle("Jaccard Similarity Score")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.y = element_text(size=6))





