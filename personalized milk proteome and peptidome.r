#library
library(readxl)
library(corrplot)
library(ggplot2)
library(ggfortify)
library(factoextra)
library(magrittr)
library(dplyr)
library(plyr)
library(reshape2)
library(ggpubr)
library(tidyverse)
library(Hmisc)
library(colorspace)
library(stringr)
library(data.table)
library(ComplexHeatmap)
library(dendextend)
library(devtools)
library(msm)
library(ggbiplot)
library(circlize)


# Proteome

#import data to analysis
setwd("D:/Projects to finish/Personalized milk/Data analysis/WM001002")
wm <- read_excel("WMproteome001and002_Proteins.xlsx",sheet = "ProteinAveAbundance")

#Normalized by sum
wm_percentage <- apply(wm[,c(18:71)], 2, function(x)(x*100/sum(x,na.rm = TRUE)))
wm_percentage <- apply(wm_percentage,2,function(x)(replace(x,is.na(x),0)))
wm_percentage_df <- data.frame(wm[,c(4:17)],wm_percentage)

#Calculated as mg/ml by Total amino acid
wm_mgml <- data.frame(wm_percentage[,c(1:3)]*19.18808/100,
                      wm_percentage[,c(4:6)]*19.66452/100,
                      wm_percentage[,c(7:9)]*15.179/100,
                      wm_percentage[,c(10:12)]*12.72446/100,
                      wm_percentage[,c(13:15)]*12.84298/100,
                      wm_percentage[,c(16:18)]*12.18363/100,
                      wm_percentage[,c(19:21)]*12.60597/100,
                      wm_percentage[,c(22:24)]*12.95751/100,
                      wm_percentage[,c(25:27)]*10.63312/100,
	     wm_percentage[,c(28:30)]*16.25398/100,
                      wm_percentage[,c(31:33)]*14.30421/100,
                      wm_percentage[,c(34:36)]*14.05865/100,
                      wm_percentage[,c(37:39)]*15.25564/100,
                      wm_percentage[,c(40:42)]*11.0365/100,
                      wm_percentage[,c(43:45)]*12.59312/100,
                      wm_percentage[,c(46:48)]*12.28856/100,
                      wm_percentage[,c(49:51)]*11.11747/100,
                      wm_percentage[,c(52:54)]*11.64257/100)
wm_mgml_df <- data.frame(wm[,c(4:17)],wm_mgml)

# Heatmap
# z score (normalized distribution transform)
protein_z <- apply(wm_mgml,1,function(x)(scale(as.numeric(x),center = TRUE,scale = TRUE)))
protein_z <- as.data.frame(t(protein_z))
colnames(protein_z) <- colnames(wm_mgml)
protein_z$Accession <- wm$Accession

protein_z_n <- protein_z[,c(1:54)]
min(protein_z_n) 
max(protein_z_n)

hrpro <- hclust(dist(protein_z_n,method="euclidean"), method="complete")
hcpro <- hclust(as.dist(1-cor(protein_z_n)), method="single")

# first Heatmap without cluster annoatation
pdf("protein_heatmap_wm001002.pdf",40,40)
Heatmap(protein_z_n,cluster_columns = as.dendrogram(hcpro),cluster_rows = as.dendrogram(hrpro),
        col = colorRamp2(c(-4,-2,0,2,4), 
                         c("#03294f","#2166ac","#f7f7f7","#c51b7d","#b50140")),
        heatmap_legend_param = list(at = c(-4,-2,0,2,4)))
dev.off()

# define cluster
cuttree <- cutree(hrpro, h = 13) #cut tree by height, alternatively use k = number to cut tree by defined number
tc_pro <- data.frame(clust_membership = as.character(cuttree))
str(tc_pro) #readout will show how many clusters are by the defined height

# export the cluter code for each imput protein
tc_pro$order <- c(1:length(hrpro$order))
proorder <- data.frame(order=hrpro$order)
proorder <- merge(proorder,tc_pro,by="order")
proteinorder <- cbind(proorder,protein_z,wm_mgml)
write.csv(proteinorder,"proteinorder_cluster5.csv")

# Heatmap with cluster annoatation
pdf("protein_heatmap_wm001002_cluster_5.pdf",20,20)
Heatmap(protein_z_n,cluster_columns = as.dendrogram(hcpro),cluster_rows = as.dendrogram(hrpro), 
        col = colorRamp2(c(-4,-2,0,2,4), 
                         c("#03294f","#2166ac","#f7f7f7","#c51b7d","#b50140")),
        heatmap_legend_param = list(at = c(-4,-2,0,2,4)))+
  rowAnnotation(tc_pro,annotation_legend_param =tc_pro$clust_membership)
dev.off()       

#overlayed lineplots by cluster membership (by each protein and emphasized median)
zprotein<-as.data.frame(t(protein_z_n))
colnames(zprotein)<-proorder$clust_membership
clusterlist <- function(x,z){
  y<-x[,names(x) %in% colnames(x)[grepl(z,colnames(x))]]%>%
   mutate( Donor=c(rep("Donor1",27),rep("Donor2",27)), 
        Weeks=c(rep("Week01",3),rep("Week02",3),rep("Week03",3),rep("Week04",3),rep("Week06",3),
                      rep("Week08",3),rep("Week10",3),rep("Week12",3),rep("Week16",3),
                        rep("Week01",3),rep("Week02",3),rep("Week03",3),rep("Week04",3),rep("Week06",3),
                        rep("Week08",3),rep("Week10",3),rep("Week12",3),rep("Week16",3)))
}

clusterplot <- function (x){
x_long <- melt (x, id = c ("Donor","Weeks"))
ggplot (data = x_long, aes(x=Weeks, y=value,color=Donor)) +
  stat_summary(data=subset(x_long,Donor=="Donor1"),fun.y=median, aes(group=variable), geom="line", size=0.5,color="#f3746e",alpha=0.1)+
  stat_summary(data=subset(x_long,Donor=="Donor2"),fun.y=median, aes(group=variable), geom="line", size=0.5,color="#1fbdc2",alpha=0.1)+
  stat_summary(data=subset(x_long,Donor=="Donor1"),fun.y=median, group=1, geom="line", size=2,color="#b30000")+
  stat_summary(data=subset(x_long,Donor=="Donor2"),fun.y=median, group=1, geom="line", size=2,color="#045a8d")+
  theme_classic(base_size=16)+
  scale_x_discrete(name ="Weeks",labels=c("1", "2","3","4", "6", "8", "10", "12","16"))+
  ylab("Z-score")+
  theme(axis.text.x = element_text(size=30,color="black"),
        axis.text.y = element_text(hjust = 1, size=30,color="black"),
        axis.title.x = element_text(size=30,color="black"),
        axis.title.y = element_text(size=30,color="black"),
        plot.title = element_text(size=40),
        axis.line = element_line(colour = "black",size=1),
        legend.justification = c(0, 0),
        legend.title = element_blank(),
        legend.position="right",
        legend.key.size = unit(5,"mm"),
        legend.text=element_blank())+
  ggsave(paste(x,".pdf"),width = 14, height = 7)
}

cluster1 <<- clusterlist(zprotein,"1")
clusterplot(cluster1)

cluster2 <<- clusterlist(zprotein,"2")
clusterplot(cluster2)

cluster3 <<- clusterlist(zprotein,"3")
clusterplot(cluster3)

cluster4 <<- clusterlist(zprotein,"4")
clusterplot(cluster4)

cluster5 <<- clusterlist(zprotein,"5")
clusterplot(cluster5)

# plot for each protein
# Table for plot
wmboth <- cbind(wm[,c(4:5)],wm_mgml)
wmboth$Gene <- gsub(".*GN=| PE=.*", "", wmboth$Description)
wmboth$Description <- sub("\\OS.*","",wmboth$Description)
write.csv(wmboth,"wmboth.csv")

# save each plot
setwd("D:/Projects to finish/Personalized milk/Data analysis/WM001002/WM001002Proteome") # place to save
plotlist <- wmboth$Description
plotdata <- c()

for(i in 1:length(plotlist)){
  plotdata[[i]] <- as.data.frame(as.numeric(c(t(wmboth[i,c(3:29)]),
                                            t(wmboth[i,c(30:56)]))))%>%
    mutate(Donor=c(rep("Donor1",27),rep("Donor2",27)),
           Weeks=c(rep("Week01",3),rep("Week02",3),rep("Week03",3),rep("Week04",3),rep("Week06",3),
           rep("Week08",3),rep("Week10",3),rep("Week12",3),rep("Week16",3),rep("Week01",3),rep("Week02",3),rep("Week03",3),rep("Week04",3),rep("Week06",3),
           rep("Week08",3),rep("Week10",3),rep("Week12",3),rep("Week16",3)))
  colnames(plotdata[[i]])[1] <- "mgml"
  ggplot(plotdata, aes(x = factor(Weeks), y = mgml, colour = Donor)) + 
         geom_point(aes(shape=Donor),stroke = 2, size = 4) + 
         scale_shape_manual(values=c(17,9)) + 
         stat_summary(fun.y = median, na.rm = T,geom = "line",mapping = aes(group=Donor),size = 2)+
         scale_x_discrete(name ="Weeks",labels=c("Week01" = "1", "Week02" = "2","Week03" = "3",
                                                    "Week04" = "4", "Week06" = "6","Week08" = "8",
                                                    "Week10" = "10", "Week12" = "12","Week16" = "16"))+
		 ggtitle(wmboth$Description[i])+ 
         ylim(0,1.1*max(plotdata[[i]]$mgml))+
         ylab("mg/ml")+
         theme_bw(base_size = 30)+
         theme(axis.text.x = element_text(size=30,color="black"),
          axis.text.y = element_text(hjust = 1, size=30,color="black"),
          axis.title.y = element_text(size=30,color="black"),
          plot.title = element_text(size=40),
          axis.line = element_line(colour = "black",size=1),
          panel.grid = element_blank(),
          panel.border = element_rect(colour=NA),
          legend.justification = c(0, 0),
          legend.title = element_blank(),
          legend.position=c(0.8, 0.8),
          legend.text = element_text(size=18),
          legend.key.size = unit(15,"mm"))+
     ggsave(paste(i,".pdf"),width = 14, height = 7)
}


# Peptidome
# import data to analysis
setwd("D:/Projects to finish/Personalized milk/Peptidome")
pep <- read_excel("Peptidome_WM001_002.xlsx", sheet = "ForR")
pep <- pep[pep$`Ions Score (by Search Engine): Mascot`>=20,] #Mascot score>=20
pep <- pep[,c(1:70)]#with value

#Normalized by sum
pep_percentage <- apply(pep[,c(17:70)], 2, function(x)(x*100/sum(x,na.rm = TRUE)))
pep_percentage <- apply(pep_percentage,2,function(x)(replace(x,is.na(x),0)))
pep_percentage_df <- data.frame(pep[,c(1:16)],pep_percentage)

#Calculated as ug/ml by Total amino acid
pep_ugml <- data.frame(pep_percentage[,c(1:3)]*19.18808*0.025*10,
                    pep_percentage[,c(4:6)]*19.66452*0.025*10,
                    pep_percentage[,c(7:9)]*15.179*0.025*10,
                    pep_percentage[,c(10:12)]*12.72446*0.025*10,
                    pep_percentage[,c(13:15)]*12.84298*0.025*10,
                    pep_percentage[,c(16:18)]*12.18363*0.025*10,
                    pep_percentage[,c(19:21)]*12.60597*0.025*10,
                    pep_percentage[,c(22:24)]*12.95751*0.025*10,
                    pep_percentage[,c(25:27)]*10.63312*0.025*10,
                    pep_percentage[,c(28:30)]*16.25398*0.025*10,
                    pep_percentage[,c(31:33)]*14.30421*0.025*10,
                    pep_percentage[,c(34:36)]*14.05865*0.025*10,
                    pep_percentage[,c(37:39)]*15.25564*0.025*10,
                    pep_percentage[,c(40:42)]*11.0365*0.025*10,
                    pep_percentage[,c(43:45)]*12.59312*0.025*10,
                    pep_percentage[,c(46:48)]*12.28856*0.025*10,
                    pep_percentage[,c(49:51)]*11.11747*0.025*10,
                    pep_percentage[,c(52:54)]*11.64257*0.025*10)
pep_ugml_df <- data.frame(pep[,c(1:16)],pep_ugml)
write.csv(pep_ugml_df,"pep_ugml_all_20190628.csv")

# Heatmap 
pep_z <- apply(pep_ugml,1,function(x)(scale(as.numeric(x),center = TRUE,scale = TRUE)))
pep_z <- as.data.frame(t(pep_z))
colnames(pep_z) <- colnames(pep_ugml)
pep_z$ID <- pep$ID


pep_z_n <- pep_z[!is.na(pep_z$X001_Week.1_1),c(1:54)]
min(pep_z_n)
max(pep_z_n)


hrpep <- hclust(dist(pep_z_n,method="euclidean"), method="complete")
hcpep <- hclust(as.dist(1-cor(pep_z_n)), method="single")

# Heatmap without annoatation
pdf("peptides_heatmap_wm001002.pdf",40,40)
Heatmap(pep_z_n,cluster_columns = as.dendrogram(hcpep),cluster_rows = as.dendrogram(hrpep),
        col = colorRamp2(c(-4,-2,0,2,4), 
                         c("#03294f","#2166ac","#f7f7f7","#c51b7d","#b50140")),
        heatmap_legend_param = list(at = c(-4,-2,0,2,4)))
dev.off()

# cluster and cluster for each peptide
pcuttree <- cutree(hrpep, h = 12.8)
tc_pep <- data.frame(clust_membership = as.character(pcuttree))
str(tc_pep)
tc_pep$order <- c(1:length(hrpep$order))
peporder <- data.frame(order=hrpep$order)
peporder <- merge(peporder,tc_pep,by="order")

pep_ug_1 <- pep_ugml_df[!rowSums(pep_ugml_df[,c(17:70)])==0,]
pep_z_n_id <- data.frame(pep_z[!is.na(pep_z$X001_Week.1_1),55],pep_z_n)

peptideorder <- cbind(peporder,pep_z_n_id,pep_ug_1)
write.csv(peptideorder,"peptideorder_cluster6.csv")

# Heatmap with annoatation
pdf("peptide_heatmap_wm001002_cluster_6.pdf",20,20)
Heatmap(pep_z_n,cluster_columns = as.dendrogram(hcpep),cluster_rows = as.dendrogram(hrpep),
        col = colorRamp2(c(-4,-2,0,2,4), 
                         c("#03294f","#2166ac","#f7f7f7","#c51b7d","#b50140")),
        heatmap_legend_param = list(at = c(-4,-2,0,2,4)))+
  rowAnnotation(tc_pep,annotation_legend_param =tc_pep$clust_membership)
dev.off()    

#overlayed lineplots by cluster membership (by each peptide and emphasized median)
zpep <- as.data.frame(t(pep_z_n))
colnames(zpep) <- peporder$clust_membership
pep_cluster1 <<- clusterlist(zpep,"1")
clusterplot([pep_cluster1)

pep_cluster2 <<- clusterlist(zpep,"2")
clusterplot([pep_cluster2)

pep_cluster3 <<- clusterlist(zpep,"3")
clusterplot([pep_cluster3)

pep_cluster4 <<- clusterlist(zpep,"4")
clusterplot([pep_cluster4)

pep_cluster5 <<- clusterlist(zpep,"5")
clusterplot([pep_cluster5)

pep_cluster6 <<- clusterlist(zpep,"6")
clusterplot([pep_cluster6)

# plot for each peptide

# Table for plot
pep_ugml_plot <- pep_ugml_df%>%
  mutate(Description=sub("\\OS.*","",Description),Nameforplot=paste(Description,ID))
write.csv(pep_ugml_plot,"pep_ugml_plot.csv")

# save each plot
setwd("D:/Projects to finish/Personalized milk/Peptidome/individual")
plotlistpep<-pep_ugml_plot$Nameforplot
plotdatapep<-c()
for(i in 1:length(plotlistpep)){
  plotdatapep[[i]]<-as.data.frame(as.numeric(c(t(pep_ugml_df[i,c(17:70)]))))%>%
    mutate(Donor=c(rep("Donor1",27),rep("Donor2",27)),
           Weeks=c(rep("Week01",3),rep("Week02",3),rep("Week03",3),rep("Week04",3),rep("Week06",3),
                   rep("Week08",3),rep("Week10",3),rep("Week12",3),rep("Week16",3),rep("Week01",3),rep("Week02",3),rep("Week03",3),rep("Week04",3),rep("Week06",3),
                   rep("Week08",3),rep("Week10",3),rep("Week12",3),rep("Week16",3)))
  colnames(plotdatapep[[i]])[1]<-"ugml"
  ggplot(plotdatapep[[i]], aes(x = factor(Weeks), y = ugml, colour = Donor)) + 
    geom_point(aes(shape=Donor),stroke=2,size=4) + 
    scale_shape_manual(values=c(17,9)) +
    scale_color_manual(values=c("#fbaa19", "#3867b1"))+
    stat_summary(fun.y = median, na.rm = T,geom = "line",mapping = aes(group=Donor),size=2)+
    ggtitle(pep_ugml_plot$Nameforplot[i])+ 
    ylim(0,1.1*max(plotdatapep[[i]]$ugml))+
    scale_x_discrete(name ="Weeks",labels=c("Week01" = "1", "Week02" = "2","Week03" = "3",
                                            "Week04" = "4", "Week06" = "6","Week08" = "8",
                                            "Week10" = "10", "Week12" = "12","Week16" = "16"))+
    ylab("ug/ml")+
    theme_bw(base_size = 30)+
    theme(axis.text.x = element_text(size=30,color="black"),
          axis.text.y = element_text(hjust = 1, size=30,color="black"),
          axis.title.y = element_text(size=30,color="black"),
          plot.title = element_text(size=12),
          axis.line = element_line(colour = "black",size=1),
          panel.grid = element_blank(),
          panel.border = element_rect(colour=NA),
          legend.justification = c(0, 0),
          legend.title = element_blank(),
          legend.position=c(0.8, 0.8),
          legend.text = element_text(size=18),
          legend.key.size = unit(15,"mm"))+
    ggsave(paste(i,".pdf"),width = 14, height = 7,limitsize = FALSE)
}



