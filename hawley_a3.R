# Applied Machine Learning for Health Data
# Assignment 3
# Author: Steve Hawley
# Date: Dec 1, 2019

##########################################
########### Data Preparation #############

#read in data set
nuc <- read.csv("Data_Cortex_Nuclear.csv")

#take a look at the data
str(nuc)
summary(nuc)

#set up data for cleaning
nuc.cl <- nuc

#remove rows and columns with >10% missing data
nuc.cl <- nuc.cl[,-which(colMeans(is.na(nuc.cl))>0.1)]
nuc.cl <- nuc.cl[-which(rowMeans(is.na(nuc.cl))>0.1),]

#replace remaining NAs with column means *within same class*
x <- unname(which(unlist(lapply(nuc.cl, is.numeric))))
for(i in x){
  nuc.cl[,i]<-ave(nuc.cl[,i],nuc.cl$class,FUN=function(y) 
    ifelse(is.na(y), mean(y,na.rm=TRUE), y))
}

summary(nuc.cl)

#scale the feature space (i.e., gene expression) for cluster analysis
sc <- scale(nuc.cl[,2:73])

##########################################
########## k means clustering ############

library(factoextra) #for cluster plots

# K-Means Cluster Analysis. Starting with 8 clusters to align with the 8 classes
fit.8 <- kmeans(sc, 8)
#plot clusters using principal components (fviz_cluster takes first to PCs by default)
fviz_cluster(list(data = sc, cluster = fit.8$cluster))
table(fit.8$cluster, nuc.cl$class)

#Will try a 2 cluster solution to look at other experimental conditions
fit.2 <- kmeans(sc, 2)
fviz_cluster(list(data = sc, cluster = fit.2$cluster))
table(fit.2$cluster,nuc.cl$Genotype)
table(fit.2$cluster,nuc.cl$Treatment)
table(fit.2$cluster,nuc.cl$Behavior)

#Will try a 2 cluster solution on the genes instead of the mice
fit.2.g <- kmeans(t(sc), 2)
fviz_cluster(list(data = t(sc), cluster = fit.2.g$cluster))


##########################################
################## PAM ###################

library(cluster) #for PAM

#will start with 8 cluster solution
pam.8 <- pam(sc,k = 5)
fviz_cluster(list(data = sc, cluster = pam.8$clustering))
table(pam.8$cluster, nuc.cl$class)

#Will try a 2 cluster solution to look at other experimental consitions
pam.2 <- pam(sc, 2)
fviz_cluster(list(data = sc, cluster = pam.2$cluster))
table(pam.2$cluster,nuc.cl$Genotype)
table(pam.2$cluster,nuc.cl$Treatment)
table(pam.2$cluster,nuc.cl$Behavior)


##########################################
##### Hierarchical clustering ############

library(gplots) #for heatmap.2

#Will attempt 8 cluster solution with complete linkage, then average, then Ward's
#Complete
hc.complete.8 <- hclust(dist(sc), method="complete")
plot(hc.complete.8, main="Complete Linkage", xlab="", sub="",cex =.9)
cut.tree.com.8 <- cutree(hc.complete.8,k=8)
fviz_cluster(list(data = sc, cluster = cut.tree.com.8),ellipse.type = "norm")
table(cut.tree.com.8,nuc.cl$class)

#completing one example heatmap
mycolhc <- rainbow(length(unique(cut.tree.com.8)), start=0.1, end=0.9); 
mycolhc <- mycolhc[as.vector(cut.tree.com.8)]
heatmap.2(sc, Rowv=as.dendrogram(hc.complete.8), trace="none",RowSideColors=mycolhc) 

#Average
hc.average.8 <- hclust(dist(sc), method="average")
plot(hc.average.8, main="Average Linkage", xlab="", sub="",cex =.9)
cut.tree.avg.8 <- cutree(hc.average.8,k=8)
fviz_cluster(list(data = sc, cluster = cut.tree.avg.8))
table(cut.tree.avg.8,nuc.cl$class)

#Wards
hc.ward.8 <- hclust(dist(sc), method="ward.D")
plot(hc.ward.8, main="Ward's Linkage", xlab="", sub="",cex =.9)
cut.tree.ward.8 <- cutree(hc.ward.8,k=8)
fviz_cluster(list(data = sc, cluster = cut.tree.ward.8))
table(cut.tree.ward.8,nuc.cl$class)

#Will repeat the above with a 2 cluster solution 
#Complete
hc.complete.2 <- hclust(dist(sc), method="complete")
plot(hc.complete.2, main="Complete Linkage", xlab="", sub="",cex =.9)
cut.tree.com.2 <- cutree(hc.complete.2,k=2)
fviz_cluster(list(data = sc, cluster = cut.tree.com.2))
table(cut.tree.com.2,nuc.cl$Treatment)
table(cut.tree.com.2,nuc.cl$Genotype)
table(cut.tree.com.2,nuc.cl$Behavior)

#Average
hc.average.2 <- hclust(dist(sc), method="average")
plot(hc.average.2, main="Average Linkage", xlab="", sub="",cex =.9)
cut.tree.avg.2 <- cutree(hc.average.2,k=2)
fviz_cluster(list(data = sc, cluster = cut.tree.avg.2))
table(cut.tree.avg.2,nuc.cl$Treatment)
table(cut.tree.avg.2,nuc.cl$Genotype)
table(cut.tree.avg.2,nuc.cl$Behavior)

#Wards
hc.ward.2 <- hclust(dist(sc), method="ward.D")
plot(hc.ward.2, main="Ward's Linkage", xlab="", sub="",cex =.9)
cut.tree.ward.2 <- cutree(hc.ward.2,k=2)
fviz_cluster(list(data = sc, cluster = cut.tree.ward.2))
table(cut.tree.ward.2,nuc.cl$Treatment)
table(cut.tree.ward.2,nuc.cl$Genotype)
table(cut.tree.ward.2,nuc.cl$Behavior)
