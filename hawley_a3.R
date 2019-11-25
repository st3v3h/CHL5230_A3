# Applied Machine Learning for Health Data
# Assignment 3
# Author: Steve Hawley
# Date: Dec 1, 2019

nuc <- read.csv("Data_Cortex_Nuclear.csv")

str(nuc)
summary(nuc)

#set up data for cleaning
nuc.cl <- nuc

#remove rows and columns with >10% missing data
nuc.cl <- nuc.cl[,-which(colMeans(is.na(nuc.cl))>0.1)]
nuc.cl <- nuc.cl[-which(rowMeans(is.na(nuc.cl))>0.1),]

#replace remaining NAs with column means within same class
x <- unname(which(unlist(lapply(nuc.cl, is.numeric))))
for(i in x){
  nuc.cl[,i]<-ave(nuc.cl[,i],nuc.cl$class,FUN=function(y) 
    ifelse(is.na(y), mean(y,na.rm=TRUE), y))
}

summary(nuc.cl)

#scale the feature space for cluster analysis
sc <- scale(nuc.cl[,2:73])

##########################################
########## k means clustering ############

# K-Means Cluster Analysis
fit <- kmeans(sc, 8, nstart = 20) # 8 cluster solution
# get cluster means
aggregate(sc,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(sc, fit$cluster) 

clusplot(sc, fit$cluster, color=TRUE, shade=TRUE,
         labels=2, lines=0)

table(fit$cluster, nuc.cl$class)

fit.2 <- kmeans(sc, 2)
table(fit.2$cluster,nuc.cl$Genotype)
table(fit.2$cluster,nuc.cl$Treatment)
table(fit.2$cluster,nuc.cl$Behavior)


##########################################
##### Hierarchical clustering ############

hc.complete <- hclust(dist(sc), method="complete")
plot(hc.complete,main="Complete Linkage", xlab="", sub="",cex =.9)
cut.tree.com <- cutree(hc.complete,k=8)
table(cut.tree.com,nuc.cl$class)


# iris.hc.complete <- hclust(dist(scale(nuc.cl[,2:73])), method="complete")
# plot(iris.hc.complete,main="Complete Linkage", xlab="", sub="", cex =.9)
# cut.tree.com <- cutree(iris.hc.complete,k=8)
# table(cut.tree.com,nuc.cl$Treatment)
# table(cut.tree.com,nuc.cl$Genotype)
# table(cut.tree.com,nuc.cl$class)


## Row- and column-wise clustering 
hr <- hclust(as.dist(1-cor(t(sc), method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(sc, method="spearman")), method="complete") 
## Tree cutting
# mycl <- cutree(hr, h=max(hr$height)/1.5); 
mycl <- cutree(hr, k=2) 
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); 
mycolhc <- mycolhc[as.vector(mycl)] 
## Plot heatmap 
mycol <- colorpanel(40, "darkblue", "yellow", "white") # or try redgreen(75)
heatmap.2(sc, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol, scale="row", density.info="none", trace="none", RowSideColors=mycolhc) 



##########################################
################## PCA ###################

pr.out <- prcomp(sc, scale=F)
pr.out$rotation

# bivariate plot, plotting together the points and the features based on the first 2 pc's
biplot(pr.out, scale=0)
pr.out$sdev
pr.var <- pr.out$sdev^2
pr.var
pve <- pr.var/sum(pr.var) #proportion that each pc explains variance
pve
plot(pve, xlab="Principal Component", ylab="Proportion of
     Variance Explained ", ylim=c(0,1) ,type="b")
plot(cumsum(pve), xlab="Principal Component", ylab="
     Cumulative Proportion of Variance Explained ", ylim=c(0,1) ,
     type="b")

##########################################
################## PAM ###################

library(cluster)
library(fpc)

pk <- pamk(nuc.cl[,2:73])
pks <- pamk(sc)
pks$nc

pr4 <- pam(sc,k = 8) #from cluster package; pam is assuming you want to use euclidian distance
si <- silhouette(pr4)
ssi <- summary(si)
ssi #cluster sizes add up to 64 units
plot(si,col = 1:8)
