# Applied Machine Learning for Health Data
# Assignment 3
# Author: Steve Hawley
# Date: Dec 1, 2019

nuc <- read.csv("Data_Cortex_Nuclear.csv")

str(nuc)
summary(nuc)

#remove rows and columns with >10% missing data
nuc <- nuc[,-which(colMeans(is.na(nuc))>0.1)]
nuc <- nuc[-which(rowMeans(is.na(nuc))>0.1),]

#replace remaining NAs with column means
x <- unname(which(unlist(lapply(nuc, is.numeric))))
for(i in x){
  nuc[is.na(nuc[,i]), i] <- mean(nuc[,i], na.rm = TRUE)
}

