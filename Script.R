library(readxl)
Russet_ineqdata <- read.table("-Data-File-and-Script-posting-Task/Russet_ineqdata.txt", sep= '\t')
Russet_ineqdata$demo <- as.factor(Russet_ineqdata$demo)
library(mice)

#Imputation By Mice
X <- mice(Russet_ineqdata)
X <- complete(X)
row.names(X) <- row.names(Russet_ineqdata)
X <- X[,0:8]
