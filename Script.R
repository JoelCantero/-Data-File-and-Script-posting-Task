library(readxl)
Russet_ineqdata <- read.table("Data-File-and-Script-posting-Task/Russet_ineqdata.txt", sep= '\t')
Russet_ineqdata$demo <- as.factor(Russet_ineqdata$demo)
levels(Russet_ineqdata$demo) = c("Stable", "Instable", "Dictatorship")
library(mice)

#Imputation By Mice And defined X as matrix with continous values.
imputed_data <- mice(Russet_ineqdata)
imputed_data <- complete(imputed_data)
row.names(imputed_data) <- row.names(Russet_ineqdata)
X <- as.matrix(imputed_data[,0:8])



#a. Define the matrix N of weights of individuals (with uniform weights).
weights <- rep(1,nrow(X))
N <- diag(weights/sum(weights))
#b. Compute the centroid G of individuals.
G <-  colMeans(X)
#c. Compute the covariance or correlation matrix of X (be aware of dividing by sum(weights_i)).
X.cor <- cor(X)
X.centered <- 
X.cor.manual <- t(X)%*%N%*%X
X.cov <- cov(X)
X.cov.manual
#d. Compute the centered X matrix and standardized X matrix.
#e. Diagonalize X centered.
#f. Do the screeplot of the eigenvalues and define the number of significant dimensions. How much is the retained information?
#g. Compute the projections of individuals in the significant dimensions.
#h. Compute the projection of variables in the significant dimensions
#i. Plot the individuals in the first factorial plane of Rp. Color the individuals according the ???demo??? variable.
#j. Plot the variables (as arrows) in the first factorial plane of Rn.
#k. According to the Russet complete data, justify which metric M is appropriate for this problem.
#l. Compute the correlation of the variables with the significant principal components and interpret them.
 
 
 
 