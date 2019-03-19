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
X.cov <- cov(X)/sum(weights)
X.centered <- scale(X, scale=FALSE)
X.cov.manual <- (t(X.centered)%*%N%*%X.centered)/sum(weights)
#View(X.cov - X.cov.manual)

X.cor <- cor(X)/sum(weights)
X.standarized <- scale(X, center = FALSE, scale=TRUE)
X.cor.manual <- (t(X.standarized)%*%N%*%X.standarized)/sum(weights)
#View(X.cor - X.cor.manual)

#d. Compute the centered X matrix and standardized X matrix.
X.centered <- scale(X, scale=FALSE)
X.standarized <- scale(X, center = FALSE, scale=TRUE)
#e. Diagonalize X centered and standarized .

X.centered.eig <- eigen(t(X.centered)%*%N%*%X.centered)
X.centered.eig$values

X.standarized.eig <- eigen(t(X.standarized)%*%N%*%X.standarized)
X.standarized.eig$values
#f. Do the screeplot of the eigenvalues and define the number of significant dimensions. 
# How much is the retained information?
plot(X.centered.eig$values, type = "l")
plot(X.standarized.eig$values, type = "l")
#g. Compute the projections of individuals in the significant dimensions.

#h. Compute the projection of variables in the significant dimensions.

#i. Plot the individuals in the first factorial plane of Rp. Color the individuals according the ???demo??? variable.
#j. Plot the variables (as arrows) in the first factorial plane of Rn.
#k. According to the Russet complete data, justify which metric M is appropriate for this problem.
#l. Compute the correlation of the variables with the significant principal components and interpret them.
 
 
 
 