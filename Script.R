library(mice)
source("PCA.R")

X <- read.table('Russet_ineqdata.txt', header=T, sep='\t', row.names=1)
X$demo <- as.factor(X$demo)
levels(X$demo) = c("Stable", "Instable", "Dictatorship")
md.pattern(X)
imputedX <- complete(mice(X))
row.names(imputedX) <- row.names(X)
X.matrix <- as.matrix(imputedX[,0:8])
weights <- rep(1,nrow(X.matrix))
weights <- weights/sum(weights) #Normalized 
N <- diag(weights) 
G = vector()
for (i in 1:ncol(X.matrix)){
  G = c(G,weighted.mean(X.matrix[,i],weights))
}
G
X.centered <- X.matrix - rep(G, rep.int(nrow(X.matrix), ncol(X.matrix)))
X.covariance = (t(X.centered)%*%N%*%X.centered)/sum(weights)

X.standarized <- scale(X.matrix)
X.correlation <- (t(X.standarized)%*%N%*%X.standarized)/sum(weights)

normalized <- TRUE

if (normalized) {
  X <- X.standarized
  S <- X.correlation
} else {
  X <- X.centered
  S <- X.covariance
}
X.centered <- X - rep(G, rep.int(nrow(X), ncol(X)))
X.standarized <- scale(X, center = FALSE, scale=TRUE)
X.svd <- svd(S)
X.values <- X.svd$d
X.vector <- X.svd$u
plot(X.values, type = "l", main="Screeplot of Russet")
percentage <- ((X.values[1]+X.values[2])/sum(X.values))
percentage
PSI <- X %*% X.vector
PHI <- sqrt(X.values)*X.vector

plot(PSI, col=c('green', 'blue', 'red')[imputedX$demo])
legend(x="topright", legend=c('Stable', 'Instable', 'Dictatorship'), fill=c('green', 'blue', 'red'))

originx <- rep(0, ncol(PHI))
originy <- rep(0, ncol(PHI))
X = imputedX[,1:8]
etiq = names(X)
plot(PHI, xlim=c(0,0), ylim=c(0,0))
text(PHI,labels=etiq, col="black")
arrows(originx, originy, PHI[,1], PHI[,2])
abline(h=0,v=0,col="gray")


COR <-cor(X, PSI)
COR


weights <- rep(1,nrow(X.matrix))
weights[11] <-0
weights <- weights/sum(weights) #Normalized
Cor.Cuba <- PCA(X.matrix, weights, eucledian=TRUE)
Cor.Cuba



X.FactoMineR <- imputedX
pca <- FactoMineR::PCA(graph=T, X.FactoMineR, ncp=8, quali.sup=9)
plot(pca, cex=0.8, habillage='demo')


firstFactorialPlane <- sort(abs(pca$ind$cos2[,1])+abs(pca$ind$cos2[,2]), decreasing=TRUE)
bestCountry <- firstFactorialPlane[1]
worseCountry <- firstFactorialPlane[length(firstFactorialPlane)]
bestCountry
worseCountry

contribPCA <- pca$ind$contrib
contribPCA.1 <- sort(abs(contribPCA[,1]), decreasing=TRUE)
contribPCA.2 <- sort(abs(contribPCA[,2]), decreasing=TRUE)
bestCountries.1 <- contribPCA.1[1:3]
bestCountries.2 <- contribPCA.2[1:3]
bestCountries.1
bestCountries.2


variableRepresented <- sort(abs(pca$var$cos2[,1])+abs(pca$var$cos2[,2]), decreasing=TRUE)
best <- variableRepresented[1]
worse <- tail(variableRepresented, n=1)
best
worse

firstPrincipalComponent <- sort(abs(pca$var$contrib[,1]), decreasing=TRUE)
secondPrincipalComponent <- sort(abs(pca$var$contrib[,2]), decreasing=TRUE)
mostInfluencingVariablesFirst <- firstPrincipalComponent[1:3]
mostInfluencingVariablesSecond <- secondPrincipalComponent[1:3]
mostInfluencingVariablesFirst
mostInfluencingVariablesSecond

significantModalitiesDemo <- sort(abs(pca$quali.sup$cos2[,1]) + abs(pca$quali.sup$cos2[,2]), decreasing=TRUE)
significantModalitiesDemo

library(nipals)
nipals <- nipals(X.standarized, 3)
scores <- nipals$scores
loadings <- nipals$loadings
biplot(scores, loadings)


PSI = pca$ind$coord[,1:2]
PHI = pca$var$coord[,1:2]
PC.ROT <- varimax(PHI)
PHI.ROT = PC.ROT$loadings[1:ncol(X),]
X = imputedX[,1:8]
Xs = scale(X)
iden = row.names(X)
etiq = names(X)
lmb.rot = diag(t(PC.ROT$loadings) %*% PC.ROT$loadings)
sum(lmb.rot)
sum(pca$eig[1:2,])
Psi_stan.rot = Xs %*% solve(cor(X)) %*% PHI.ROT
PSI.ROT = Psi_stan.rot %*% diag(sqrt(lmb.rot))
library(calibrate)
ze = rep(0,ncol(X))
plot(PHI.ROT,type="n",xlim=c(-1,1),ylim=c(-1,1))
text(PHI.ROT,labels=etiq, col="blue")
arrows(ze, ze, PHI.ROT[,1], PHI.ROT[,2], length = 0.07,col="blue")
abline(h=0,v=0,col="gray")
circle(1)

plot(PSI.ROT,type="n")
text(PSI.ROT,labels=iden,col=as.numeric(imputedX$demo))
abline(h=0,v=0,col="gray")

pca$ind$coord <- pca$ind$coord[,1:2] %*% PC.ROT$rotmat
FactoMineR::dimdesc(pca, axes=1:2)