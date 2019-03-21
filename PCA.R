PCA<-function(X, weights, eucledian=TRUE) {
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
  
  euclidean <- TRUE
  
  if (euclidean) {
    X <- X.standarized
    S <- X.correlation
  } else {
    X <- X.centered
    S <- X.covariance
  }
  X.centered <- X - rep(G, rep.int(nrow(X), ncol(X)))
  X.standarized <- scale(X, center = FALSE, scale=TRUE)
  X.svd <- svd(S)
  X.values = X.svd$d
  X.vector <- X.svd$u
  plot(X.values, type = "l", main="Screeplot of Russet")
  PSI <- X %*% X.vector
  PHI <- sqrt(X.values)*X.vector
  plot(PSI, col=c('green', 'orange', 'red')[imputedX$demo])
  legend(x="topright", legend=c('Stable', 'Instable', 'Dictatorship'), fill=c('green', 'orange', 'red'))
  originx <- rep(0, ncol(PHI))
  originy <- rep(0, ncol(PHI))
  plot(PHI, xlim=c(min(PHI), max(PHI)), ylim=c(min(PHI), max(PHI)))
  arrows(originx, originy, PHI[,1], PHI[,2])
  COR <-cor(X, PSI)
  plot(COR)
  arrows(originx, originy, COR[,1], COR[,2])
  return(COR)
}