#Problem 1
ginv <- function(X, tol=1.0E-08){
  udv <- svd(X, nu = nrow(X), nv = ncol(X))
  udv$d[abs(udv$d) < tol] <- 0
  k1 <- which(udv$d != 0)
  
  V1 <- udv$v[,k1, drop=FALSE]
  U1 <- udv$u[,k1, drop=FALSE]
  
  # Moore-Penrose inverse
  mp.X <- V1 %*% diag( 1/udv$d[k1], nrow=length(k1)) %*% t(U1)
  mp.X
}

A <- matrix(c(1,1,1,1), nrow=2, ncol=2)
ginv(A)


B <- matrix(c(1,2,3,4), nrow=2, ncol=2)
ginv(B)


C <- matrix(c(1,1,1,1,1,1,1,1,1), nrow=3, ncol=3)
ginv(C)


D <- matrix(c(1,2,3,4,5,6,7,8,9), nrow=3, ncol=3)
ginv(D)


#Problem 2
cholesky <- function(x){
  L = matrix(0, nrow = 2, ncol=2)
  L[1,1] = sqrt(x[1,1])
  sum = L[1,1] * L[2,1]
  L[1,2] = (1.0 / L[1,1] * (x[2,1] - sum))
  sum = L[2,1] * L[2,1]
  L[2,2] = sqrt(x[2,2] - sum)
  
  L
}

E <- matrix(c(4,2,2,4), nrow=2, ncol=2)
cholesky(E)


G <- matrix(c(1,0,0,1), nrow=2, ncol=2)
cholesky(G)


H <- matrix(c(3,1,1,2), nrow =2, ncol = 2)
cholesky(H)

#Problem 5
invdet <- function(x){
  library(matlib)
  inv <- -1*swp(x,1:ncol(x))
  det <- x[1,1]
  for(i in 2:ncol(x)){
    det <- det*swp(x,1:i-1)[i,i]
  }
  list(inv, det)
}

invdet(E)
invdet(G)

J <- matrix(c(2,-1,0,-1,2,-1,0,-1,2), ncol = 3, nrow =3)
invdet(J)
