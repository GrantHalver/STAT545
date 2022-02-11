#Problem 1
gsweep = function (S, K) {
  # sweep operator in STAT 553
  # calculated using "operation in place"
  
  det = attr(S, "det")
  ln.det = if(is.null(det)) 0.0 else log(det)
  
  for(k in K){
    flag = sign(k)  # 1: sweep, -1: reverse sweep
    k = abs(k)      # pivotal index
    a=S[k,k]        # pivotal value
    if(flag*a <= 0) stop("invalid pivotal index k = ", k) 
    ln.det = ln.det + log(flag*a)
    S[k,k]= -1/a
    S[-k,-k]=S[-k,-k]-S[-k,k,drop=F]%*%S[k,-k, drop=F]/a
    a = a*flag
    S[k,-k]=S[k,-k]/a
    S[-k,k]=S[-k,k]/a
  }
  
  structure(S, det = exp(ln.det))
}

em.norm <- function(data, max.it=1){
  
  data <- data.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  
  # Start values
  mu <- rep(0, p)
  Sigma <- diag(p)
  
  
  for(it in 1:max.it){
    # E-step: compute the complete-data sufficient
    # statistics
    
    S.x <- rep(0, p)
    S.xx <- matrix(0, nrow=p, ncol=p)
    
    for(i in 1:n){
      X.i <- data[i,]
      
      mis <- which(is.na(X.i))
      # cat("[", i, "] X.i", X.i, " mis: ", mis,"\n")
      obs <- which(!is.na(X.i))
      # cat("[", i, "] X.i", X.i, " obs: ", obs,"\n")
      
      S <- Sigma
      S <- gsweep(S, obs)
      
      X.i[mis] <- mu[mis] + S[mis, obs, drop=F] %*%( X.i[obs] - mu[obs])
      S.x <- S.x + X.i
      
      CS <- matrix(0, nrow=p, ncol=p)
      CS[mis, mis] <- S[mis, mis]
      
      CX.i <- rep(0, p)
      CX.i[mis] <- X.i[mis]
      
      S.xx <- S.xx +  (CS + outer(X.i, X.i))
      
      
      
    }
    
    # M-step
    
    mu <- S.x/n
    Sigma <- S.xx/n - (S.x/n) %*% t(S.x/n)
    
    #    cat("[", it, "]", mu, "\n")
    #    cat("Sigma:\n") print(Sigma)
    #    cat("S.xx:\n") print(S.xx)
  }
  
  list(mu = mu, Sigma=Sigma)
  
}