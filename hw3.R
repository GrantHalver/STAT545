#Question 1
#  thelog-likelihood function
llCauchy <- function(mu, x){
  sum(dt(x-mu, df=1, log=TRUE))
}

# gradient function of the log-likelihood function
g.llCauchy <- function(mu, x){
  2 * sum((x-mu)/(1+(x-mu)^2))
}

# hessian function of the log-likelihood function
h.llCauchy <- function(mu, x){
  2 * sum( (-(1+(x-mu)^2) + 2*(x-mu)^2) /(1+(x-mu)^2)^2)
}

newton.raphsonCauchy <- function(X, eps = 1.0E-10) {
  X <- sort(X)
  mut <- median(X)
  while(abs(g.llCauchy(mu=mut,x=X)) > eps){
    mut <- mut - (g.llCauchy(mu=mut,x=X)/h.llCauchy(mu=mut,x=X))
  }
  mut
}

#Question 2
mu0 <- 8
iter <-100
muhat <- numeric(iter)
crLower <- numeric(iter)
for(i in 1:iter){
  X <- rt(n=5, df=1) + mu0
  muhat[i] <- newton.raphsonCauchy(X)
}

biashat = mean(muhat)-8
varmuhat = var(muhat)
msemuhat = biashat^2+varmuhat

biashat
varmuhat
msemuhat


#Question 3
bisect <- function(fn, a, b, eps = 1.0E-10, ...){
  fa <- fn(a, ...) 
  fb <- fn(b, ...) 
  if(fa*fb > 0) {
    warning("invalid arguments: fn(a) * fn(b) > 0")
    return(NULL)
  }
  while(abs(b-a) > eps) {
    m <- (a+b)/2 
    fm <- fn(m, ...)
    if(fb*fm > 0){
      b <- m
      fb <- fm
    } else {
      a <- m
      fa <- fm
    }
  }
  m <- (a+b)/2
  fm <- fn(m, ...) 
  list(a=a, b=b, m=m, value=fm)
}

global.optim <- function(x){
  n <- 5 # for now
  if(length(x) != n) stop("length(x) != n")
  x <- sort(x)
  
  # try interval x[i:i+1]
  rs <- data.frame(a = NULL, b = NULL, m = NULL, value = NULL)
  for(i in 1:(n-1)){
    if(g.llCauchy(x[i], x)*g.llCauchy(x[i+1], x) <= 0) {
      rs <- rbind(rs, bisect(g.llCauchy, a=x[i], b=x[i+1], eps = 1.0E-10, x=x))
    }
  }
  
  rs <- rs[complete.cases(rs),]
  rs <- rs[order(rs[ ,4],decreasing=FALSE),]
  list(a=rs[1,1], b=rs[1,2], m=rs[1,3], value=rs[1,4])
}

mu0 <- 8
iter <-100
muhat <- data.frame(a = NULL, b = NULL, m = NULL, value = NULL)
for(i in 1:iter){
  X <- rt(n=5, df=1) + mu0
  muhat <- rbind(muhat, global.optim(X))
}

biashat = mean(muhat$m)-8
varmuhat = var(muhat$m)
msemuhat = biashat^2+varmuhat

biashat
varmuhat
msemuhat

#Question 4
# gradient function of the log-likelihood function
g.llCauchy2 <- function(theta, x){
  r <- matrix(0, nrow = 2, ncol = 1)
  u <- theta[1,1]
  s <- sqrt(theta[2,1])
  r[1,1] <- sum(1/pi*(2 * s*(-u + x))/(s^2 + (u - x)^2)^2)
  r[2,1] <- sum(1/pi*(-s^2 + (-u + x)^2)/(s^2 + (-u + x)^2)^2)
  r
}

# hessian function of the log-likelihood function
h.llCauchy2 <- function(theta, x){
  z = (x-theta[1,1])/theta[2,1]
  r <- matrix(0, nrow = 2, ncol = 2)
  u <- theta[1,1]
  s <- sqrt(theta[2,1])
  r[1,1] <- sum(1/pi*((8 *(-u + x)^2)/(s^4 *(1 + (-u + x)^2/s^2)^3) - 2/(s^2 * (1 + (-u + x)^2/s^2)^2))/s)
  r[1,2] <- sum(1/pi*(2 *(-3 *s^2 + (u - x)^2) *(-u + x))/(s^2 + (u - x)^2)^3)
  r[2,1] <- r[1,2]
  r[2,2] <- sum(1/pi*(2* s* (s^2 - 3 *(u - x)^2))/(s^2 + (u - x)^2)^3)
  r
}

newton.raphsonCauchy2 <- function(X, eps = 1.0E-10) {
  require(MASS)
  X <- sort(X)
  mut <- median(X)
  sigma2t <- IQR(X)/2
  thetat <- matrix(c(mut,sigma2t), nrow = 2, ncol = 1)
  while(abs(g.llCauchy2(theta=thetat,x=X)[1,1]) > eps || abs(g.llCauchy2(theta=thetat,x=X)[2,1]) > eps){
    thetat <- thetat - (ginv(h.llCauchy2(theta=thetat, x=X))%*%g.llCauchy2(theta=thetat, x=X))
  }
  thetat
}

library(MASS)
X <- rt(n=5, df=1) + mu0
newton.raphsonCauchy2(X)
