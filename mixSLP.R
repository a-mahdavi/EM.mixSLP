# A function to perform EM estimation for the Sigmoid-Transformed Laplace mixtuers (SLPM)
# g is the number of components and family type can be "LP" for Laplace family or "LN" for logit-normal family
# w is a vector of the weighes for each components
mix.SLP <- function(y, g=1, w=1, mu, s, family="LP", iter.max=100, tol=10^-6, get.init = TRUE, group=T){  
  begin <- proc.time()[3]  ; y[which(y==0|y==1)]=10^(-30)
  dLLP<-function (x, mu = 0, sigma = 1) 
  {
    ql <- qlogis(x)
    ifelse(x <= 0 | x >= 1, 0, extraDistr::dlaplace(ql, mu, sigma)/x/(1 - x))
  }
  dmixLLP <- function(y, w, mu, s){
    d <- 0 ; g <- length(w)
    for ( j in 1:g)
      d <- d + w[j]*dLLP(y, mu[j], s[j])
    return(d) }
  dLN<-function (x, mu = 0, sigma = 1) 
  {
    ql <- qlogis(x)
    d <-  dnorm(ql, mu, sigma)/x/(1 - x)
	d[which(d==0)] <- 10^-30
	return(d)
  }
  dmixLN <- function(y, w, mu, s){
    d <- 0 ; g <- length(w)
    for ( j in 1:g)
      d <- d + w[j]*dLN(y, mu[j], s[j])
    return(d) }
  n <- length(y)   ;     dif <- 1;        count <- 0 
  if (get.init == TRUE) {
    init <- kmeans(qlogis(y), g,  algorithm="Hartigan-Wong")
    w <- init$size/n ;mu <- s <- NULL
	for( j in 1:g){
	x <- y[init$cluster==j]
	mu[j] <- mean(qlogis(x))
	s[j] <- 1/length(x)*sum((qlogis(x)-mu[j])^2)
		}
    nu <- runif(g, 1,5)
 	 }
  if(family=="LP"){
    LL <- 1 
    while ((dif > tol) && (count <= iter.max)) {
      z.hat <- matrix(0,n,g)
      # E step
      for (j in 1:g){
        z.hat[,j] <- w[j]*dLLP(y,mu[j],s[j])/dmixLLP(y, w, mu, s)
        eta=(qlogis(y)-mu[j])/s[j]
        tau.hat <- 1/abs(eta)
        # MCE steps
        w[j] <- sum(z.hat[,j])/n
        mu[j] <- sum(z.hat[,j]*tau.hat*qlogis(y))/sum(z.hat[,j]*tau.hat)
        s[j] <- sqrt(sum(z.hat[,j]*(qlogis(y)-mu[j])^2*tau.hat)/sum(z.hat[,j]))
       }
      LL.new <- sum(log(dmixLLP(y,w,mu,s))) # log-likelihood function
      count <- count +1 
      dif <- abs(LL.new/LL-1)
      LL <- LL.new
      cat('iter =', count, '\tloglike =', LL.new, '\n')
    } 
    aic <- -2 * LL.new + 2 * (2*g+g-1)
    bic <- -2 * LL.new + log(n) * (2*g+g-1)
    edc <- -2 * LL.new + 0.2*sqrt(n) * (2*g+g-1)
    end <- proc.time()[3]
    time <- end-begin
    obj.out <- list(family=family,w=w, mu=mu, sigma=s , loglik=LL.new, aic=aic, bic=bic, edc=edc, iter=count,elapsed=as.numeric(time))
  }
  if(family=="LN"){
    LL <- 1 
    while ((dif > tol) && (count <= iter.max)) {
      z.hat  <- matrix(0,n,g)
      # E step
      for (j in 1:g){
        z.hat[,j] <- w[j]*dLN(y,mu[j],s[j])/dmixLN(y, w, mu, s)
        # MCE steps
        w[j] <- sum(z.hat[,j])/n
        mu[j] <- sum(z.hat[,j]*qlogis(y))/sum(z.hat[,j])
        s[j] <- sqrt(sum(z.hat[,j]*(qlogis(y)-mu[j])^2)/sum(z.hat[,j]))
      }
      LL.new <- sum(log(dmixLN(y,w,mu,s))) # log-likelihood function
      count <- count +1 
      dif <- abs(LL.new/LL-1)
      LL <- LL.new
     cat('iter =', count, '\tloglike =', LL.new, '\n')
    } 
    aic <- -2 * LL.new + 2 * (2*g+g-1)
    bic <- -2 * LL.new + log(n) * (2*g+g-1)
    edc <- -2 * LL.new + 0.2*sqrt(n) * (2*g+g-1)
    end <- proc.time()[3]
    time <- end-begin
    obj.out <- list(family=family,w=w, mu=mu, sigma=s , loglik=LL.new, aic=aic, bic=bic, edc=edc, iter=count,elapsed=as.numeric(time))
  }
   if (group)
    obj.out$group <- apply(z.hat, 1, which.max)
  obj.out
}





