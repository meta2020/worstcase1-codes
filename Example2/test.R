curve((1-dnorm(x))*dnorm(x), -qnorm(0.5), 10, xlim = c(-10,10), ylim=c(0,1))
curve(dnorm(x)*dnorm(x), -10, -qnorm(0.5), add = TRUE)
abline(v=-qnorm(0.5))

integrate(function(x) {(1-dnorm(x))*dnorm(x)},-Inf, -qnorm(0.7))$value

integrate({function(x) dnorm(x)*dnorm(x)}, -qnorm(0.7),Inf)$value


library(dtametasa)
fit <- dtametasa.fc(data=example2, p=1, correct.type = "all", c1.square = 1)
ldata <- fit$l.data
mlpar <- fit$par
Omega <- matrix(c(mlpar[3]^2, prod(mlpar[3:5]), prod(mlpar[3:5]), mlpar[4]^2),2,2)

Sigma1.list <- lapply(1:33, function(i){
  
  Sigmai <- matrix(c(ldata$v1[i], 0, 0, ldata$v2[i]), 2,2) + Omega
  mpower(Sigmai, -1)
}) 

Sigma2.list <- lapply(1:33, function(i){
  
  Sigmai <- matrix(c(ldata$v1[i], 0, 0, ldata$v2[i]), 2,2) + Omega
  mpower(Sigmai, -1/2)
}) 

E1 <- Reduce("+", Sigma1.list) 
E2 <- Reduce("+", Sigma2.list) 

p <- seq(1,0.1,-0.1)
C <- matrix(c(1,1), 2,1)
C <- matrix(c(0,1), 2,1)
as.vector(t(C) %*% solve(E1) %*% E2 %*% C)*(dnorm(qnorm(p))/p)
