##------------------------------------------------------------------------------
##
## HECKMANTYPE SELECTION FUNCTION IN Li et al. (2021)
##
##------------------------------------------------------------------------------


h1 <- function(N, alpha, n){
  
  h1 <- log(choose(N, n)) + (N - n)*log(1 - alpha)
  
  return(h1)
}

h2 <- function(eta, y1, y2, v1, v2){
  
  beta1  <- eta[1]; beta2  <- eta[2]
  tau1   <- eta[3]; tau2   <- eta[4]
  rho.b  <- eta[5]
  rho1   <- eta[6]; rho2   <- eta[7];
  gamma1 <- eta[8]; gamma2 <- eta[9]; gamma3 <- eta[10]
  
  
  h2i <- sapply(1:length(y1), function(i){
    
    Sigma.i12 <- c(rho1 * sqrt(v1[i]), rho2 * sqrt(v2[i]))
    Sigma.i22 <- matrix(c(tau1^2 + v1[i], rho.b * tau1 * tau2, rho.b * tau1 * tau2, tau2^2 + v2[i]), 2, 2)
    D.i       <- c(y1[i] - beta1, y2[i] - beta2) 
    
    inv.Sigma.i22 <- solve(Sigma.i22)
    
    numer <- gamma1 + gamma2/sqrt(v1[i]) + gamma3/sqrt(v2[i]) + Sigma.i12 %*% inv.Sigma.i22 %*% D.i
    denom <- sqrt(1 - Sigma.i12 %*% inv.Sigma.i22 %*% Sigma.i12)
    
    v <- numer/denom
    v[v < -37] <- -37
    v[v > 37]  <- 37
    term1 <- log(pnorm(v))
    term2 <- -0.5*log(det(Sigma.i22))
    term3 <- -0.5*t(D.i) %*% inv.Sigma.i22 %*% D.i
    
    term1 + term2 + term3
    
  })
  
  h2 <- sum(h2i, na.rm = TRUE)
  
  return(h2)
}


h3 <- function(lambda, alpha, eta1, n, v1, v2){
  
  gamma1 <- eta1[1]; gamma2 <- eta1[2]; gamma3 <- eta1[3]
  probit <- pnorm(gamma1 + gamma2/sqrt(v1) + gamma3/sqrt(v2))
  z      <- 1 + lambda*(probit - alpha)
  
  logz   <- rep(0, length(z)) 
  logz[z>1/n] <- log(z[z>1/n])
  logz[z<= 1/n] <- -log(n) - 1.5 + 2*z[z<= 1/n]*n - 0.5*z[z<= 1/n]^2*n^2
  
  h3 <- -sum(logz, na.rm = TRUE)
  
  return(h3)
}


h23.new <- function(eta, y1, y2, v1, v2){
  
  beta1  <- eta[1]; beta2  <- eta[2]
  tau1   <- eta[3]; tau2   <- eta[4]
  rho.b  <- eta[5]
  rho1   <- eta[6]; rho2   <- eta[7];
  gamma1 <- eta[8]; gamma2 <- eta[9]; gamma3 <- eta[10]
  
  
  h2i <- sapply(1:length(y1), function(i){
    
    Sigma.i12 <- c(rho1 * sqrt(v1[i]), rho2 * sqrt(v2[i]))
    Sigma.i22 <- matrix(c(tau1^2 + v1[i], rho.b * tau1 * tau2, rho.b * tau1 * tau2, tau2^2 + v2[i]), 2, 2)
    D.i       <- c(y1[i] - beta1, y2[i] - beta2) 
    
    inv.Sigma.i22 <- solve(Sigma.i22)
    
    numer <- gamma1 + gamma2/sqrt(v1[i]) + gamma3/sqrt(v2[i]) + Sigma.i12 %*% inv.Sigma.i22 %*% D.i
    denom <- sqrt(1 - Sigma.i12 %*% inv.Sigma.i22 %*% Sigma.i12)
    
    v <- numer/denom
    v[v < -37] <- -37
    v[v > 37]  <- 37
    term1 <- log(pnorm(v))
    term2 <- -0.5*log(det(Sigma.i22))
    term3 <- -0.5*t(D.i) %*% inv.Sigma.i22 %*% D.i
    
    term1 + term2 + term3
    
  })
  
  h2 <- sum(h2i, na.rm = TRUE)
  
  ## h3
  
  probit <- pnorm(gamma1 + gamma2/sqrt(v1) + gamma3/sqrt(v2))
  z      <- 1 + lambda*(probit - alpha)
  
  logz   <- rep(0, length(z)) 
  logz[z>1/n] <- log(z[z>1/n])
  logz[z<= 1/n] <- -log(n) - 1.5 + 2*z[z<= 1/n]*n - 0.5*z[z<= 1/n]^2*n^2
  
  h3 <- -sum(logz, na.rm = TRUE)
  
  return(-(h2 + h3))
}


Li2021 <- function(
    data,
    eta1.tilde, eta0){
  
## optimize alpha
##   

dy1 <- data$y1
dy2 <- data$y2
dv1 <- data$v1
dv2 <- data$v2
dn  <- nrow(data)

eta1.piao <- eta1.tilde
eta.init <- eta0
names(eta.init) <- c("mu1", "mu2", "tau1", "tau2", "rho",
                     "rho1", "rho2", "gamma0", "gamma1", "gamma2")

h123 <- function(alpha.given) {
  
  
  h1.opt <- optimize(h1, c(dn, 100), maximum = TRUE, alpha = alpha.given, n = dn)
  h1.tilde <- h1.opt$maximum
  
  
  h3.opt <- optimize(h3, c(-10, 10), maximum = FALSE, alpha = alpha.given, eta1 = eta1.piao, n = dn, v1 = dv1, v2 = dv2)
  h3.tilde <- h3.opt$minimum
  
  h23 <- function(eta, alpha = alpha.given) {
    -(h2(eta, y1 = dy1, y2 = dy2, v1 = dv1, v2 = dv2) + h3.tilde)
    }  # h23 = -h2 - h3
  
  
  h23.opt <- nlminb(
    eta.init,
    h23,
    lower=c(-5, -5, 0, 0, -1, -1, -1, -5, 0, 0), 
    upper=c( 5,  5, 3, 3,  1,  1,  1,  5, 5, 5))
  h23.tilde <- h23.opt$objective
  
  return(-(h1.tilde + h23.tilde))
  
} # - h1 - h2 - h3


alpha.opt <- optimize(h123, c(0, 1), maximum = TRUE)
alpha.hat <- alpha.opt$maximum

## step 4
# optimize N and eta

N.hat.opt <- optimize(h1, c(dn, 100), maximum = TRUE, alpha = alpha.hat, n = dn)
N.hat <- round(N.hat.opt$maximum)

h3.opt <- optimize(h3, c(-10, 10), maximum = FALSE, alpha = alpha.hat, eta1 = eta1.piao, n = dn, v1 = dv1, v2 = dv2)
h3.tilde <- h3.opt$minimum

h23 <- function(eta, alpha = alpha.hat) {
  -(h2(eta, y1 = dy1, y2 = dy2, v1 = dv1, v2 = dv2) + h3.tilde)
}

h23.opt <- nlminb(
  eta.init,
  h23,
  lower=c(-5, -5, 0, 0, -1, -1, -1, -5, 0, 0),
  upper=c( 5,  5, 3, 3,  1,  1,  1,  5, 5, 5))


res <- list(par = h23.opt$par,
            alpha = alpha.hat,
            N = N.hat,
            conv = h23.opt$convergence
  
)
return(res)
}


