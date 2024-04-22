.DID.sroc <- function(u1, u2, t1, t2, r, var.matrix){
  
  Q1 <- function(x) {
    
    g <- plogis(u1 - (t1*t2*r/(t2^2)) * (qlogis(x) + u2))
    g*(1-g)
  }
  
  Q2 <- function(x) {
    
    g <- plogis(u1 - (t1*t2*r/(t2^2)) * (qlogis(x) + u2))
    p.u2 <- (-r*t1/t2)*g*(1-g)
  }
  
  Q3 <- function(x) {
    
    g <- plogis(u1 - (t1*t2*r/(t2^2)) * (qlogis(x) + u2))
    p.t1 <- (-r/t2*(qlogis(x)+u2))*g*(1-g)
    
  }
  
  Q4 <- function(x) {
    
    g <- plogis(u1 - (t1*t2*r/(t2^2)) * (qlogis(x) + u2))
    p.t2 <- r*t1/t2^2*( qlogis(x)+u2)*g*(1-g)
    
  }
  
  Q5 <- function(x) {
    
    g <- plogis(u1 - (t1*t2*r/(t2^2)) * (qlogis(x) + u2))
    p.r  <- (-t1)/t2*(qlogis(x) + u2)*g*(1-g)
  }
  
  fd <- c(integrate(Q1, 0, 1)$value,
          integrate(Q2, 0, 1)$value,
          integrate(Q3, 0, 1)$value,
          integrate(Q4, 0, 1)$value,
          integrate(Q5, 0, 1)$value
  )
  
  (fd %*% var.matrix %*% fd)
  
}


saucci <- function(par, var.matrix, sauc){
  
  u1 <- par$mu1
  u2 <- par$mu2
  t1 <- sqrt(par$tau11)
  t2 <- sqrt(par$tau22)
  r  <- par$tau12/(t1*t2)
  # matv <- fit$var.ml
  saucv<- .DID.sroc(u1, u2, t1, t2, r, var.matrix)
  lb <- plogis(qlogis(sauc)+qnorm((1-0.95)/2)*sqrt(saucv)/(sauc*(1-sauc)))
  ub <- plogis(qlogis(sauc)+qnorm((1-0.95)/2, lower.tail = F)*sqrt(saucv)/(sauc*(1-sauc)))
  
  return(c(lb, ub))
  
}






