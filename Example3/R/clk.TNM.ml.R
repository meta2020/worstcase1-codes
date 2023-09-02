##******************************************************************************
##
## LIKELIHOOD OF THE OBSERVED (NEGATIVE)
##
##******************************************************************************

clk.TNM.ml <- function(
  par, y1, y2, y3,
  v1, v2, v3, v12, v13, v23,
  p, a.interval
  
){

  u1  <- par[1]      ## logit-se
  u2  <- par[2]      ## logit-sp
  u3  <- par[3]      ## ln-HR
  
  t1  <- par[4]      ## logit-se
  t2  <- par[5]    
  t3  <- par[6]
  
  r1  <- par[7]      ## logit-sp
  r2  <- par[8]
  r3  <- par[9]      ## ln-HR
  
  b   <- par[10] 
  # b <- 1

  t11 <- t1^2
  t22 <- t2^2
  t33 <- t3^2

  t12 <- t1*t2*r1
  t23 <- t2*t3*r2
  t13 <- t1*t3*r3
  
  v11 <- v1  + t11
  v22 <- v2  + t22
  v33 <- v3  + t33
  
  v12 <- v12 + t12 
  v13 <- v13 + t13
  v23 <- v23 + t23

  t_lnHR   <- abs(y3/sqrt(v3))

  ##
  ## FUNCTOIN b(Sigma) ----
  ##

  f.b <- function(a){
    
    sq <- suppressWarnings(sqrt(1 + b^2 * (1 + t33 / v3)))
    
    pnorm( (a + b * u3/sqrt(v3)) / sq )
    
  }


  ##
  ## FIND THE ROOT OF a = a.opt ----
  ##
  
  n <- length(y1)

  a.p <- function(a) {sum(1/f.b(a), na.rm = TRUE) - n/p}

  a.opt.try <- suppressWarnings(try(uniroot(a.p, a.interval, extendInt="no"), silent = TRUE)) 

  a.opt <- a.opt.try$root


  ##
  ##  LOGLIKELIHOOD-1 OF y|Sigma ----
  ##
  
  det1  <- v22 * v33 - v23 * v23
  det2  <- v12 * v33 - v13 * v23
  det3  <- v12 * v23 - v13 * v22
  
  det   <- v11 * det1 - v12 * det2 + v13 * det3
  
  log.dev <- suppressWarnings(log(det)) 
  
  a11 <-  (v22 * v33 - v23 * v23)
  a12 <- -(v12 * v33 - v13 * v23)
  a13 <-  (v12 * v23 - v13 * v22)
  
  a21 <- -(v12 * v33 - v23 * v13)
  a22 <-  (v11 * v33 - v13 * v13)
  a23 <- -(v11 * v23 - v12 * v13)
  
  a31 <-  (v12 * v23 - v13 * v22)
  a32 <- -(v11 * v23 - v12 * v13)
  a33 <-  (v11 * v22 - v12 * v12)
  
  Y1  <- y1 - u1
  Y2  <- y2 - u2
  Y3  <- y3 - u3
  
  y_Adj_y <- Y1 * (Y1 * a11 + Y2 * a21 + Y3 * a31) + 
    Y2 * (Y1 * a12 + Y2 * a22 + Y3 * a32) +
    Y3 * (Y1 * a13 + Y2 * a23 + Y3 * a33)
  
  
  
  s.l1 <- -0.5 * (sum(log.dev + y_Adj_y / det, na.rm = TRUE))
  

  ##
  ##  LOGLIKELIHOOD-2 OF a(a.opt) ----
  ##


  f.l2 <- pnorm(a.opt + b * t_lnHR)

  s.l2 <- sum( log(f.l2), na.rm = TRUE )


  ##
  ##  LOGLIKELIHOOD-3 OF b(a.opt) ----
  ##

  f.l3 <- f.b(a.opt)

  s.l3 <- sum( log(f.l3), na.rm = TRUE )


  ##
  ##  FINAL LOGLIKELIHOOD ----
  ##

  -(s.l1 + s.l2 - s.l3) ## NEGATIVE

  
}
