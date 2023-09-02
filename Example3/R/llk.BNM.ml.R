##******************************************************************************
##
## LIKELIHOOD OF THE OBSERVED (NEGATIVE)
##
##******************************************************************************


llk.BNM.ml <- function(par,y1, y2, v1, v2, v12){

  u1 <- par[1]
  u2 <- par[2]
  
  t1 <- par[3]
  t2 <- par[4]
  
  r  <- par[5]

  t11 <- t1^2
  t22 <- t2^2
  t12 <- t1*t2*r
  
  v11 <- v1  + t11
  v22 <- v2  + t22
  v12 <- v12 + t12


  ##
  ##  LOGLIKELIHOOD-1 OF y|Sigma ----
  ##

  det.vec <- v11 * v22 -v12^2

  log.det.vec <- suppressWarnings(log(det.vec))

  f.l1  <- ((y1-u1)^2*(v22) - 2*(y2-u2)*(y1-u1)*v12 + (y2-u2)^2*(v11)) / det.vec + log.det.vec

  s.l1  <- 0.5*sum(f.l1, na.rm = TRUE)

  return(s.l1)
}
