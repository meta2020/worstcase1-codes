##******************************************************************************
##
## DATA PRE-PROCESS
##
##******************************************************************************

##
## CONTINUITY CORRECTION -------------------------------------------------------
##

correction <- function(
  data, 
  value = 0.5,
  type = c("single", "all")
){

  type <- match.arg(type)

  if(type == "all"){

    if(any(c(data$TP,data$FN,data$FP,data$TN) == 0)){

      data$TP <- data$TP + value
      data$FN <- data$FN + value
      data$FP <- data$FP + value
      data$TN <- data$TN + value
    }
  }

  if(type == "single"){

    correction = ((((data$TP == 0)|(data$FN == 0))|(data$FP == 0))| (data$TN == 0)) * value

    data$TP <- correction + data$TP
    data$FN <- correction + data$FN
    data$FP <- correction + data$FP
    data$TN <- correction + data$TN

  }

  return(data)

}


##
## TRANSFORM DATA: TO GENERATE y1 y2 -------------------------------------------
##

logit.data <- function(data){

  sens <- data$TP/(data$TP+data$FN)
  spec <- data$TN/(data$TN+data$FP)

  y1 <- qlogis(sens)  ##y1 <- log(sens/(1-sens))
  y2 <- qlogis(spec)  ##y2 <- log(spec/(1-spec))

  v1 <- (1/data$TP)+(1/data$FN)
  v2 <- (1/data$TN)+(1/data$FP)

  data.frame(sens = sens,
             spec = spec,
             y1 = y1,
             y2 = y2,
             v1 = v1,
             v2 = v2)

}


sim.pdata <- function(par){
  
  S  <- par[1]
  u1 <- par[2]
  u2 <- par[3]
  
  t11 <- par[4]
  t22 <- par[5]
  t12 <- par[6]
  
  v1 <- rnorm(S, 0.1, 0.5)^2
  v2 <- rnorm(S, 0.5, 0.5)^2
  
  u   <- c(u1, u2)
  
  ## Sigma+Omega
  
  SO <- lapply(1:S, function(s){
    
    matrix(c(v1[s]+t11, t12, t12, v2[s]+t22),2,2)
    
  })
  
  ## CHECK PD
  
  check.PD <- sapply(1:S, function(s){
    
    eigen(SO[[s]])$values
    
  })
  
  if (any(as.vector(check.PD)<= 0)) stop("VAR-COV MATRIX (S+O) is NOT PD")
  
  ## y FROM N(u, SO)
  
  Y <- t(sapply(1:S, function(s) mvtnorm::rmvnorm(1, u, SO[[s]])))
  
  ## SENS AND SPEC
  
  X <- plogis(Y)
  
  DT <- cbind.data.frame(Y, v1, v2, X)
  
  colnames(DT) <- c("y1", "y2", "v1", "v2", "se", "sp")
  
  DT
}
