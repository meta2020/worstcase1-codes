##
## EXAMPLE 2 TROPNIN DATA
##

library(kableExtra)
library(mada)
# devtools::install_github("meta2020/dtametasa")
library(dtametasa)
library(tidyr)
source("data.pre.R")
source("sauc.R")
source("sauc.ci.R")
source("Piao2019.R")
source("Li2021.R")

## GENERATE ANALYTICAL DATA FOR SAS PROGRAMMING---------------------------------
## DATA PREPROCESS
example2 <- read.csv("../example-data/eg2-IVD.csv")
example2.1 <- correction(data = example2, type = "single") 
example2.2 <- logit.data(example2.1)
ma2 <- reitsma(data=example2, method = "ml", correction.control = "single")
par <- data.frame(
  mu1 = ma2$coefficients[1], mu2 = -ma2$coefficients[2],
  tau11 = ma2$Psi[1], tau22= ma2$Psi[4], tau12= -ma2$Psi[2])

## METHOD OF ZHOU ET AL
p.10 <- seq(1, 0.1, -0.1)

## SAUC
## c1=1,c2=0
zhou.sauc1 <- sapply(p.10, function(p) {
  opt <- try(dtametasa.fc(data=example2, p=p, c1.square = 1, beta.interval = c(0,2)))
  if(inherits(opt, "try-error")) rep(NA,3) else c(opt$sauc.ci)  
})

## c1=0,c2=1
zhou.sauc2 <- sapply(p.10, function(p) {
  opt <- try(dtametasa.fc(data=example2, p=p, c1.square = 0, beta.interval = c(0,2)))
  if(inherits(opt, "try-error")) rep(NA,3) else c(opt$sauc.ci)  
})

## c1=0.5,c2=0.5
zhou.sauc3 <- sapply(p.10, function(p) {
  opt <- try(dtametasa.fc(data=example2, p=p, c1.square = 0.5, beta.interval = c(0,2)))
  if(inherits(opt, "try-error")) rep(NA,3) else c(opt$sauc.ci)  
})

## ESTIMATE c1 c2
zhou.sauc4 <- sapply(p.10, function(p) {
  opt <- try(dtametasa.rc(data=example2, p=p, c1.square0 = 0.5, beta.interval = c(0,2)))
  if(inherits(opt, "try-error")) rep(NA,3) else c(opt$sauc.ci)
})

## PARAMETERS
## c1=1,c2=0
zhou.par1 <- sapply(p.10, function(p) {
  opt <- try(dtametasa.fc(data=example2, p=p, c1.square = 1, beta.interval = c(0,2)))
  if(inherits(opt, "try-error")) rep(NA,5) else c(opt$par.all[1:5])  
})

## c1=0,c2=1
zhou.par2 <- sapply(p.10, function(p) {
  opt <- try(dtametasa.fc(data=example2, p=p, c1.square = 0, beta.interval = c(0,2)))
  if(inherits(opt, "try-error")) rep(NA,5) else c(opt$par.all[1:5])  
})

## c1=0.5,c2=0.5
zhou.par3 <- sapply(p.10, function(p) {
  opt <- try(dtametasa.fc(data=example2, p=p, c1.square = 0.5, beta.interval = c(0,2)))
  if(inherits(opt, "try-error")) rep(NA,5) else c(opt$par.all[1:5])  
})

## ESTIMATE c1 c2
zhou.par4 <- sapply(p.10, function(p) {
  opt <- try(dtametasa.rc(data=example2, p=p, c1.square0 = 0.5, beta.interval = c(0,2)))
  if(inherits(opt, "try-error")) rep(NA,5) else c(opt$par.all[1:5])  
})



## METHOD OF PIAO ET AL
Piao.mod <- Piao2019(data = example2.2, init = c(unlist(par), -0.1, -0.1, rep(-0.1, 3)) )
Piao.sauc <- SAUC(Piao.mod$par[1:5])
Piao.p <- nrow(example2)/(nrow(example2)+ round(max(Piao.mod$m)))
df.piao <- data.frame(p = Piao.p, sauc = Piao.sauc)

## MEHTHOD OF LI ET AL
Li.mod <- Li2021(data = example2.2, eta1.tilde = Piao.mod$par[8:10], eta0 = Piao.mod$par)
Li.sauc <- SAUC(Li.mod$par[1:5])
Li.p <- Li.mod$alpha
df.li <- data.frame(p = Li.p, sauc = Li.sauc)

save(zhou.sauc1,zhou.sauc2,zhou.sauc3,zhou.sauc4,
     zhou.par1,zhou.par2,zhou.par3,zhou.par4,
     df.piao,df.li,
     file = "comp-ivd-dt.RData")
