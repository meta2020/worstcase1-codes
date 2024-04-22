##
## EXAMPLE 2 IVD DATA
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
zhou.mod1 <- sapply(p.10, function(p) {
  opt2 <- try(dtametasa.rc(data=example2, p=p, c1.square0 = 0.5, beta.interval = c(0,2)))
  if(inherits(opt2, "try-error")) rep(NA,3) else c(opt2$sauc.ci)
})

zhou.mod2 <- sapply(p.10, function(p) {
  opt2 <- try(dtametasa.fc(data=example2, p=p, c1.square = 0.5, beta.interval = c(0,2)))
  if(inherits(opt2, "try-error")) rep(NA,3) else c(opt2$sauc.ci)  
})


zhou.mod3 <- sapply(p.10, function(p) {
  opt2 <- try(dtametasa.fc(data=example2, p=p, c1.square = 1, beta.interval = c(0,2)))
  if(inherits(opt2, "try-error")) rep(NA,3) else c(opt2$sauc.ci)  
})

zhou.mod4 <- sapply(p.10, function(p) {
  opt2 <- try(dtametasa.fc(data=example2, p=p, c1.square = 0, beta.interval = c(0,2)))
  if(inherits(opt2, "try-error")) rep(NA,3) else c(opt2$sauc.ci)  
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

save(zhou.mod1,zhou.mod2,zhou.mod3,zhou.mod4,df.piao,df.li,
     file = "comp-ivd-dt.RData")
