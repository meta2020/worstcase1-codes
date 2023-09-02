##
## LOAD DATA
##
library(readxl)
files.sources <- list.files(path = "R/")
x <- sapply(paste0("R/", files.sources), source)


## MEDIAN FU 
med.data <- read_excel("Ki67.xlsx", sheet = "MCT")
med.data$mty <- med.data$mct_mo/12

## OBTAIN ETA
eta <- cens.eta(data = med.data, med.year = mty,  n1 = n1, n0 = n0, s1_med = s1_mct, s0_med = s0_mct)$par


## OS DATA

os.data  <- read_excel("Ki67.xlsx", sheet = "OS")
os.data$ty <- os.data$t/12

## HR DATA

dataHR <-read_excel("Ki67.xlsx", sheet = "HR")
dataHR$u_lnHR <- log(dataHR$HR)
ln_ci.low   <- log(dataHR$ci.low)
ln_ci.up    <- log(dataHR$ci.up)
se          <- (ln_ci.up - ln_ci.low)/(2 * qnorm(0.975))
dataHR$v_lnHR <- se^2

## MERGED DATA (OS+HR)
os_hr <- merge(os.data, dataHR, all.x = TRUE)

## GENERATE 3RD YEAR DATA
data3 <- convert.dt(
  data = os_hr, tK = 3, study = study, ty = ty, 
  n1 = n1, n0 = n0, s1 = s1, s0 = s0, u_lnHR = u_lnHR, v_lnHR = v_lnHR, 
  eta = eta)

## REMOVE OBSERVATIONS WITH MISSING VALUES
data3 <- na.omit(data3)

## GENERATE 5RD YEAR DATA
data5 <- convert.dt(
  data = os_hr, tK = 5, study = study, ty = ty, 
  n1 = n1, n0 = n0, s1 = s1, s0 = s0, u_lnHR = u_lnHR, v_lnHR = v_lnHR, 
  eta = eta)

## REMOVE OBSERVATIONS WITH MISSING VALUES
data5 <- na.omit(data5)

data <- rbind.data.frame(data3,data5)

write.csv(data, "data-35y.csv", row.names = F)

## META-ANALYSIS USING THE HZ MODEL WITHOUT CONSIDERING PUBLICATION BIAS

data <- data3
y1  <- data$u_sen
y2  <- data$u_spe 
v1  <- data$v_sen
v2  <- data$v_spe 
v12 <- data$v_senspe 

ma3.llk  <- function(par) llk.BNM.ml(par, y1, y2, v1, v2, v12)
ma3 <- nlminb(c(rep(0.1,4), -0.1), ma3.llk)
par3 <- ma3$par
## THE ESTIMATES
par3 <- data.frame(
  t = 3,
  theta1 = par3[1], theta2 = par3[2],
  tau11 = par3[3]^2, tau22= par3[4]^2, tau12= prod(par3[3:5]))
write.csv(par3, "ml-par3.csv", row.names = F)

data <- data5
y1  <- data$u_sen
y2  <- data$u_spe 
v1  <- data$v_sen
v2  <- data$v_spe 
v12 <- data$v_senspe 

ma3.llk  <- function(par) llk.BNM.ml(par, y1, y2, v1, v2, v12)
ma3 <- nlminb(c(rep(0.1,4), -0.1), ma3.llk)

par3 <- ma3$par
## THE ESTIMATES
par3 <- data.frame(
  t = 5,
  theta1 = par3[1], theta2 = par3[2],
  tau11 = par3[3]^2, tau22= par3[4]^2, tau12= prod(par3[3:5]))
write.csv(par3, "ml-par5.csv", row.names = F)




## LOAD ESTIMATED BOUNDS
bound1 <- read.csv("MCbound3-5.1.csv")
bound2 <- read.csv("MCbound3-5.2.csv")
bound3 <- read.csv("MCbound3-5.3.csv")
bound4 <- read.csv("MCbound5-5.1.csv")
bound5 <- read.csv("MCbound5-5.2.csv")
bound6 <- read.csv("MCbound5-5.3.csv")

p <- seq(1,0.1,-0.1)

## COMPARISON PLOT
setEPS(width = 12, height = 8); postscript("eg3.eps")

par(mfrow = c(2,3))
matplot(bound1[order(-bound1$P),2:3], ylab = "Worst-case bounds of the SAUC(3)", xlab = "Overall selction probability (p)", ylim = c(0.4,0.9), xaxt = "n", type = "b", pch = 1, col = 1, lty = 1)
abline(h = 0.5, col = "grey", lty = 2)
axis(1, at = 1:10, labels = p)
title("(A) Constraint (5.1)", adj = 0, font.main = 1, cex.main = 1.5)
matplot(bound2[order(-bound2$P),2:3], ylab = "Worst-case bounds of the SAUC(3)", xlab = "Overall selction probability (p)", ylim = c(0.4,0.9), xaxt = "n", type = "b", pch = 2, col = 1, lty = 1)
abline(h = 0.5, col = "grey", lty = 2)
axis(1, at = 1:10, labels = p)
title("(B) Constraint (5.2)", adj = 0, font.main = 1, cex.main = 1.5)
matplot(bound3[order(-bound3$P),2:3], ylab = "Worst-case bounds of the SAUC(3)", xlab = "Overall selction probability (p)", ylim = c(0.4,0.9), xaxt = "n", type = "b", pch = 5, col = 1, lty = 1)
abline(h = 0.5, col = "grey", lty = 2)
axis(1, at = 1:10, labels = p)
title("(C) Constraint (5.3)", adj = 0, font.main = 1, cex.main = 1.5)
matplot(bound4[order(-bound1$P),2:3], ylab = "Worst-case bounds of the SAUC(5)", xlab = "Overall selction probability (p)", ylim = c(0.4,0.9), xaxt = "n", type = "b", pch = 1, col = 1, lty = 1)
abline(h = 0.5, col = "grey", lty = 2)
axis(1, at = 1:10, labels = p)
title("(D) Constraint (5.1)", adj = 0, font.main = 1, cex.main = 1.5)
matplot(bound5[order(-bound2$P),2:3], ylab = "Worst-case bounds of the SAUC(5)", xlab = "Overall selction probability (p)", ylim = c(0.4,0.9), xaxt = "n", type = "b", pch = 2, col = 1, lty = 1)
abline(h = 0.5, col = "grey", lty = 2)
axis(1, at = 1:10, labels = p)
title("(E) Constraint (5.2)", adj = 0, font.main = 1, cex.main = 1.5)
matplot(bound6[order(-bound3$P),2:3], ylab = "Worst-case bounds of the SAUC(5)", xlab = "Overall selction probability (p)", ylim = c(0.4,0.9), xaxt = "n", type = "b", pch = 5, col = 1, lty = 1)
abline(h = 0.5, col = "grey", lty = 2)
axis(1, at = 1:10, labels = p)
title("(F) Constraint (5.3)", adj = 0, font.main = 1, cex.main = 1.5)
par(mfrow = c(1,1))

dev.off()

## TABLE2

p <- seq(1,0.1,-0.1)
tab1 <- cbind.data.frame(p = p, 
                         A = sprintf("[%.3f, %.3f]", bound1[order(-bound1$P),2], bound1[order(-bound1$P),3]), 
                         B = sprintf("[%.3f, %.3f]", bound2[order(-bound2$P),2], bound2[order(-bound2$P),3]),
                         C = sprintf("[%.3f, %.3f]", bound3[order(-bound3$P),2], bound3[order(-bound3$P),3]))
colnames(tab1) <- c("p", "Constraint 4.1", "Constraint 4.2", "Constraint 4.3")

sink("tab2.tex")
kbl(tab1, 
    format = "latex",
    longtable = F, 
    booktabs = T, 
    linesep = "\\hline",
    digits = 3,
    align = "r",
    escape = FALSE,
    caption = "The values of the MC bound and the CJ bound in Example 1",
    label = "tab1",
    row.names = NA)
sink()
