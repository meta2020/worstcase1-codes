## NOTES
## DATA IS SAVED IN "eg2-IVD.csv"
## 

library(kableExtra)

## DATA PREPROCESS
example2 <- read.csv("eg2-IVD.csv")

## CONTINUITY CORRECTION BY ADDING 0.5 TO ALL THE CELLS
example2.1 <- example2 
example2.1[,-1] <- example2[,-1]+0.5
write.csv(example2.1, "eg2-IVD-cc.csv", row.names = F)

## META-ANALYSIS USING THE REITSMA MODEL WITHOUT CONSIDERING PUBLICATION BIAS
library(mada)
ma2 <- reitsma(data=example2, method = "ml",correction.control = "all")

## THE ESTIMATES
par <- data.frame(
  mu1 = ma2$coefficients[1], mu2 = -ma2$coefficients[2],
  tau11 = ma2$Psi[1], tau22= ma2$Psi[4], tau12= -ma2$Psi[2])
write.csv(par, "ml-par.csv", row.names = F)

## SAUC
sroc <- function(x) plogis(par$mu1 - (par$tau12/par$tau22) * (qlogis(x) + par$mu2))
integrate(sroc, 0, 1)


## LOAD ESTIMATED BOUNDS
bound1 <- read.csv("MCbound4.1.csv")
bound2 <- read.csv("MCbound4.2.csv")
bound3 <- read.csv("MCbound4.3.csv")

p <- seq(1,0.1,-0.1)

## COMPARISON PLOT
setEPS(width = 12, height = 4); postscript("eg2.eps")

par(mfrow = c(1,3))
matplot(bound1[order(-bound1$P),2:3], ylab = "Worst-case bounds of the SAUC", xlab = "Overall selction probability (p)", ylim = c(0.5,1), xaxt = "n", type = "b", pch = 1, col = 1, lty = 1)
abline(h = 0.5, col = "grey", lty = 2)
axis(1, at = 1:10, labels = p)
title("(A) Constraint (4.1)", adj = 0, font.main = 1, cex.main = 1.5)
matplot(bound2[order(-bound2$P),2:3], ylab = "Worst-case bounds of the SAUC", xlab = "Overall selction probability (p)", ylim = c(0.5,1), xaxt = "n", type = "b", pch = 2, col = 1, lty = 1)
abline(h = 0.5, col = "grey", lty = 2)
axis(1, at = 1:10, labels = p)
title("(B) Constraint (4.2)", adj = 0, font.main = 1, cex.main = 1.5)
matplot(bound3[order(-bound3$P),2:3], ylab = "Worst-case bounds of the SAUC", xlab = "Overall selction probability (p)", ylim = c(0.5,1), xaxt = "n", type = "b", pch = 5, col = 1, lty = 1)
abline(h = 0.5, col = "grey", lty = 2)
axis(1, at = 1:10, labels = p)
title("(C) Constraint (4.3)", adj = 0, font.main = 1, cex.main = 1.5)
par(mfrow = c(1,1))

dev.off()

## TABLE2


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
