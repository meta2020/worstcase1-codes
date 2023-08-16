## NOTES
## DATA IS SAVED IN "eg1-raw.csv"
## 

library(kableExtra)

## DATA PREPROCESS
example1 <- read.csv("eg1-raw.csv")
example1$sigmasq <- (1/example1$precision)^2

## META-ANALYSIS WITHOUT CONSIDERING PUBLICATION BIAS
library(metafor)
ma1 <- rma.uni(yi = y, vi = sigmasq, data = example1, method = "ML")
ma1$b
ma1$se
ma1$ci.ub+b.MC
## ANALYTICAL COPAS-JACKSON BOUND (tau=0)
p <- seq(1, 0.1, -0.1)
u <- qnorm(p)
b.CJ <- with(example1, {1/p*dnorm(u)*sum(precision)/sum(precision^2)})
names(b.CJ) <- sprintf("p=%.2f", p)

b.CJ
b.CJ.or <- as.vector(ma1$b)+b.CJ
b.CJ.ub <- b.CJ.or+2*pnorm(0.975)*as.vector(ma1$se)
b.CJ.lb <- b.CJ.or-2*pnorm(0.975)*as.vector(ma1$se)

## IMPORT THE MC BOUND
b.MC.result <- read.csv("MCbounds.csv")
b.MC <- rev(b.MC.result$maxb)
names(b.MC) <- sprintf("p=%.2f", p)

b.MC
b.MC.or <- as.vector(ma1$b)+b.MC
b.MC.ub <- b.MC.or+2*pnorm(0.975)*as.vector(ma1$se)
b.MC.lb <- b.MC.or-2*pnorm(0.975)*as.vector(ma1$se)

## COMPARISON PLOT
setEPS(width = 6, height = 6); postscript("eg1.eps")

plot(b.MC, ylab = "Maximum bias (b)", xlab = "Overall selction probability (p)", ylim = c(0,0.7), xaxt = "n", type = "b")
points(b.CJ, pch = 2, type = "b")
axis(1, at = 1:10, labels = p)
legend("bottomright", legend = c("The MC bound", "The CJ bound"), 
       pch=c(1,2), cex = 1.2)
dev.off()

## TABLE1

tab1 <- cbind.data.frame(p = p, sprintf("%.3f [%.3f, %.3f]",b.MC.or, b.MC.lb, b.MC.ub), 
                         sprintf("%.3f [%.3f, %.3f]",b.CJ.or, b.CJ.lb, b.CJ.ub))
colnames(tab1) <- c("p", "The MC bound [95% CI]", "The MC bound [95% CI]")
sink("tab1.tex")
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
