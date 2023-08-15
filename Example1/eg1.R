## NOTES
## DATA IS SAVED IN "eg1-raw.csv"
## 


## DATA PREPROCESS
example1 <- read.csv("eg1-raw.csv")
example1$sigmasq <- (1/example1$precision)^2

## META-ANALYSIS WITHOUT CONSIDERING PUBLICATION BIAS
library(metafor)
rma.uni(yi = y, vi = sigmasq, data = example1, method = "ML")

## ANALYTICAL COPAS-JACKSON BOUND (tau=0)
p <- seq(1, 0.1, -0.1)
u <- qnorm(p)
b.CJ <- with(example1, {1/p*dnorm(u)*sum(precision)/sum(precision^2)})
names(b.CJ) <- sprintf("p=%.2f", p)
b.CJ

## IMPORT THE MC BOUND
b.MC.result <- read.csv("MCbounds.csv")
b.MC <- rev(b.MC.result$maxb)
names(b.MC) <- sprintf("p=%.2f", p)

## COMPARISON PLOT
setEPS(width = 6, height = 6); postscript("eg1.eps")
plot(b.MC, ylab = "Maximum bias (b)", xlab = "Overall selction probability (p)", ylim = c(0,0.7), xaxt = "n")
points(b.CJ, pch = 2)
axis(1, at = 1:10, labels = p)
legend("bottomright", legend = c("The MC bound", "The CJ bound"), 
       pch=c(1,2), cex = 1.2)
dev.off()