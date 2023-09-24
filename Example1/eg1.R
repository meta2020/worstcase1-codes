## NOTES
## DATA IS SAVED IN "eg1-raw.csv"
## 

library(kableExtra)
library(ggplot2)
## DATA PREPROCESS
example1 <- read.csv("eg1-raw.csv")
example1$sigmasq <- (1/example1$precision)^2

## META-ANALYSIS WITHOUT CONSIDERING PUBLICATION BIAS
library(metafor)
ma1 <- rma.uni(yi = y, vi = sigmasq, data = example1, method = "ML")
ma1$b
ma1$se
ma1$ci.ub
## ANALYTICAL COPAS-JACKSON BOUND (tau=0)
p <- seq(1, 0.1, -0.1)
u <- qnorm(p)
b.CJ <- with(example1, {1/p*dnorm(u)*sum(precision)/sum(precision^2)})
names(b.CJ) <- sprintf("p=%.2f", p)

b.CJ.or <- as.vector(ma1$b)+b.CJ
b.CJ.ub <- b.CJ.or+qnorm(0.975)*as.vector(ma1$se)
b.CJ.lb <- b.CJ.or-qnorm(0.975)*as.vector(ma1$se)

## IMPORT THE MC BOUND
# b.MC.result <- read.csv("SASresult/result1-R10-K2000-19SEP202352272.csv")
b.MC.result <- read.csv("SASresult/result1-R1-K2000-19SEP202353299.csv")

df <- aggregate(b.MC.result[,c(4,5,6,9,10)], list(b.MC.result$p), FUN=mean) 

# b.MC <- rev(b.MC.result$admaxb)
b.MC <- rev(df$maxb2)
names(b.MC) <- sprintf("p=%.2f", p)

b.MC
b.MC.or <- as.vector(ma1$b)+b.MC
b.MC.ub <- b.MC.or+qnorm(0.975)*as.vector(ma1$se)
b.MC.lb <- b.MC.or-qnorm(0.975)*as.vector(ma1$se)

## PLOT 1
df <- data.frame(p = seq(1, 0.1, -0.1),
                 b.MC.or, b.MC.ub, b.MC.lb,
                 b.CJ.or, b.CJ.ub, b.CJ.lb)
plot1 <- ggplot(df, aes(x = p)) + 
  geom_ribbon(mapping=aes(ymin=b.MC.ub,ymax=b.MC.lb, fill="Simulation-based"), colour="white", alpha=0.1)+
  geom_ribbon(mapping=aes(ymin=b.CJ.ub,ymax=b.CJ.lb, fill="Copas-Jackson"), colour="white", alpha=0.1)+
  geom_line(aes(y = b.MC.or, colour="Simulation-based")) +
  geom_point(aes(y = b.MC.or), colour =2, pch=19, size=2) + 
  geom_line(aes(y = b.CJ.or, colour="Copas-Jackson"), lty=2) +
  geom_point(aes(y = b.CJ.or), colour =1, pch=1, size=2) + 
  scale_x_reverse(n.breaks = 10, name="Overall selection probability") + 
  scale_y_continuous(limits = c(-0.75,0.5), name = "Worst-case upper bound") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"))+
  scale_colour_manual("The upper bound", 
                      breaks = c("Simulation-based", "Copas-Jackson"),
                      values = c("red", "black"),
                      guide = guide_legend(override.aes = list(lty = c(1,2))))+
  scale_fill_manual("95% CI region", 
                      breaks = c("Simulation-based", "Copas-Jackson"),
                      values = c("red", "black"),
                    guide = guide_legend(override.aes = list(color = c("white", "white"))))
plot1
ggsave(filename = "eg1.eps", plot = plot1, device = cairo_ps, 
       width = 8, height = 6) 

# setEPS(width = 6, height = 6); postscript("eg1.eps")
# plot(b.MC.or, ylab = "Maximum bias (b)", xlab = "Overall selction probability (p)", ylim = c(-0.7,0.5), xaxt = "n", type = "b")
# lines(b.MC.ub, lty=2)
# lines(b.MC.lb, lty=2)
# points(b.CJ.or, pch = 2, type = "b")
# lines(b.CJ.ub)
# lines(b.CJ.lb)
# axis(1, at = 1:10, labels = p)
# legend("bottomright", legend = c("The MC bound", "The CJ bound"), pch=c(1,2), cex = 1.2)
# dev.off()

## TABLE1

tab1 <- cbind.data.frame(p = p, sprintf("%.3f [%.3f, %.3f]",b.MC.or, b.MC.lb, b.MC.ub), 
                         sprintf("%.3f [%.3f, %.3f]",b.CJ.or, b.CJ.lb, b.CJ.ub))
colnames(tab1) <- c("p", "The MC bound [95% CI]", "The CJ bound [95% CI]")
sink("tab1.tex")
kbl(tab1, 
    format = "latex",
    longtable = F, 
    booktabs = T, 
    linesep = "",
    digits = 3,
    align = "r",
    escape = FALSE,
    caption = "The values of the MC bound and the CJ bound in Example 1",
    label = "tab1",
    row.names = NA)
sink()
