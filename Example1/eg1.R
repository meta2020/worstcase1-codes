## NOTES
## DATA IS SAVED IN "eg1-raw.csv"
## 

library(kableExtra)
library(ggplot2)
library(gridExtra)
library(tidyr)

## DATA PREPROCESS
example1 <- read.csv("eg1-raw.csv")
example1$sigmasq <- (1/example1$precision)^2

## META-ANALYSIS WITHOUT CONSIDERING PUBLICATION BIAS (PB)
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
b.MC.result <- read.csv("SASresult/result1-R10-K2000-15OCT202371943.csv")

## PLOT 1 
df1 <- b.MC.result[,c(1,2,4,5)]
df1$maxb <- as.vector(ma1$b)+df1$maxb
df1$bound <- as.vector(ma1$b)+df1$bound
df1$group <- as.factor(df1$group)
df_wide <- spread(df1, key = group, value = maxb)
names(df_wide)[-c(1,2)] <- paste0("sim", 1:10)

p1 <- ggplot(df_wide, aes(x=p)) + 
  geom_line(aes(y = sim1, color = "Simulation-based")) +
  geom_line(aes(y = sim2, color = "Simulation-based")) +
  geom_line(aes(y = sim3, color = "Simulation-based")) +
  geom_line(aes(y = sim4, color = "Simulation-based")) +
  geom_line(aes(y = sim5, color = "Simulation-based")) +
  geom_line(aes(y = sim6, color = "Simulation-based")) +
  geom_line(aes(y = sim7, color = "Simulation-based")) +
  geom_line(aes(y = sim8, color = "Simulation-based")) +
  geom_line(aes(y = sim9, color = "Simulation-based")) +
  geom_line(aes(y = sim10, color = "Simulation-based")) +
  geom_line(aes(y = bound, color = "Copas-Jackson"), lty=2, size=1) +
  scale_x_reverse(n.breaks = 10, name="Marginal selection probability") +
  scale_y_continuous(limits = c(-0.75,0.5), name = "Worst-case upper bound") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.2, 0.8), legend.background = element_rect(fill = "white", color = "black"))+
  scale_colour_manual("Worst-case bound",
                      breaks = c("Simulation-based", "Copas-Jackson"),
                      values = c("steelblue", "red"),
                      guide = guide_legend(override.aes = list(lty = c(1,2), size = c(0.5, 1)))) +
  ggtitle("A. 10 simulation-based and the Copas-Jackson bounds")

## PLOT 2 
df <- aggregate(b.MC.result[,c(4,5)], list(b.MC.result$p), FUN=median) 
b.MC <- rev(df$maxb)
names(b.MC) <- sprintf("p=%.2f", p)

b.MC.or <- as.vector(ma1$b)+b.MC
b.MC.ub <- b.MC.or+qnorm(0.975)*as.vector(ma1$se)
b.MC.lb <- b.MC.or-qnorm(0.975)*as.vector(ma1$se)

## PLOT 2
df <- data.frame(p = seq(1, 0.1, -0.1),
                 b.MC.or, b.MC.ub, b.MC.lb,
                 b.CJ.or, b.CJ.ub, b.CJ.lb)
p2 <- ggplot(df, aes(x = p)) +
  geom_ribbon(mapping=aes(ymin=b.CJ.ub,ymax=b.CJ.lb, fill="Copas-Jackson"), colour = "white", alpha=0.1)+
  geom_ribbon(mapping=aes(ymin=b.MC.ub,ymax=b.MC.lb, fill="Simulation-based"), colour = "white", alpha=0.1)+
  geom_line(aes(y = b.MC.or, colour="Simulation-based"), lty=1, size=1) +
  # geom_point(aes(y = b.MC.or), colour =2, pch=19, size=2) +
  geom_line(aes(y = b.CJ.or, colour="Copas-Jackson"), lty=2, size=1) +
  # geom_point(aes(y = b.CJ.or), colour =1, pch=1, size=2) +
  scale_x_reverse(n.breaks = 10, name="Marginal selection probability") +
  scale_y_continuous(limits = c(-0.75,0.5), name = "Worst-case upper bound") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.2, 0.8), legend.background = element_rect(fill = "white", color = "black"))+
  scale_colour_manual("Worst-case bound",
                      breaks = c("Simulation-based", "Copas-Jackson"),
                      values = c("steelblue", "red"),
                      guide = guide_legend(override.aes = list(lty = c(1,2))))+
  scale_fill_manual("95% CI region",
                      breaks = c("Simulation-based", "Copas-Jackson"),
                      values = c("steelblue", "red"),
                    guide = guide_legend(override.aes = list(color = c("white", "white")))) +
  ggtitle("B. Median of simulation-based and the Copas-Jackson bounds")

plot <- grid.arrange(p1, p2, ncol=2)
ggsave(filename = "eg1.eps", plot = plot, device = cairo_ps, width = 12, height = 6) 


## TABLE 1

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
