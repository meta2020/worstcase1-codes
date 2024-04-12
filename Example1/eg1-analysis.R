##
## EXAMPLE 1
##


library(kableExtra)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(metafor)

## NEED TO SET LOCAL PATH-------------------------------------------------------
# setwd()

## DATA PREPROCESS
example1 <- read.csv("eg1-data.csv")
example1$si2 <- (1/example1$precision)^2

## RANDOM-EFFECTS META-ANALYSIS WITHOUT CONSIDERING PUBLICATION BIAS (PB)
ma1 <- rma.uni(yi = y, vi = si2, data = example1, method = "ML")

## PRINT RESULTS WITHOUT PB
sprintf("mu (SE; 95.CI) = %.3f (%.3f; %.3f to %.3f); tau2 (SE) = %.3f (%.3f)", 
        ma1$b, ma1$se, ma1$ci.lb, ma1$ci.ub, ma1$tau2, ma1$se.tau2)

## CALCULATE ANALYTICAL COPAS-JACKSON BOUND (tau=0)
p <- seq(1, 0.1, -0.1)
u <- qnorm(p)
b.CJ.max <- with(example1, {1/p*dnorm(u)*sum(precision)/sum(precision^2)})
b.CJ.min <- with(example1, -{1/p*dnorm(u)*sum(precision)/sum(precision^2)})
names(b.CJ.max) <- sprintf("p=%.2f", p)
names(b.CJ.min) <- sprintf("p=%.2f", p)

## COPAS-JACKSONG BOUNG OF PB
## MAX BIAS FOR THE UPPER BOUND
b.CJ.max
## MIN BIAS FOR THE LOWER BOUND
b.CJ.min

## ANALYTICAL COPAS-JACKSON BOUND (tau=0) 95%CI BASED ON SE WITHOUT PB
b.CJ.ub <- b.CJ.max+qnorm(0.975)*as.vector(ma1$se)
b.CJ.lb <- b.CJ.min-qnorm(0.975)*as.vector(ma1$se)

## CONFIDENCE LOWER AND UPPER BANDS
rbind(b.CJ.ub, b.CJ.lb)


##------------------------------------------------------------------------------
## IMPORT THE MC BOUND; NEED RESULTS FROM SAS**!!-------------------------------
b.MC.result <- read.csv("SAS/result1-2000.csv")
##------------------------------------------------------------------------------


##------------------------------------------------------------------------------
## GENERATE PLOTS AND TABLES----------------------------------------------------
## 
## PLOT 1
##  
df <- b.MC.result[,c(1,2,4)]

df$theta <- rep(as.vector(ma1$b), 10)
df$theta.lb <- rep(as.vector(ma1$ci.lb), 10)
df$theta.ub <- rep(as.vector(ma1$ci.ub), 10)

df$theta.max <- as.vector(ma1$b)+df$maxb
df$theta.min <- as.vector(ma1$b)+df$minb
df$theta.max.ub <- as.vector(ma1$b)+df$maxb+qnorm(0.975)*as.vector(ma1$se)
df$theta.min.lb <- as.vector(ma1$b)+df$minb-qnorm(0.975)*as.vector(ma1$se)

df$theta.CJ.max <- as.vector(ma1$b)+b.CJ.max
df$theta.CJ.min <- as.vector(ma1$b)+b.CJ.min
df$theta.CJ.max.ub <- as.vector(ma1$b)+b.CJ.ub
df$theta.CJ.min.lb <- as.vector(ma1$b)+b.CJ.lb


p1 <- ggplot(df, aes(x=p)) + 
  geom_ribbon(mapping=aes(ymin=theta.lb,ymax=theta.ub, color="Estimation without PB"), fill = "grey50", alpha=0.1, lty=3)+
  geom_ribbon(mapping=aes(ymin=theta.CJ.max,ymax=theta.CJ.max.ub, color="Copas-Jackson worst-case bound"), fill = "steelblue", alpha=0.1, lty=3)+
  geom_ribbon(mapping=aes(ymin=theta.CJ.min,ymax=theta.CJ.min.lb, color="Copas-Jackson worst-case bound"), fill = "steelblue", alpha=0.1, lty=3)+
  geom_ribbon(mapping=aes(ymin=theta.max,ymax=theta.max.ub, color="Simulation-based worst-case bound"), fill = "red", alpha=0.1, lty=3)+
  geom_ribbon(mapping=aes(ymin=theta.min,ymax=theta.min.lb, color="Simulation-based worst-case bound"), fill = "red", alpha=0.1, lty=3)+
  geom_line(aes(y = theta, color = "Estimation without PB"), lty=2, size=1) +
  geom_line(aes(y = theta.CJ.max, color = "Copas-Jackson worst-case bound"), lty=1, size=1) +
  geom_line(aes(y = theta.CJ.min, color = "Copas-Jackson worst-case bound"), lty=1, size=1) +
  geom_line(aes(y = theta.max, color = "Simulation-based worst-case bound"), lty=2, size=1) +
  geom_line(aes(y = theta.min, color = "Simulation-based worst-case bound"), lty=2, size=1) +
  scale_x_reverse(n.breaks = 10, name="Marginal selection probability") +
  scale_y_continuous(limits = c(-2,0.5), name = "lnOR") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.2, 0.2), 
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank())+
  scale_colour_manual("",
                      breaks = c("Simulation-based worst-case bound",
                                 "Copas-Jackson worst-case bound",
                                 "Estimation without PB"),
                      values = c("red","steelblue", "grey50"),
                      guide = guide_legend(override.aes = list(lty = c(1,2,2), 
                                                               fill = c("red", "steelblue", "grey50"))))

p1

ggsave(filename = "fig-tab/eg1.eps", plot = p1, device = cairo_ps, width = 8, height = 8)

## TABLE 1

tab1 <- cbind.data.frame(p = p, 
                         sprintf("%.3f (%.3f)", df$theta.min, df$theta.min.lb),
                         sprintf("%.3f (%.3f)", df$theta.max, df$theta.max.ub),
                         sprintf("%.3f (%.3f)", df$theta.CJ.min, df$theta.CJ.min.lb),
                         sprintf("%.3f (%.3f)", df$theta.CJ.max, df$theta.CJ.max.ub))
colnames(tab1) <- c("p", "LB (95\\% CLB)", "UB (95\\% CUB)", "LB (95\\% CLB)", "UB (95\\% CUB)")
sink("fig-tab/tab1.tex")
kbl(tab1, 
    format = "latex",
    longtable = F, 
    booktabs = T, 
    linesep = "",
    digits = 3,
    align = "r",
    escape = FALSE,
    caption = "Example 1: the upper and lower bounds of the lnOR by the simulation-based and the Copas-Jackson bounds.",
    label = "tab1",
    row.names = NA) %>%
  add_header_above(c("", "Simulation-based bounds"=2, "Copas-Jackson bounds"=2), escape = FALSE) %>% 
  footnote(general = "
           LB indicates lower bound; CLB indicates the confidence lower band of the LB; 
           UB indicates upper bound; CUB indicates the confidence upper band of the UB.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")
sink()


## PLOT 2: COMPARISON OF SETTING DIFFERENT K
b.MC.result1 <- read.csv("SASresult/result1-1000.csv")
b.MC.result2 <- read.csv("SASresult/result1-2000.csv")
b.MC.result3 <- read.csv("SASresult/result1-20000.csv")

df <- b.MC.result1[,c(1), drop = F]
df$theta.max1 <- as.vector(ma1$b)+b.MC.result1$maxb
df$theta.min1 <- as.vector(ma1$b)+b.MC.result1$minb
df$theta.max2 <- as.vector(ma1$b)+b.MC.result2$maxb
df$theta.min2 <- as.vector(ma1$b)+b.MC.result2$minb
df$theta.max3 <- as.vector(ma1$b)+b.MC.result3$maxb
df$theta.min3 <- as.vector(ma1$b)+b.MC.result3$minb
df$theta.CJ.max <- as.vector(ma1$b)+b.CJ.max
df$theta.CJ.min <- as.vector(ma1$b)+b.CJ.min
df$theta <- rep(as.vector(ma1$b), 10)

p2 <- ggplot(df, aes(x=p)) + 
  geom_line(aes(y = theta, color = "Estimation without PB"), lty=2, size=1) +
  geom_line(aes(y = theta.CJ.max, color = "Copas-Jackson worst-case bound"), lty=1, size=1) +
  geom_line(aes(y = theta.CJ.min, color = "Copas-Jackson worst-case bound"), lty=1, size=1) +
  geom_line(aes(y = theta.max1, color = "Simulation-based worst-case bound (K=1000)"), lty=1, size=1) +
  geom_line(aes(y = theta.min1, color = "Simulation-based worst-case bound (K=1000)"), lty=1, size=1) +
  geom_line(aes(y = theta.max2, color = "Simulation-based worst-case bound (K=2000)"), lty=2, size=1) +
  geom_line(aes(y = theta.min2, color = "Simulation-based worst-case bound (K=2000)"), lty=2, size=1) +
  geom_line(aes(y = theta.max3, color = "Simulation-based worst-case bound (K=20000)"), lty=3, size=1) +
  geom_line(aes(y = theta.min3, color = "Simulation-based worst-case bound (K=20000)"), lty=3, size=1) +
  scale_x_reverse(n.breaks = 10, name="Marginal selection probability") +
  scale_y_continuous(limits = c(-2,0.5), name = "lnOR") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.3, 0.2), 
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank())+
  scale_colour_manual("",
                      breaks = c("Simulation-based worst-case bound (K=1000)",
                                 "Simulation-based worst-case bound (K=2000)",
                                 "Simulation-based worst-case bound (K=20000)",
                                 "Copas-Jackson worst-case bound",
                                 "Estimation without PB"),
                      values = c("#808800","#dc143c","#60100b","steelblue","grey50"),
                      guide = guide_legend(override.aes = list(lty = c(1,2,3,1,2))))
p2
ggsave(filename = "fig-tab/eg1-k.eps", plot = p2, device = cairo_ps, width = 8, height = 8)



