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

## IMPORT THE MC BOUND; NEED RESULTS FROM SAS**!!-------------------------------
b.MC.result <- read.csv("SAS/NEQ/result1-K2000-R1.csv")
##

## GENERATE PLOTS AND TABLES----------------------------------------------------
## 
## 
df <- b.MC.result[,c(1,2,3,7)]
df$theta <- rep(ma1$b,10)
df$theta.lb <- rep(ma1$ci.lb,10)
df$theta.ub <- rep(ma1$ci.ub,10)

df$MC.ub <- as.vector(ma1$b)+df$maxb
df$MC.lb <- as.vector(ma1$b)+df$minb
# df$MC.cub <- as.vector(ma1$b)+df$maxb+qnorm(0.975)*as.vector(ma1$se)
# df$MC.clb <- as.vector(ma1$b)+df$minb-qnorm(0.975)*as.vector(ma1$se)

df$CJ.ub <- as.vector(ma1$b)+b.CJ.max
df$CJ.lb <- as.vector(ma1$b)+b.CJ.min
# df$CJ.cub <- as.vector(ma1$b)+b.CJ.ub
# df$CJ.clb <- as.vector(ma1$b)+b.CJ.lb


## PLOT S1 ----
## 
lgd2 <- c("Estimate without PB",
         "Copas-Jackson bound",
         "Simulation-based bound")
dfb <- data.frame(
  p = rep(df$p, 3),
  est1 = c(df$theta, df$CJ.ub, df$MC.ub),
  est2 = c(df$theta, df$CJ.lb, df$MC.lb),
  grp = rep(lgd2, each=10)
)
dfl <- data.frame(
  p = rep(df$p, 2),
  est = c(df$CJ.lb, df$MC.lb),
  grp = rep(lgd2[-1], each=10)
)
p2 <- ggplot(
  dfb, aes(x = p, y = est1, group = grp, color = grp, linetype=grp)) +
  geom_line(size=1) +
  geom_line(aes(x = p, y = est2, group = grp, color = grp, linetype=grp), linewidth=1)+
  scale_y_continuous(limits = c(-1.5,0.5), name = "lnOR") +
  scale_x_reverse(n.breaks = 10, name="Marginal selection probability") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.2, 0.2), 
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank(),
        legend.text=element_text(size=15),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))+
  scale_colour_manual("",
                      breaks = lgd2,
                      values = c("grey50","steelblue","red","steelblue","red"))+
  scale_linetype_manual("",
                        breaks = lgd2,
                        values = c(2,1,2,2,3))

ggsave(filename = "fig-tab/eg1-bothbounds.eps", plot = p2, device = cairo_ps, width = 8, height = 8)

## TABLE 2: UPPER AND LOWER BOUNDS ----
## 
tab1 <- cbind.data.frame(
  p = p,
  m = round(nrow(example1)*(1/p-1)),
  sprintf("%.3f", df$MC.lb),
  sprintf("%.3f", df$MC.ub),
  sprintf("%.3f", df$CJ.lb),
  sprintf("%.3f", df$CJ.ub))
colnames(tab1) <- c("$p$", "$M$","Lower bound", "Upper bound", "Lower bound", "Upper bound")

kbl(tab1, 
    format = "latex",
    longtable = F, 
    booktabs = T, 
    linesep = "",
    digits = 3,
    align = "r",
    escape = FALSE,
    caption = "Example 1: worst-case bounds of the lnOR.",
    label = "tab1",
    position = "h",
    row.names = NA) %>%
  add_header_above(c("","", "Simulation-based bounds"=2, "Copas-Jackson bounds"=2), escape = FALSE)%>% 
  footnote(general = "
           $M$ indicates the number of unpublished studies.", 
           escape = FALSE, threeparttable = T,  general_title = "")



## PLOT 3A: COMPARISON OF SETTING DIFFERENT K (UPPER BOUNDS, ONE TIME) ----
##
b.MC.result1 <- read.csv("SAS/NEQ/result1-K1000-R1.csv")
b.MC.result2 <- read.csv("SAS/NEQ/result1-K2000-R1.csv")
b.MC.result3 <- read.csv("SAS/NEQ/result1-K5000-R1.csv")
b.MC.result4 <- read.csv("SAS/NEQ/result1-K20000-R1.csv")

ldg3 <- c("Copas-Jackson bound","K=1000","K=2000","K=5000","K=20000")
dfk <- data.frame(
  p = rep(df$p,5),
  est = c(df$CJ.ub, as.vector(ma1$b)+c(b.MC.result1$maxb,b.MC.result2$maxb,b.MC.result3$maxb,b.MC.result4$maxb)),
  grp = rep(ldg3,each = 10)
)

p3a <- ggplot(
  dfk, aes(x = p, y = est, group = grp, color = grp, linetype=grp)) +
  geom_line(size=1) +
  scale_y_continuous(limits = c(-0.5,0.27), name = "lnOR") +
  scale_x_reverse(n.breaks = 10, name="Marginal selection probability") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.25, 0.8), 
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank(),
        legend.text=element_text(size=15))+
  scale_colour_manual("",
                      breaks = ldg3,
                      values = c("#377eb8","#4daf4a","#e41a1c","#984ea3","#ff7f00"))+
  scale_linetype_manual("",
                        breaks = ldg3,
                        values = c(1,2,2,2,2))+
  ggtitle("A. One-time simulation-based upper bounds")

## TABLE 3A: COMPARISON OF SETTING DIFFERENT K (UPPER BOUNDS, ONE TIME) ----
##

tab21 <- dfk %>% spread(key = grp, value = est) %>% as.data.frame()
tab21.p <- tab21[order(tab21$p, decreasing = T),]
rownames(tab21.p) <- c()

kbl(tab21.p[,c(1,2,3,4,6,5)], 
    format = "latex",
    longtable = F, 
    booktabs = T, 
    linesep = "",
    digits = 3,
    align = "r",
    escape = FALSE,
    caption = "Example 1: the upper bounds of the lnOR by the simulation-based given different $K$.",
    label = "tab1",
    row.names = NA) %>%
  add_header_above(c("", "", "Simulation-based upper bounds"=4), escape = FALSE)  

## PLOT 3B: COMPARISON OF SETTING DIFFERENT K AND 10 SAMPLES ----------------------
##

b.MC.rep1 <- read.csv("SAS/NEQ/result1-K1000-R10.csv")
b.MC.rep2 <- read.csv("SAS/NEQ/result1-K2000-R10.csv")
b.MC.rep3 <- read.csv("SAS/NEQ/result1-K5000-R10.csv")
b.MC.rep4 <- read.csv("SAS/NEQ/result1-K20000-R10.csv")

lgd4 <- c("K=1000","K=2000","K=5000","K=20000") 
bpdf0 <- rbind(b.MC.rep1,b.MC.rep2,b.MC.rep3,b.MC.rep4)
bpdf <- data.frame(
  pp = as.factor(c(bpdf0$p)),
  est = bpdf0$maxb+as.vector(ma1$b),
  grp = rep(lgd4,each=100)
  ) 
bpdf$grp2 <- factor(bpdf$grp, levels = lgd4, ordered = TRUE)
p3b <- ggplot(bpdf, aes(x=pp, y=est, fill=grp2)) + 
  geom_boxplot(alpha = 0.5)+
  scale_x_discrete(limits = rev(levels(bpdf$pp)),name="Marginal selection probability") +
  scale_y_continuous(limits = c(-0.5,0.27), name = "lnOR") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.2, 0.8), legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank(),
        legend.text=element_text(size=15),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(size = 13))+
  scale_fill_manual("",
                      breaks = lgd4,
                      values = c("#4daf4a","#e41a1c","#984ea3","#ff7f00")) +
  ggtitle("B. Boxplot of 10-time simulation-based upper bounds")

p3 <- grid.arrange(p3a, p3b, ncol=2)
ggsave(filename = "fig-tab/eg1-k.eps", plot = p3, device = cairo_ps, width = 12, height = 6)

