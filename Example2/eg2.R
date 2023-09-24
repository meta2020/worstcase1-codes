## NOTES
## DATA IS SAVED IN "eg2-IVD.csv"
## 

library(kableExtra)
library(mada)
# devtools::install_github("meta2020/dtametasa")
library(dtametasa)
library(gridExtra)
library(ggplot2)

## DATA PREPROCESS
example2 <- read.csv("eg2-IVD.csv")

## CONTINUITY CORRECTION BY ADDING 0.5 TO ALL THE CELLS
example2.1 <- example2 
example2.1[,-1] <- example2[,-1]+0.5
write.csv(example2.1, "eg2-IVD-cc.csv", row.names = F)

## META-ANALYSIS USING THE REITSMA MODEL WITHOUT CONSIDERING PUBLICATION BIAS

ma2 <- reitsma(data=example2, method = "ml",correction.control = "all")


## THE ESTIMATES
par <- data.frame(
  mu1 = ma2$coefficients[1], mu2 = -ma2$coefficients[2],
  tau11 = ma2$Psi[1], tau22= ma2$Psi[4], tau12= -ma2$Psi[2])
write.csv(par, "ml-par.csv", row.names = F)

## SOP
mu1 <- ma2$coefficients[1]
mu2 <- -ma2$coefficients[2]
se <- plogis(mu1)
sp <- plogis(mu2)

## SROC


## BOUNDS UNDER CONDITION D1-D3
y1d1.dat <- read.csv("SASresult/result2-R10-K2000-c1-as1-23SEP202351067.csv")
y1d2.dat <- read.csv("SASresult/result2-R10-K2000-c1-as2-23SEP202353513.csv")
y1d3.dat <- read.csv("SASresult/result2-R10-K2000-c1-as3-23SEP202354471.csv")

y2d1.dat <- read.csv("SASresult/result2-R10-K2000-c2-as1-23SEP202356785.csv")
y2d2.dat <- read.csv("SASresult/result2-R10-K2000-c2-as2-23SEP202357572.csv")
y2d3.dat <- read.csv("SASresult/result2-R10-K2000-c2-as3-23SEP202382616.csv")

y3d1.dat <- read.csv("SASresult/result2-R10-K2000-c3-as1-23SEP202360218.csv")
y3d2.dat <- read.csv("SASresult/result2-R10-K2000-c3-as2-23SEP202375429.csv")
y3d3.dat <- read.csv("SASresult/result2-R10-K2000-c3-as3-23SEP202362479.csv")


# y1d1.dat <- read.csv("SASresult/result2-R1-K2000-c1-as1-24SEP202368307.csv")
# y1d2.dat <- read.csv("SASresult/result2-R1-K2000-c1-as2-24SEP202368518.csv")
# y1d3.dat <- read.csv("SASresult/result2-R1-K2000-c1-as3-24SEP202368919.csv")
# 
# y2d1.dat <- read.csv("SASresult/result2-R1-K2000-c2-as1-24SEP202369179.csv")
# y2d2.dat <- read.csv("SASresult/result2-R1-K2000-c2-as2-24SEP202369560.csv")
# y2d3.dat <- read.csv("SASresult/result2-R1-K2000-c2-as3-24SEP202369651.csv")
# 
# y3d1.dat <- read.csv("SASresult/result2-R1-K2000-c3-as1-24SEP202369757.csv")
# y3d2.dat <- read.csv("SASresult/result2-R1-K2000-c3-as2-24SEP202370013.csv")
# y3d3.dat <- read.csv("SASresult/result2-R1-K2000-c3-as3-24SEP202370822.csv")

## the average
y1d1 <- aggregate(y1d1.dat[,c("maxb", "minb")], list(p=y1d1.dat$p), FUN=mean) 
y1d2 <- aggregate(y1d2.dat[,c("maxb", "minb")], list(p=y1d2.dat$p), FUN=mean) 
y1d3 <- aggregate(y1d3.dat[,c("maxb", "minb")], list(p=y1d3.dat$p), FUN=mean) 

y2d1 <- aggregate(y2d1.dat[,c("maxb", "minb")], list(p=y2d1.dat$p), FUN=mean) 
y2d2 <- aggregate(y2d2.dat[,c("maxb", "minb")], list(p=y2d2.dat$p), FUN=mean) 
y2d3 <- aggregate(y2d3.dat[,c("maxb", "minb")], list(p=y2d3.dat$p), FUN=mean) 

y3d1 <- aggregate(y3d1.dat[,c("maxb", "minb")], list(p=y3d1.dat$p), FUN=mean) 
y3d2 <- aggregate(y3d2.dat[,c("maxb", "minb")], list(p=y3d2.dat$p), FUN=mean) 
y3d3 <- aggregate(y3d3.dat[,c("maxb", "minb")], list(p=y3d3.dat$p), FUN=mean) 

## POINTS OF STUDIES

sei <- with(example2.1, TP/(TP+FN))
spi <- with(example2.1, TN/(TN+FP))

## SOP WITH BIAS
sed1 <- plogis(mu1+y1d1[,"minb"])
spd1 <- plogis(mu2+y2d1[,"minb"])
sed2 <- plogis(mu1+y1d2[,"minb"])
spd2 <- plogis(mu2+y2d2[,"minb"])
sed3 <- plogis(mu1+y1d3[,"minb"])
spd3 <- plogis(mu2+y2d3[,"minb"])



df <- data.frame(y = sei, x = 1-spi)

sroc.func <- function(x, par, bias){
  
  u1  <- par[1]
  u2  <- par[2]
  t12 <- par[5]
  t22 <- par[4]
  
  plogis(u1 - (t12/t22) * (qlogis(x) + u2) + bias)
  
}

sroc.func2 <- function(x, par, b1, b2){
  
  u1  <- par[1]
  u2  <- par[2]
  t12 <- par[5]
  t22 <- par[4]
  
  plogis(u1+b1 - (t12/t22) * (qlogis(x)+b2 + u2) )
  
}


## D1
srocd1 <- ggplot(df, aes(x, y)) + geom_point(shape = "x") + 
  scale_x_continuous(limits = c(0,1), name = "1 - Specificity") +
  scale_y_continuous(limits = c(0,1), name = "Sensitivity") +
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1[y3d1$p==1,"minb"]), aes(color ="p=1"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1[y3d1$p==0.9,"minb"]), aes(color ="p=0.9"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1[y3d1$p==0.8,"minb"]), aes(color ="p=0.8"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1[y3d1$p==0.7,"minb"]), aes(color ="p=0.7"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1[y3d1$p==0.6,"minb"]), aes(color ="p=0.6"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1[y3d1$p==0.5,"minb"]), aes(color ="p=0.5"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1[y3d1$p==0.4,"minb"]), aes(color ="p=0.4"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1[y3d1$p==0.3,"minb"]), aes(color ="p=0.3"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1[y3d1$p==0.2,"minb"]), aes(color ="p=0.2"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1[y3d1$p==0.1,"minb"]), aes(color ="p=0.1"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position="none")+
  ggtitle("A. SROC curves under Condition (D1)") +
  scale_color_brewer(palette = "Spectral")

## D2 
srocd2 <- ggplot(df, aes(x, y)) + geom_point(shape = "x") + 
  scale_x_continuous(limits = c(0,1), name = "1 - Specificity") +
  scale_y_continuous(limits = c(0,1), name = "Sensitivity") +
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d2[y3d2$p==1,"minb"]), aes(color ="p=1"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d2[y3d2$p==0.9,"minb"]), aes(color ="p=0.9"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d2[y3d2$p==0.8,"minb"]), aes(color ="p=0.8"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d2[y3d2$p==0.7,"minb"]), aes(color ="p=0.7"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d2[y3d2$p==0.6,"minb"]), aes(color ="p=0.6"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d2[y3d2$p==0.5,"minb"]), aes(color ="p=0.5"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d2[y3d2$p==0.4,"minb"]), aes(color ="p=0.4"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d2[y3d2$p==0.3,"minb"]), aes(color ="p=0.3"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d2[y3d2$p==0.2,"minb"]), aes(color ="p=0.2"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d2[y3d2$p==0.1,"minb"]), aes(color ="p=0.1"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position="none")+
  ggtitle("B. SROC curves under Condition (D2)") +
  scale_color_brewer(palette = "Spectral")

## D3
srocd3 <- ggplot(df, aes(x, y)) + geom_point(shape = "x") + 
  scale_x_continuous(limits = c(0,1), name = "1 - Specificity") +
  scale_y_continuous(limits = c(0,1), name = "Sensitivity") +
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d3[y3d3$p==1,"minb"]), aes(color ="p=1"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d3[y3d3$p==0.9,"minb"]), aes(color ="p=0.9"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d3[y3d3$p==0.8,"minb"]), aes(color ="p=0.8"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d3[y3d3$p==0.7,"minb"]), aes(color ="p=0.7"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d3[y3d3$p==0.6,"minb"]), aes(color ="p=0.6"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d3[y3d3$p==0.5,"minb"]), aes(color ="p=0.5"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d3[y3d3$p==0.4,"minb"]), aes(color ="p=0.4"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d3[y3d3$p==0.3,"minb"]), aes(color ="p=0.3"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d3[y3d3$p==0.2,"minb"]), aes(color ="p=0.2"))+
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d3[y3d3$p==0.1,"minb"]), aes(color ="p=0.1"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position="none")+
  ggtitle("C. SROC curves under Condition (D3)") +
  scale_color_brewer(palette = "Spectral")


pplot <- ggplot(df, aes(x, y)) + geom_point(shape = "x") + 
  scale_x_continuous(limits = c(0,1), name = "1 - Specificity") +
  scale_y_continuous(limits = c(0,1), name = "Sensitivity") +
  geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1[y3d1$p==1,"minb"]), aes(color ="p=1"))+
  geom_line(data = data.frame(x = 1-spd1, y = sed1), color = "black", alpha =0.5)+ 
  annotate(geom="text", x=0.5, y=0.6, label="Condition (D1)", color="black")+
  geom_point(data = data.frame(x = 1-spd1, y = sed1), aes(color = paste0("p=",seq(0.1, 1, 0.1))))+
  geom_line(data = data.frame(x = 1-spd2, y = sed2), color = "red", alpha =0.5)+
  annotate(geom="text", x=0.5, y=0.25, label="Condition (D2)", color="red")+
  geom_point(data = data.frame(x = 1-spd2, y = sed2), aes(color = paste0("p=",seq(0.1, 1, 0.1))))+
  geom_line(data = data.frame(x = 1-spd3, y = sed3), color = "green", alpha =0.5)+
  annotate(geom="text", x=0.5, y=0.4, label="Condition (D3)", color="green")+
  geom_point(data = data.frame(x = 1-spd3, y = sed3), aes(color = paste0("p=",seq(0.1, 1, 0.1))))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"))+
  labs(color=' ', title = "D. SOPs under 3 conditions")+
  scale_color_brewer(palette = "Spectral")

p1 <- grid.arrange(srocd1, srocd2, srocd3, pplot, ncol=2)
ggsave(filename = "eg2-1-rep10.eps", plot = p1, device = cairo_ps, width = 12, height = 12) 

# ma21<- dtametasa.fc(data=example2, p=1, correct.type = "all")
# l1 <- ma21$sauc.ci[2]-ma21$sauc.ci[1]
# l2 <- ma21$sauc.ci[3]-ma21$sauc.ci[1]

## Under CONDITION D1-D3
## PLOT 2

# df1 <- data.frame(p=p, bound = rev(bound1$SAUC_LB), lb = rev(bound1$SAUC_LB+l1), ub = rev(bound1$SAUC_LB+l2))
# df2 <- data.frame(p=p, bound = rev(bound2$SAUC_LB), lb = rev(bound2$SAUC_LB+l1), ub = rev(bound2$SAUC_LB+l2))
# df3 <- data.frame(p=p, bound = rev(bound3$SAUC_LB), lb = rev(bound3$SAUC_LB+l1), ub = rev(bound3$SAUC_LB+l2))
# 
# plot1 <- ggplot(df1, aes(x = p)) + 
#   geom_ribbon(mapping=aes(ymin=lb,ymax=ub), fill="black", colour="white", alpha=0.1)+
#   geom_line(aes(y = bound), colour="black") +
#   geom_point(aes(y = bound), colour ="black", pch=19, size=2) + 
#   scale_x_reverse(n.breaks = 10, name="Overall selection probability") + 
#   scale_y_continuous(limits = c(0.45,1), n.breaks = 7, name = "Worst-case lower bound") +
#   theme(panel.background = element_rect(fill = "white", colour = "grey50"),
#         panel.grid.major = element_line(colour = "grey87"))+
#   ggtitle("(A) Condition (D1)") 
# 
# plot2 <- ggplot(df2, aes(x = p)) + 
#   geom_ribbon(mapping=aes(ymin=lb,ymax=ub), fill="black", colour="white", alpha=0.1)+
#   geom_line(aes(y = bound), colour="black") +
#   geom_point(aes(y = bound), colour ="black", pch=19, size=2) + 
#   scale_x_reverse(n.breaks = 10, name="Overall selection probability") + 
#   scale_y_continuous(limits = c(0.45,1), n.breaks = 7, name = "Worst-case lower bound") +
#   theme(panel.background = element_rect(fill = "white", colour = "grey50"),
#         panel.grid.major = element_line(colour = "grey87"))+
#   ggtitle("(B) Condition (D2)") 
# 
# plot3 <- ggplot(df3, aes(x = p)) + 
#   geom_ribbon(mapping=aes(ymin=lb,ymax=ub), fill="black", colour="white", alpha=0.1)+
#   geom_line(aes(y = bound), colour="black") +
#   geom_point(aes(y = bound), colour ="black", pch=19, size=2) + 
#   scale_x_reverse(n.breaks = 10, name="Overall selection probability") + 
#   scale_y_continuous(limits = c(0.45,1), n.breaks = 7, name = "Worst-case lower bound") +
#   theme(panel.background = element_rect(fill = "white", colour = "grey50"),
#         panel.grid.major = element_line(colour = "grey87"))+
#   ggtitle("(C) Condition (D3)") 
# plot <- grid.arrange(plot1, plot2, ncol=2)
# ggsave(filename = "eg2.eps", plot = plot, device = cairo_ps, 
#        width = 12, height = 4) 


## SAUC

sauc.value <- function(par, bias){
  
  sauc <- integrate({function (x) sroc.func(x, par, bias)}, 0, 1)
  sauc$value
  
}

##
fit <- dtametasa.fc(example2.1, p = 1)
fit$sauc.ci
saucl <- fit$sauc.ci[1] - fit$sauc.ci[2]
saucu <- fit$sauc.ci[3] - fit$sauc.ci[1]

##
saucd1 <- sapply(10:1, function(i) sauc.value(par= unlist(par), bias=y3d1[i,"minb"]) )
saucd1.lb <- saucd1-saucl
saucd1.ub <- saucd1+saucu
saucd2 <- sapply(10:1, function(i) sauc.value(par= unlist(par), bias=y3d2[i,"minb"]) )
saucd2.lb <- saucd2-saucl
saucd2.ub <- saucd2+saucu
saucd3 <- sapply(10:1, function(i) sauc.value(par= unlist(par), bias=y3d3[i,"minb"]) )
saucd3.lb <- saucd3-saucl
saucd3.ub <- saucd3+saucu



## SAUC PLOT
df <- data.frame(p = seq(1, 0.1, -0.1),
                 saucd1, saucd1.lb, saucd1.ub,
                 saucd2, saucd2.lb, saucd2.ub)

plot1 <- ggplot(df, aes(x = p)) + 
  geom_ribbon(mapping=aes(ymin=saucd1.lb,ymax=saucd1.ub, fill="D1/D3"), colour="white", alpha=0.1)+
  geom_ribbon(mapping=aes(ymin=saucd2.lb,ymax=saucd2.ub, fill="D2"), colour="white", alpha=0.1)+
  geom_line( aes(y = saucd1, colour="D1/D3")) +
  geom_point(aes(y = saucd1, colour="D1/D3"), pch=19, size=2) + 
  geom_line( aes(y = saucd2, colour="D2"), lty=2) +
  geom_point(aes(y = saucd2, colour="D2"), pch=19, size=2)+ 
  scale_x_reverse(n.breaks = 10, name="Overall selection probability") + 
  scale_y_continuous(limits = c(0,1), name = "Worst-case lower bound") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"))+
  labs(color='Conditions', fill='95% CI region', title = "")+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")
ggsave(filename = "eg2-rep10.eps", plot = plot1, device = cairo_ps, width = 8, height = 6)


## TABLE2


tab1 <- cbind.data.frame(p = p, 
                         A = sprintf("%.3f [%.3f, %.3f]", df[,2], df[,3], df[,4]), 
                         B = sprintf("%.3f [%.3f, %.3f]", df[,5], df[,6], df[,7]))
colnames(tab1) <- c("p", "Condition (D1)/(D3)", "Condition (D2)")

sink("tab2.tex")
kbl(tab1, 
    format = "latex",
    longtable = F, 
    booktabs = T, 
    linesep = "",
    digits = 3,
    align = "r",
    escape = FALSE,
    caption = "The values of the MC bound and the Zhou et al. bound in Example 1",
    label = "tab1",
    row.names = NA)
sink()
