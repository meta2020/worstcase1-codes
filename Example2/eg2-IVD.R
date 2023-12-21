##
## Example 2 Troponin data
##

library(kableExtra)
library(mada)
# devtools::install_github("meta2020/dtametasa")
library(dtametasa)
library(gridExtra)
library(ggplot2)
library(tidyr)
source("data.pre.R")
source("plotfunc.R")
source("sauc.R")
source("sauc.ci.R")
source("sroc.R")

## DATA PREPROCESS
example2 <- read.csv("eg2-IVD.csv")

## CONTINUITY CORRECTION BY ADDING 0.5 TO THE STUDIES WITH 0 CELLS
example2.1 <- correction(data = example2, type = "single") 
write.csv(example2.1, "eg2-IVD-cc.csv", row.names = F)

## META-ANALYSIS USING THE REITSMA MODEL WITHOUT CONSIDERING PUBLICATION BIAS

ma2 <- reitsma(data=example2, method = "ml", correction.control = "single")

## THE ESTIMATES
par <- data.frame(
  mu1 = ma2$coefficients[1], mu2 = -ma2$coefficients[2],
  tau11 = ma2$Psi[1], tau22= ma2$Psi[4], tau12= -ma2$Psi[2])
write.csv(par, "ml-par-IVD.csv", row.names = F)




## THE SOP
mu1 <- ma2$coefficients[1]
mu2 <- -ma2$coefficients[2]
se <- plogis(mu1)
sp <- plogis(mu2)


## WORST-CASE BOUNDS UNDER CONDITION D1-D3, REPEAT 10 TIMES, PRODUCED BY SAS
## WHEN c = (1, 0)
y1d1.dat <- read.csv("SASresult/ivd/result2-R10-K2000-c1-as1.csv")
y1d2.dat <- read.csv("SASresult/ivd/result2-R10-K2000-c1-as2.csv")
y1d3.dat <- read.csv("SASresult/ivd/result2-R10-K2000-c1-as3.csv")
## WHEN c = (0, 1)
y2d1.dat <- read.csv("SASresult/ivd/result2-R10-K2000-c2-as1.csv")
y2d2.dat <- read.csv("SASresult/ivd/result2-R10-K2000-c2-as2.csv")
y2d3.dat <- read.csv("SASresult/ivd/result2-R10-K2000-c2-as3.csv")
## WHEN tilde c IS USED
y3d1.dat <- read.csv("SASresult/ivd/result2-R10-K2000-c3-as1.csv")
y3d2.dat <- read.csv("SASresult/ivd/result2-R10-K2000-c3-as2.csv")
y3d3.dat <- read.csv("SASresult/ivd/result2-R10-K2000-c3-as3.csv")




## THE MEDIAN OF THE 10 REPEATS
y1d1 <- aggregate(y1d1.dat[,"minb"], list(p=y1d1.dat$p), FUN=median) 
y1d2 <- aggregate(y1d2.dat[,"minb"], list(p=y1d2.dat$p), FUN=median) 
y1d3 <- aggregate(y1d3.dat[,"minb"], list(p=y1d3.dat$p), FUN=median) 

y2d1 <- aggregate(y2d1.dat[,"minb"], list(p=y2d1.dat$p), FUN=median) 
y2d2 <- aggregate(y2d2.dat[,"minb"], list(p=y2d2.dat$p), FUN=median) 
y2d3 <- aggregate(y2d3.dat[,"minb"], list(p=y2d3.dat$p), FUN=median) 

y3d1 <- aggregate(y3d1.dat[,"minb"], list(p=y3d1.dat$p), FUN=median) 
y3d2 <- aggregate(y3d2.dat[,"minb"], list(p=y3d2.dat$p), FUN=median) 
y3d3 <- aggregate(y3d3.dat[,"minb"], list(p=y3d3.dat$p), FUN=median) 

## 10 SIMULATION-BASED WORST-CASE BOUNDS (KEEP ONLY THE LOWER BOUNDS)
y3d1_wide <- spread(y3d1.dat[,c(1,5,4)], key = group, value = minb)
y3d2_wide <- spread(y3d2.dat[,c(1,5,4)], key = group, value = minb)
y3d3_wide <- spread(y3d3.dat[,c(1,5,4)], key = group, value = minb)

names(y3d1_wide)[-1] <- paste0("sim", 1:10)
names(y3d2_wide)[-1] <- paste0("sim", 1:10)
names(y3d3_wide)[-1] <- paste0("sim", 1:10)


## THE POINTS OF STUDIES
sei <- with(example2.1, TP/(TP+FN))
spi <- with(example2.1, TN/(TN+FP))
df <- data.frame(y = sei, x = 1-spi)

## THE SOP WITH BIAS
sed1 <- plogis(mu1+y1d1[,"x"])
spd1 <- plogis(mu2+y2d1[,"x"])
sed2 <- plogis(mu1+y1d2[,"x"])
spd2 <- plogis(mu2+y2d2[,"x"])
sed3 <- plogis(mu1+y1d3[,"x"])
spd3 <- plogis(mu2+y2d3[,"x"])


## THE SAUC
fit <- dtametasa.fc(example2.1, p = 1)
fit$sauc.ci
saucl <- fit$sauc.ci[1] - fit$sauc.ci[2]
saucu <- fit$sauc.ci[3] - fit$sauc.ci[1]

## THE SAUC OF THE MEDIAN OF RESULTS
saucd1 <- sapply(10:1, function(i) sauc.value(par= unlist(par), bias=y3d1[i,"x"]) )
sauc.ci1 <- sapply(1:10, function(i) saucci(par, fit$var.ml, saucd1[i]))
saucd1.lb <- sauc.ci1[1,]
saucd1.ub <- sauc.ci1[2,]
saucd2 <- sapply(10:1, function(i) sauc.value(par= unlist(par), bias=y3d2[i,"x"]) )
sauc.ci2 <- sapply(1:10, function(i) saucci(par, fit$var.ml, saucd2[i]))
saucd2.lb <- sauc.ci2[1,]
saucd2.ub <- sauc.ci2[2,]
saucd3 <- sapply(10:1, function(i) sauc.value(par= unlist(par), bias=y3d3[i,"x"]) )
sauc.ci3 <- sapply(1:10, function(i) saucci(par, fit$var.ml, saucd3[i]))
saucd3.lb <- sauc.ci3[1,]
saucd3.ub <- sauc.ci3[2,]

## THE SAUC OF 10 REPEAT
sauc.matd1 <- sauc.matd2 <- sauc.matd3 <- NULL
for(j in 2:11){
  mat <- as.matrix(sapply(10:1, function(i) sauc.value(par= unlist(par), bias=y3d1_wide[i,j]) ))
  sauc.matd1 <- cbind(sauc.matd1, mat)
}

for(j in 2:11){
  mat <- as.matrix(sapply(10:1, function(i) sauc.value(par= unlist(par), bias=y3d2_wide[i,j]) ))
  sauc.matd2 <- cbind(sauc.matd2, mat)
}

for(j in 2:11){
  mat <- as.matrix(sapply(10:1, function(i) sauc.value(par= unlist(par), bias=y3d3_wide[i,j]) ))
  sauc.matd3 <- cbind(sauc.matd3, mat)
}

dfall <- cbind.data.frame(p = rev(y3d1_wide$p), sauc.matd1, sauc.matd2, sauc.matd3)
colnames(dfall)[-1] <- c(paste0("simd1", 1:10), paste0("simd2", 1:10), paste0("simd3", 1:10))


p1 <- ggplot(dfall, aes(x=p)) + 
  geom_line(aes(y = simd11, color = "D4.1"), lty=1) +
  geom_line(aes(y = simd12, color = "D4.1"), lty=1) +
  geom_line(aes(y = simd13, color = "D4.1"), lty=1) +
  geom_line(aes(y = simd14, color = "D4.1"), lty=1) +
  geom_line(aes(y = simd15, color = "D4.1"), lty=1) +
  geom_line(aes(y = simd16, color = "D4.1"), lty=1) +
  geom_line(aes(y = simd17, color = "D4.1"), lty=1) +
  geom_line(aes(y = simd18, color = "D4.1"), lty=1) +
  geom_line(aes(y = simd19, color = "D4.1"), lty=1) +
  geom_line(aes(y = simd110, color = "D4.1"),lty=1)+
  geom_line(aes(y = simd21, color = "D4.2"), lty=5) +
  geom_line(aes(y = simd22, color = "D4.2"), lty=5) +
  geom_line(aes(y = simd23, color = "D4.2"), lty=5) +
  geom_line(aes(y = simd24, color = "D4.2"), lty=5) +
  geom_line(aes(y = simd25, color = "D4.2"), lty=5) +
  geom_line(aes(y = simd26, color = "D4.2"), lty=5) +
  geom_line(aes(y = simd27, color = "D4.2"), lty=5) +
  geom_line(aes(y = simd28, color = "D4.2"), lty=5) +
  geom_line(aes(y = simd29, color = "D4.2"), lty=5) +
  geom_line(aes(y = simd210,color = "D4.2"),lty=5) +
  geom_line(aes(y = simd31, color = "D4.3"), lty=3) +
  geom_line(aes(y = simd32, color = "D4.3"), lty=3) +
  geom_line(aes(y = simd33, color = "D4.3"), lty=3) +
  geom_line(aes(y = simd34, color = "D4.3"), lty=3) +
  geom_line(aes(y = simd35, color = "D4.3"), lty=3) +
  geom_line(aes(y = simd36, color = "D4.3"), lty=3) +
  geom_line(aes(y = simd37, color = "D4.3"), lty=3) +
  geom_line(aes(y = simd38, color = "D4.3"), lty=3) +
  geom_line(aes(y = simd39, color = "D4.3"), lty=3) +
  geom_line(aes(y = simd310,color = "D4.3"),lty=3) +
  scale_x_reverse(n.breaks = 10, name="Marginal selection probability") +
  scale_y_continuous(limits = c(0,1), name = "Worst-case upper bound") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.2, 0.2), legend.background = element_rect(fill = "white", color = "black"))+
  scale_colour_manual("Condition",
                      breaks = c("D4.1", "D4.2", "D4.3"),
                      values = c("#4daf4a", "#e41a1c", "#377eb8"),
                      guide_legend(override.aes = list(lty = c(1,5,3)))) +
  ggtitle("D. Lower bounds of the SAUC under 3 conditions")

## SROC PLOT A-C WHEN P =1, 0.8, 0.6, 0.4, 0.2
srocd1 <- plot10(y3d1_wide = y3d1_wide, title = "A. Lower bounds of the SROC curves under Condition (D4.1)", spd1 = spd1, sed1 = sed1)
srocd2 <- plot10(y3d1_wide = y3d2_wide, title = "B. Lower bounds of the SROC curves under Condition (D4.2)", spd1 = spd2, sed1 = sed2)
srocd3 <- plot10(y3d1_wide = y3d3_wide, title = "C. Lower bounds of the SROC curves under Condition (D4.3)", spd1 = spd3, sed1 = sed3)
p11 <- grid.arrange(srocd1, srocd2, srocd3, p1, ncol=2)
ggsave(filename = "eg2-1-app.eps", plot = p11, device = cairo_ps, width = 12, height = 12)



## PLOT 2
## PLOT D. MEDIAN OF THE SAUC
df.sauc <- data.frame(p = seq(1, 0.1, -0.1),
                      saucd1, saucd1.lb, saucd1.ub,
                      saucd2, saucd2.lb, saucd2.ub,
                      saucd3, saucd3.lb, saucd3.ub)

plot2D <- ggplot(df.sauc, aes(x = p)) + 
  geom_ribbon(mapping=aes(ymin=saucd1.lb,ymax=saucd1.ub, fill="D4.1"), colour="white", alpha=0.1)+
  geom_ribbon(mapping=aes(ymin=saucd2.lb,ymax=saucd2.ub, fill="D4.2"), colour="white", alpha=0.1)+
  geom_ribbon(mapping=aes(ymin=saucd3.lb,ymax=saucd3.ub, fill="D4.3"), colour="white", alpha=0.1)+
  geom_line( aes(y = saucd1, colour="D4.1"), lty=1, linewidth = 1) +
  geom_line( aes(y = saucd2, colour="D4.2"), lty=5, linewidth = 1) +
  geom_line( aes(y = saucd3, colour="D4.3"), lty=3, linewidth = 1) +
  scale_x_reverse(n.breaks = 10, name="Marginal selection probability") + 
  scale_y_continuous(limits = c(0,1), name = "Worst-case lower bound") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.2, 0.25), legend.background = element_rect(fill = "white", color = "black"))+
  labs(color='Conditions', fill='95% CI region', title = "")+
  scale_colour_manual("Condition",
                      breaks = c("D4.1","D4.2","D4.3"),
                      values = c("#4daf4a", "#e41a1c", "#377eb8"),
                      guide = guide_legend(override.aes = list(lty = c(1,5,3), linewidth = c(1,1,1))))+
  scale_fill_manual("95% CI region",
                    breaks = c("D4.1","D4.2","D4.3"),
                    values = c("#4daf4a", "#e41a1c", "#377eb8"),
                    guide = guide_legend(override.aes = list(color = c("white", "white", "white")))) +
  ggtitle("D. Lower bounds of the SAUC under 3 conditions")
## PLOT A-C MEDIAN OF THE SROC WHEN P =1, 0.8, 0.6, 0.4, 0.2
srocd1.m <- plot.ave(y3d1 = y3d1, title = "A. Lower bounds of the SROC curves under Condition (D4.1)", spd1 = spd1, sed1 = sed1)
srocd2.m <- plot.ave(y3d1 = y3d1, title = "B. Lower bounds of the SROC curves under Condition (D4.2)", spd1 = spd2, sed1 = sed2)
srocd3.m <- plot.ave(y3d1 = y3d1, title = "C. Lower bounds of the SROC curves under Condition (D4.3)", spd1 = spd3, sed1 = sed3)
p21 <- grid.arrange(srocd1.m, srocd2.m, srocd3.m, plot2D, ncol=2)
ggsave(filename = "eg2-2-app.eps", plot = p21, device = cairo_ps, width = 12, height = 12)


## TABLE2
tab1 <- cbind.data.frame(p =　seq(1, 0.1, -0.1),
                         A = sprintf("%.3f [%.3f, %.3f]", df.sauc[,2], df.sauc[,3], df.sauc[,4]),
                         B = sprintf("%.3f [%.3f, %.3f]", df.sauc[,5], df.sauc[,6], df.sauc[,7]),
                         C = sprintf("%.3f [%.3f, %.3f]", df.sauc[,8], df.sauc[,9], df.sauc[,10]))
colnames(tab1) <- c("p", "Condition (D4.1)", "Condition (D4.2)", "Condition (D4.3)")

sink("tab2-appp.tex")
kbl(tab1[,1:3],
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
