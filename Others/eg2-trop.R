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
source("Piao2019.R")
source("Li2021.R")


## DATA PREPROCESS
example2 <- read.csv("rawdata/eg2-trop.csv")


## CONTINUITY CORRECTION BY ADDING 0.5 TO THE STUDIES WITH 0 CELLS
example2.1 <- correction(data = example2, type = "single") 
write.csv(example2.1, "anadata/eg2-trop-cc.csv", row.names = F)

## DATA AFTER LOGIT-TRANSFORMATION
example2.2 <- logit.data(example2.1)
example2.2$fpr <- 1-example2.2$spec

## META-ANALYSIS USING THE REITSMA MODEL WITHOUT CONSIDERING PUBLICATION BIAS

ma2 <- reitsma(data=example2, method = "ml", correction.control = "single")

## THE ESTIMATES
par <- data.frame(
  mu1 = ma2$coefficients[1], mu2 = -ma2$coefficients[2],
  tau11 = ma2$Psi[1], tau22= ma2$Psi[4], tau12= -ma2$Psi[2])
write.csv(par, "anadata/ml-par-trop.csv", row.names = F)

## THE SAUC WITHOUT PUBLICATION BIAS

fit <- dtametasa.fc(example2.1, p = 1)
fit$sauc.ci
sprintf("SAUC (95.CI) = %.3f (%.3f to %.3f)", 
        fit$sauc.ci[1], fit$sauc.ci[2], fit$sauc.ci[3])


## WORST-CASE BOUNDS UNDER CONDITION D4.1-D4.3, PRODUCED BY SAS
df <- read.csv("SASresult/result2-2000trop.csv")
df.as1 <- df[df$as==1,]
df.as2 <- df[df$as==2,]
df.as3 <- df[df$as==3,]


## THE BIAS-ADJUSTED SAUC GIVEN MARGINAL SELECTION PROBABILITY
saucd1 <- sapply(1:10, function(i) sauc.value.pb(par= unlist(par), bias=df.as1[i,"minb"]) )
sauc.ci1 <- sapply(1:10, function(i) saucci(par, fit$var.ml, saucd1[i]))
saucd1.lb <- sauc.ci1[1,]
saucd1.ub <- sauc.ci1[2,]
saucd2 <- sapply(1:10, function(i) sauc.value.pb(par= unlist(par), bias=df.as2[i,"minb"]) )
sauc.ci2 <- sapply(1:10, function(i) saucci(par, fit$var.ml, saucd2[i]))
saucd2.lb <- sauc.ci2[1,]
saucd2.ub <- sauc.ci2[2,]
saucd3 <- sapply(1:10, function(i) sauc.value.pb(par= unlist(par), bias=df.as3[i,"minb"]) )
sauc.ci3 <- sapply(1:10, function(i) saucci(par, fit$var.ml, saucd3[i]))
saucd3.lb <- sauc.ci3[1,]
saucd3.ub <- sauc.ci3[2,]



## PLOT 1
df.sauc <- data.frame(p = seq(1, 0.1, -0.1),
                      sauc = rep(fit$sauc.ci[1], 10),
                      sauc.lb = rep(fit$sauc.ci[2], 10),
                      sauc.ub = rep(fit$sauc.ci[3], 10),
                      saucd1, saucd1.lb,
                      saucd2, saucd2.lb,
                      saucd3, saucd3.lb)

plot.D <- ggplot(df.sauc, aes(x = p)) + 
  geom_ribbon(mapping=aes(ymin=sauc.lb,ymax=sauc.ub, color="Estimation without PB"), fill="grey50", alpha=0.1, lty=3)+
  geom_ribbon(mapping=aes(ymin=saucd1,ymax=saucd1.lb, color="Simulation-based bounds under Condition (D4.1)"), fill="#4daf4a", alpha=0.1, lty=3)+
  geom_ribbon(mapping=aes(ymin=saucd2,ymax=saucd2.lb, color="Simulation-based bounds under Condition (D4.2)"), fill="#e41a1c", alpha=0.1, lty=3)+
  geom_ribbon(mapping=aes(ymin=saucd3,ymax=saucd3.lb, color="Simulation-based bounds under Condition (D4.3)"), fill="#377eb8", alpha=0.1, lty=3)+
  geom_line( aes(y = sauc, colour="Estimation without PB"), lty=2, linewidth = 1) +
  geom_line( aes(y = saucd1, colour="Simulation-based bounds under Condition (D4.1)"), lty=1, linewidth = 1) +
  geom_line( aes(y = saucd2, colour="Simulation-based bounds under Condition (D4.2)"), lty=2, linewidth = 1) +
  geom_line( aes(y = saucd3, colour="Simulation-based bounds under Condition (D4.3)"), lty=3, linewidth = 1) +
  scale_x_reverse(n.breaks = 10, name="Marginal selection probability") + 
  scale_y_continuous(limits = c(0,1), name = "Worst-case lower bound") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.4, 0.25), legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank())+
  scale_colour_manual("",
                      breaks = c("Simulation-based bounds under Condition (D4.1)",
                                 "Simulation-based bounds under Condition (D4.2)",
                                 "Simulation-based bounds under Condition (D4.3)", 
                                 "Estimation without PB"),
                      values = c("#4daf4a", "#e41a1c", "#377eb8", "grey50"),
                      guide = guide_legend(override.aes = 
                                             list(lty = c(1,2,3,2), 
                                                  fill = c("#4daf4a", "#e41a1c", "#377eb8", "grey50"))))+
  ggtitle("D. Lower bounds of the SAUC under 3 conditions")
## PLOT A-C MEDIAN OF THE SROC WHEN P =1, 0.8, 0.6, 0.4, 0.2
srocd1 <- plot.sorc.lb(df = example2.2, dfas = df.as1, title = "A. Lower bounds of the SROC curves under Condition (D4.1)")
srocd2 <- plot.sorc.lb(df = example2.2, dfas = df.as2, title = "B. Lower bounds of the SROC curves under Condition (D4.2)")
srocd3 <- plot.sorc.lb(df = example2.2, dfas = df.as3, title = "C. Lower bounds of the SROC curves under Condition (D4.3)")
p21 <- grid.arrange(srocd1, srocd2, srocd3, plot.D, ncol=2)
ggsave(filename = "eg2-1.eps", plot = p21, device = cairo_ps, width = 12, height = 12)


## TABLE2
tab1 <- cbind.data.frame(p =　seq(1, 0.1, -0.1),
                         A = sprintf("%.3f (%.3f)", df.sauc$saucd1, df.sauc$saucd1.lb),
                         B = sprintf("%.3f (%.3f)", df.sauc$saucd2, df.sauc$saucd2.lb),
                         C = sprintf("%.3f (%.3f)", df.sauc$saucd3, df.sauc$saucd3.lb))
colnames(tab1) <- c("p", "LB (95\\% CLB)", "LB (95\\% CLB)", "LB (95\\% CLB)")

sink("tab2.tex")
kbl(tab1,
    format = "latex",
    longtable = F,
    booktabs = T,
    linesep = "",
    digits = 3,
    align = "r",
    escape = FALSE,
    caption = "Example 2: the lower bound of the SAUC by the simulation-based bound.",
    label = "tab2",
    row.names = NA)%>%
  add_header_above(c("", "Condition (D4.1)"=1, "Condition (D4.2)"=1, "Condition (D4.3)"=1), escape = FALSE) %>% 
  footnote(general = "
           LB indicates lower bound; CLB indicates the confidence lower band of the LB.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")
sink()


## COMPARISION WITH OTHER METHODS

## METHOD OF ZHOU ET AL
p.10 <- seq(1, 0.1, -0.1)
zhou.mod1 <- sapply(p.10, function(p) {
  
  opt2 <- try(dtametasa.rc(data=example2, p=p, c1.square0 = 0.5, beta.interval = c(0,2)))
  if(inherits(opt2, "try-error")) rep(NA,3) else c(opt2$sauc.ci)
  
})


zhou.mod2 <- sapply(p.10, function(p) {
  
  opt2 <- try(dtametasa.fc(data=example2, p=p, c1.square = 0.5, beta.interval = c(0,2)))
  if(inherits(opt2, "try-error")) rep(NA,3) else c(opt2$sauc.ci)  
  
})


zhou.mod3 <- sapply(p.10, function(p) {
  
  opt2 <- try(dtametasa.fc(data=example2, p=p, c1.square = 1, beta.interval = c(0,2)))
  if(inherits(opt2, "try-error")) rep(NA,3) else c(opt2$sauc.ci)  
  
})

zhou.mod4 <- sapply(p.10, function(p) {
  
  opt2 <- try(dtametasa.fc(data=example2, p=p, c1.square = 0, beta.interval = c(0,2)))
  if(inherits(opt2, "try-error")) rep(NA,3) else c(opt2$sauc.ci)  
  
})


## METHOD OF PIAO ET AL
Piao.mod <- Piao2019(data = example2.2, init = c(unlist(par), -0.1, -0.1, rep(-0.1, 3)) )
Piao.sauc <- SAUC(Piao.mod$par[1:5])
Piao.p <- nrow(example2)/(nrow(example2)+ round(max(Piao.mod$m)))
df.piao <- data.frame(p = Piao.p, y = Piao.sauc)

## MEHTHOD OF LI ET AL
Li.mod <- Li2021(data = example2.2, eta1.tilde = Piao.mod$par[8:10], eta0 = Piao.mod$par)
Li.sauc <- SAUC(Li.mod$par[1:5])
Li.p <- Li.mod$alpha
df.li <- data.frame(p = Li.p, y = Li.sauc)

## PLOT
df.sauc2 <- cbind.data.frame(
  p = seq(1, 0.1, -0.1), sauc = rep(fit$sauc.ci[1], 10),
  zhou.mod1 = t(zhou.mod1),
  zhou.mod2 = t(zhou.mod2),
  zhou.mod3 = t(zhou.mod3),
  zhou.mod4 = t(zhou.mod4),
  saucd1=saucd1, saucd2=saucd2, saucd3=saucd3)



p3 <- ggplot(df.sauc2, aes(x = p)) + 
  geom_line( aes(y = sauc, colour="Estimate without PB"), lty=1, linewidth = 1) +
  geom_ribbon(mapping=aes(ymin=zhou.mod1.sauc.lb,ymax=zhou.mod1.sauc.ub, color="Zhou et al. (2023) with estimated c1, c2"), fill="lightpink", alpha=0.05, lty=3)+
  geom_ribbon(mapping=aes(ymin=zhou.mod2.sauc.lb,ymax=zhou.mod2.sauc.ub, color="Zhou et al. (2023) with fixed c1 = c2"), fill="maroon", alpha=0.05, lty=3)+
  geom_ribbon(mapping=aes(ymin=zhou.mod3.sauc.lb,ymax=zhou.mod3.sauc.ub, color="Zhou et al. (2023) with fixed c1 = 1, c2 = 0"), fill="lightsalmon", alpha=0.05, lty=3)+
  geom_ribbon(mapping=aes(ymin=zhou.mod4.sauc.lb,ymax=zhou.mod4.sauc.ub, color="Zhou et al. (2023) with fixed c1 = 0, c2 = 1"), fill="mediumpurple", alpha=0.05, lty=3)+
  geom_ribbon(mapping=aes(ymin=saucd1,ymax=saucd1.lb, color="Proposal under Condition (D4.1)"), fill="#4daf4a", alpha=0.05, lty=3)+
  geom_ribbon(mapping=aes(ymin=saucd2,ymax=saucd2.lb, color="Proposal under Condition (D4.1)"), fill="#e41a1c", alpha=0.05, lty=3)+
  geom_ribbon(mapping=aes(ymin=saucd3,ymax=saucd3.lb, color="Proposal under Condition (D4.1)"), fill="#377eb8", alpha=0.05, lty=3)+
  geom_line( aes(y = zhou.mod1.sauc, colour="Zhou et al. (2023) with estimated c1, c2"), lty=1, linewidth = 1) +
  geom_line( aes(y = zhou.mod2.sauc, colour="Zhou et al. (2023) with fixed c1 = c2"), lty=1, linewidth = 1) +
  geom_line( aes(y = zhou.mod3.sauc, colour="Zhou et al. (2023) with fixed c1 = 1, c2 = 0"), lty=2, linewidth = 1) +
  geom_line( aes(y = zhou.mod4.sauc, colour="Zhou et al. (2023) with fixed c1 = 0, c2 = 1"), lty=2, linewidth = 1) +
  geom_line( aes(y = saucd1, colour="Proposal under Condition (D4.1)"), lty=1, linewidth = 1) +
  geom_line( aes(y = saucd2, colour="Proposal under Condition (D4.2)"), lty=2, linewidth = 1) +
  geom_line( aes(y = saucd3, colour="Proposal under Condition (D4.3)"), lty=3, linewidth = 1) +
  geom_point(data = df.piao, aes(p, y, color = "Piao et al. (2019)"), shape = 16, size = 1.5)+
  geom_point(data = df.li, aes(p, y, color = "Li et al. (2021)"), shape = 18, size = 2.5)+
  scale_x_reverse(n.breaks = 10, name="Marginal selection probability") + 
  scale_y_continuous(limits = c(0,1), name = "SAUC") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.25, 0.25), legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank())+
  scale_colour_manual("",
                      breaks = c("Piao et al. (2019)",
                                 "Li et al. (2021)",
                                 "Zhou et al. (2023) with estimated c1, c2",
                                 "Zhou et al. (2023) with fixed c1 = c2",
                                 "Zhou et al. (2023) with fixed c1 = 1, c2 = 0",
                                 "Zhou et al. (2023) with fixed c1 = 0, c2 = 1",
                                 "Proposal under Condition (D4.1)",
                                 "Proposal under Condition (D4.2)",
                                 "Proposal under Condition (D4.3)", 
                                 "Estimate without PB"),
                      values = c("black","black", "lightpink", "maroon", "lightsalmon", "mediumpurple", "#4daf4a", "#e41a1c", "#377eb8", "grey50"),
                      guide = guide_legend(override.aes = 
                                             list(pch = c(16,18, rep(NA,8)),
                                                  lty = c(0,0, 1,1,2,2,1,2,3,1), 
                                                  fill = c("white","white", "lightpink", "maroon", "lightsalmon", "mediumpurple", "#4daf4a", "#e41a1c", "#377eb8", "grey50"))))


ggsave(filename = "comp-trop.eps", plot = p3, device = cairo_ps, width = 8, height = 8)


## PLOT 3: COMPARISON OF SETTING DIFFERENT K
b.MC.result1 <- read.csv("SASresult/result2-1000trop.csv")
b.MC.result2 <- read.csv("SASresult/result2-2000trop.csv")
b.MC.result3 <- read.csv("SASresult/result2-20000trop.csv")

p2.1 <- plot.sauc.comp(as.num = 1, title = "A. Lower bounds of the SAUC under Condition (D4.1)")
p2.2 <- plot.sauc.comp(as.num = 2, title = "B. Lower bounds of the SAUC under Condition (D4.2)")
p2.3 <- plot.sauc.comp(as.num = 3, title = "C. Lower bounds of the SAUC under Condition (D4.3)")
p2k <- grid.arrange(p2.1, p2.2, p2.3, ncol=2)
ggsave(filename = "eg2-k.eps", plot = p2k, device = cairo_ps, width = 12, height = 12)