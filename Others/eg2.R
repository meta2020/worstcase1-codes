##
## Example 2 IVD data
##

library(kableExtra)
library(mada)
# devtools::install_github("meta2020/dtametasa")
library(dtametasa)
library(grid)
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
example2 <- read.csv("rawdata/data-schuetz.csv")
# example2 <- data("skin_tests")
# example2 <- skin_tests


## CONTINUITY CORRECTION BY ADDING 0.5 TO THE STUDIES WITH 0 CELLS
example2.1 <- correction(data = example2, type = "single") 
# write.csv(example2.1, "anadata/eg2-schuetz2-cc.csv", row.names = F)

## DATA AFTER LOGIT-TRANSFORMATION
example2.2 <- logit.data(example2.1)
(order(example2.2$v2+example2.2$v1)-order(example2.2$v2))
(order(example2.2$v2+example2.2$v1)-order(example2.2$v1))

example2.2$fpr <- 1-example2.2$spec
# plot(sort(1/example2.2$v1))
# text(x=1:nrow(example2.2), y = 0.1+sort(1/example2.2$v1), order(1/example2.2$v1))
# plot(sort(1/example2.2$v2))
# text(x=1:nrow(example2.2), y = 0.01+sort(1/example2.2$v2), order(1/example2.2$v2))
# boxplot(example2.2[, c("v1","v2")])
# table(example2.2$v2)
# 
par(oma = c(4, 2, 0.2, 0.2), mfrow = c(1, 3), mar = c(4, 2, 2, 0.2))

## Sens
ml.se <- rma.uni(yi = y1, vi = v1, data = example2.2, measure="GEN", method="ML")
tf.se <- trimfill(ml.se, estimator="R0")
funnel(tf.se, xlab = "logit-sensitivity")
title("A", adj = 0, font.main = 1)
mtext("Standar Error", side = 2, line = 2, at = 1, cex = 0.7)

## Spec
ml.sp <- rma.uni(yi = y2, vi = v2, data = example2.2, measure="GEN", method="ML")
tf.sp <- trimfill(ml.sp, estimator="L0")
funnel(tf.sp, xlab = "logit-specificity")
title("B", adj = 0, font.main = 1)

## lnDOR
ml.lnDOR <- rma.uni(yi = with(example2.2, y1+y2), vi= with(example2.2, v1+v2), method="ML")
tf.lndor <- trimfill(ml.lnDOR, estimator="L0")
funnel(tf.lndor, xlab = "lnDOR")
title("C", adj = 0, font.main = 1)

par(mfrow = c(1,1))

## META-ANALYSIS USING THE REITSMA MODEL WITHOUT CONSIDERING PUBLICATION BIAS

ma2 <- reitsma(data=example2, method = "ml", correction.control = "single")

## THE ESTIMATES
par <- data.frame(
  mu1 = ma2$coefficients[1], mu2 = -ma2$coefficients[2],
  tau11 = ma2$Psi[1], tau22= ma2$Psi[4], tau12= -ma2$Psi[2])
write.csv(par, "anadata/ml-par-schuetz2.csv", row.names = F)
plogis(ma2$coefficients)

## THE SAUC WITHOUT PUBLICATION BIAS

fit <- dtametasa.fc(example2.1, p = 1)
fit$sauc.ci
sprintf("SAUC (95.CI) = %.3f (%.3f to %.3f)", 
        fit$sauc.ci[1], fit$sauc.ci[2], fit$sauc.ci[3])


## WORST-CASE BOUNDS UNDER CONDITION D4.1-D4.3, PRODUCED BY SAS
df <- read.csv("SASresult/result2-1000schuetz2.csv")
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


## COMPARISION WITH OTHER METHODS

## METHOD OF ZHOU ET AL
p.10 <- seq(1, 0.1, -0.1)
zhou.mod1 <- sapply(p.10, function(p) {
  
  opt2 <- try(dtametasa.rc(data=example2, p=p, c1.square0 = 0.5, beta.interval = c(-2,2)))
  if(inherits(opt2, "try-error")) rep(NA,3) else c(opt2$sauc.ci)
  
})


zhou.mod2 <- sapply(p.10, function(p) {
  
  opt2 <- try(dtametasa.fc(data=example2, p=p, c1.square = 0.5, beta.interval = c(-2,2)))
  if(inherits(opt2, "try-error")) rep(NA,3) else c(opt2$sauc.ci)  
  
})


zhou.mod3 <- sapply(p.10, function(p) {
  
  opt2 <- try(dtametasa.fc(data=example2, p=p, c1.square = 1, beta.interval = c(-2,2)))
  if(inherits(opt2, "try-error")) rep(NA,3) else c(opt2$sauc.ci)  
  
})

zhou.mod4 <- sapply(p.10, function(p) {
  
  opt2 <- try(dtametasa.fc(data=example2, p=p, c1.square = 0, beta.interval = c(-2,2)))
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
  geom_ribbon(mapping=aes(ymin=saucd2,ymax=saucd2.lb, color="Proposal under Condition (D4.2)"), fill="#e41a1c", alpha=0.05, lty=3)+
  geom_ribbon(mapping=aes(ymin=saucd3,ymax=saucd3.lb, color="Proposal under Condition (D4.3)"), fill="#377eb8", alpha=0.05, lty=3)+
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
        legend.position = c(0.3, 0.3), legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank(),
        )+
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

df.var <- data.frame(variance = c(example2.2$v1, example2.2$v2), group = c(rep("v1", nrow(example2.2)),rep("v2", nrow(example2.2))))
pvar <- ggplot(df.var, aes(x=group, y=variance, fill=group)) + 
  geom_boxplot()

df.t <- data.frame(variance = c(example2.2$y1/sqrt(example2.2$v1), 
                                example2.2$y1/sqrt(example2.2$v1),
                                (example2.2$y1+example2.2$y2)/sqrt(example2.2$v1+example2.2$v2)), 
                   group = c(rep("t-sen", nrow(example2.2)),
                             rep("t-spe", nrow(example2.2)),
                             rep("t-lndor", nrow(example2.2))))
pt <- ggplot(df.t, aes(x=group, y=variance, fill=group)) + 
  geom_boxplot()

pall <- grid.arrange(pvar, pt, p3,ncol=2, top = textGrob("SKin"))
ggsave(filename = "var-skin.eps", plot = pall, device = cairo_ps, width = 12, height = 12)

