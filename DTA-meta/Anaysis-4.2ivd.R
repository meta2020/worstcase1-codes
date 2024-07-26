##
## EXAMPLE 2 IVD DATA
##

library(kableExtra)
library(mada)
# devtools::install_github("meta2020/dtametasa")
library(dtametasa)
library(gridExtra)
library(ggplot2)
library(tidyr)
library(latex2exp)
library(meta)
source("rfunc/data.pre.R")
source("rfunc/plotfunc.R")
source("rfunc/sauc.R")
source("rfunc/sauc.ci.R")
source("rfunc/sroc.R")
source("rfunc/Piao2019.R")
source("rfunc/Li2021.R")

## CREATE ANALYTICAL DATA FOR SAS PROGRAMMING ----
## DATA PREPROCESS
example2 <- read.csv("example-data/eg2-ivd.csv")


## CONTINUITY CORRECTION BY ADDING 0.5 TO THE STUDIES WITH 0 CELLS
example2.1 <- correction(data = example2, type = "single") 
# write.csv(example2.1, "anadata/eg2-ivd-cc.csv", row.names = F)

## DATA AFTER LOGIT-TRANSFORMATION
example2.2 <- logit.data(example2.1)
example2.2$fpr <- 1-example2.2$spec

## META-ANALYSIS USING THE REITSMA MODEL withoutT PUBLICATION BIAS ----

ma2 <- reitsma(data=example2, method = "ml", correction.control = "single")

## THE ESTIMATES
par <- data.frame(
  mu1 = ma2$coefficients[1], mu2 = -ma2$coefficients[2],
  tau11 = ma2$Psi[1], tau22= ma2$Psi[4], tau12= -ma2$Psi[2])
# write.csv(par, "anadata/ml-par-ivd.csv", row.names = F)


##
## THE SAUC withoutT PUBLICATION BIAS ----

fit <- dtametasa.fc(example2.1, p = 1)
sprintf("SAUC (95.CI) = %.3f (%.3f to %.3f)", fit$sauc.ci[1], fit$sauc.ci[2], fit$sauc.ci[3])

##
## PUBLICATION BIAS ANALYSIS ----

## PLOT 1: FUNNEL PLOT ----
## 

setEPS(width = 10, height = 4); postscript("fig-tab/ivd/eg2-fp-ivd.eps")
par(oma = c(4, 2, 0.2, 0.2), mfrow = c(1, 3), mar = c(4, 2, 2, 0.2))

## Sens
ml.se <- rma.uni(yi = y1, vi = v1, data = example2.2, measure="GEN", method="ML")
tf.se <- trimfill(ml.se, estimator="R0")
funnel(tf.se, xlab = "logit-sensitivity")
title("A", adj = 0, font.main = 1)
mtext("Standar Error", side = 2, line = 2, at = 1, cex = 0.7)

## Spec
ml.sp <- rma.uni(yi = y2, vi = v2, data = example2.2, measure="GEN", method="ML")
tf.sp <- trimfill(ml.sp, estimator="R0")
funnel(tf.sp, xlab = "logit-specificity")
title("B", adj = 0, font.main = 1)

## lnDOR
ml.lnDOR <- rma.uni(yi = with(example2.2, y1+y2), vi= with(example2.2, v1+v2), method="ML")
tf.lndor <- trimfill(ml.lnDOR, estimator="L0")
funnel(tf.lndor, xlab = "lnDOR")
title("C", adj = 0, font.main = 1)

par(mfrow = c(1,1))
dev.off()

## OTHER METHODS ----
load("rfunc/comp-ivd-dt.RData")

## WORST-CASE BOUNDS UNDER CONDITION D4.1-D4.3, PRODUCED BY SAS
df <- read.csv("SAS/ivd/IVDresult2-K2000-R1.csv")
## THE BIAS-ADJUSTED SAUC GIVEN MARGINAL SELECTION PROBABILITY
df$sauc.lb <- sapply(1:30, function(i) .sauc.value.pb(par= unlist(par), bias=df[i,"minb"]) )
lgd1 <- c("Estimate without PB","Condition (D4.1)","Condition (D4.2)","Condition (D4.3)")
dfl <- data.frame(
  p = c(df$p[1:10],df$p),
  est = c(rep(fit$sauc.ci[1],10),df$sauc.lb),
  grp = rep(lgd1, each=10)
)

## PLOT 2: ZHOU'S METHOD AND PROPOSED METHOD ----
## SROC 
srocc1 <- .plot.sorc.par(
  example2.2, zhou.par1,
  title = TeX("A. Zhou's method with $(a_1, a_2) = (1, 0)$"))
srocc2 <- .plot.sorc.par(
  example2.2, zhou.par2,
  title = TeX("B. Zhou's method with $(a_1, a_2) = (0, 1)$"))
srocc3 <- .plot.sorc.par(
  example2.2, zhou.par3,
  title = TeX("C. Zhou's method with $a_1=a_2$"))
srocd1 <- .plot.sorc.npar(
  example2.2, par, df[df$as==1,], 
  title = "D. Simulation-based lower bounds under Condition (D4.1)")
srocd2 <- .plot.sorc.npar(
  example2.2, par, df[df$as==2,], 
  title = "E. Simulation-based lower bounds under Condition (D4.2)")
srocd3 <- .plot.sorc.npar(
  example2.2, par, df[df$as==3,], 
  title = "F. Simulation-based lower bounds under Condition (D4.3)")

## SAUC
## a1=1, SENSITIVITY
ldg2a <- c("Estimate without PB",
           "Simulation-based bound under Condition (D4.1)",
           "Zhou et al.'s method when a1=1 and a2=0",
           "95% lower confidence interval of Zhou et al.'s method")
dfc1 <- data.frame(
  p = rep(df$p[1:10],4),
  est = c(dfl[dfl$grp %in% c("Estimate without PB","Condition (D4.1)"),"est"],
          zhou.sauc1[1,],zhou.sauc1[2,]),
  grp = rep(ldg2a, each=10)
)
p2a <- .plot.sauc.all(
  dfc = dfc1, ldg=ldg2a, 
  title = "G. Bounds for the SAUC given PB from sensitivity")

## a1=0, SPECIFICITY 
ldg2b <- c("Estimate without PB",
           "Simulation-based bound under Condition (D4.2)",
           "Zhou et al.'s method when a1=0 and a2=1",
           "95% lower confidence interval of Zhou et al.'s method")
dfc2 <- data.frame(
  p = rep(df$p[1:10],4),
  est = c(dfl[dfl$grp %in% c("Estimate without PB","Condition (D4.2)"),"est"],
          zhou.sauc2[1,],zhou.sauc2[2,]),
  grp = rep(ldg2b, each=10)
)
p2b <- .plot.sauc.all(
  dfc = dfc2, ldg=ldg2b, 
  title = "H. Bounds for the SAUC given PB from specificity")

## a1=a2, lnDOR 
ldg2c <- c("Estimate without PB",
           "Simulation-based bound under Condition (D4.3)",
           "Zhou et al.'s method when a1=a2",
           "95% lower confidence interval of Zhou et al.'s method")
dfc3 <- data.frame(
  p = rep(df$p[1:10],4),
  est = c(dfl[dfl$grp %in% c("Estimate without PB","Condition (D4.3)"),"est"],
          zhou.sauc3[1,],zhou.sauc3[2,]),
  grp = rep(ldg2c, each=10)
)
p2c <- .plot.sauc.all(
  dfc = dfc3, ldg=ldg2c, 
  title = "I. Bounds for the SAUC given PB from the lnDOR")

p1 <- grid.arrange(srocc1, srocc2, srocc3, srocd1, srocd2, srocd3, p2a, p2b, p2c, ncol=3)
ggsave(filename = "fig-tab/ivd/eg2-ivd2.eps", plot = p1, device = cairo_ps, width = 18, height = 18)



## PLOT 3 ----
srocd1a <- .plot.sorc.npar(
  example2.2, par, df[df$as==1,], 
  title = "A. Simulation-based lower bounds under Condition (D4.1)")
srocd2a <- .plot.sorc.npar(
  example2.2, par, df[df$as==2,], 
  title = "B. Simulation-based lower bounds under Condition (D4.2)")
srocd3a <- .plot.sorc.npar(
  example2.2, par, df[df$as==3,], 
  title = "C. Simulation-based lower bounds under Condition (D4.3)")
## SAUC
## a1=1, SENSITIVITY
ldg2a <- c("Estimate without PB",
           "Simulation-based bound under Condition (D4.1) ")
dfc1 <- data.frame(
  p = rep(df$p[1:10],2),
  est = c(dfl[dfl$grp %in% c("Estimate without PB","Condition (D4.1)"),"est"]),
  grp = rep(ldg2a, each=10)
)
p2aa <- .plot.sauc.single(
  dfc=dfc1, ldg=ldg2a,
  title="D. Simulation-based lower bounds under Condition (D4.1)")

## a1=0, SPECIFICITY 
ldg2b <- c("Estimate without PB",
           "Simulation-based bound under Condition (D4.2) ")
dfc2 <- data.frame(
  p = rep(df$p[1:10],2),
  est = c(dfl[dfl$grp %in% c("Estimate without PB","Condition (D4.2)"),"est"]),
  grp = rep(ldg2b, each=10)
)
p2ba <- .plot.sauc.single(
  dfc=dfc2, ldg=ldg2b,
  title="E. Simulation-based lower bounds under Condition (D4.2)")

## a1=a2, lnDOR 
ldg2c <- c("Estimate without PB",
           "Simulation-based bound under Condition (D4.3) ")
dfc3 <- data.frame(
  p = rep(df$p[1:10],4),
  est = c(dfl[dfl$grp %in% c("Estimate without PB","Condition (D4.3)"),"est"]),
  grp = rep(ldg2c, each=10)
)
p2ca <- .plot.sauc.single(
  dfc=dfc3, ldg=ldg2c,
  title="F. Simulation-based lower bounds under Condition (D4.3)")

p2 <- grid.arrange(srocd1a, srocd2a, srocd3a, p2aa, p2ba, p2ca, ncol=3)
ggsave(filename = "fig-tab/ivd/eg2-ivd3.eps", plot = p2, device = cairo_ps, width = 18, height = 12)

## TABLE 1: LOWER BOUND OF THE SAUC ----

tab12 <- dfl %>% spread(key = grp, value = est)
tab1 <- cbind.data.frame(p =　seq(0.1, 1, 0.1),
                         M = round(nrow(example2)*(1/seq(0.1, 1, 0.1)-1)),
                         Z1= sprintf("%.3f (%.3f)", rev(zhou.sauc1[1,]), rev(zhou.sauc1[3,])),
                         Z2= sprintf("%.3f (%.3f)", rev(zhou.sauc2[1,]), rev(zhou.sauc2[3,])),
                         Z3= sprintf("%.3f (%.3f)", rev(zhou.sauc3[1,]), rev(zhou.sauc3[3,])),
                         A = sprintf("%.3f", tab12$`Condition (D4.1)`),
                         B = sprintf("%.3f", tab12$`Condition (D4.2)`),
                         C = sprintf("%.3f", tab12$`Condition (D4.3)`))
colnames(tab1) <- c("$p$", "$M$","$(a_1^2,a_2^2)=(1,0)$", "$(a_1^2,a_2^2)=(0,1)$", "$(a_1^2,a_2^2)=(0.5,0.5)$",
                    "Condition (D4.1)", "Condition (D4.2)", "Condition (D4.3)")
tab11 <- tab1[order(tab1[,1], decreasing = T),]
rownames(tab11) <- NULL

kbl(tab11,
    format = "latex",
    longtable = F,
    booktabs = T,
    linesep = "",
    digits = 3,
    align = "r",
    escape = FALSE,
    caption = "Example 1: the lower bound of the SAUC by the method of Zhou et al. and simulation-based bound.",
    label = "tab2-trop",
    row.names = NA)%>%
  add_header_above(c("","","Method of Zhou et al."=3, "Lower simulation-based bound"=3), escape = FALSE) %>% 
  footnote(general = "
           $M$ indicates the number of unpublished studies.
           Data in the paretheses are the lower band of the 95% CI", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")


## COMPARISON OF SETTING DIFFERENT K ----
## 
b.MC.result1 <- read.csv("SAS/ivd/Ivdresult2-K1000-R1.csv")
b.MC.result2 <- read.csv("SAS/ivd/IVDresult2-K2000-R1.csv")
b.MC.result3 <- read.csv("SAS/ivd/Ivdresult2-K5000-R1.csv")
b.MC.result4 <- read.csv("SAS/ivd/Ivdresult2-K20000-R1.csv")
dfk <- rbind.data.frame(b.MC.result1,b.MC.result2,b.MC.result3,b.MC.result4)
dfk.sauc <- sapply(1:nrow(dfk), function(i) .sauc.value.pb(par= unlist(par), bias=dfk[i,"minb"]) )

lgd3 <- c("K=1000","K=2000","K=5000","K=20000")
dfk2 <- data.frame(
  p = dfk$p,
  pp = as.factor(dfk$p),
  est = dfk.sauc,
  as = dfk$as,
  grp = rep(lgd3,each = 30)
)
## TABLE 2: COMPARISON TABLE ----
## 
tab21 <- dfk2 %>% spread(key = grp, value = est)
tab22 <- tab21 %>% dplyr::arrange(as, desc(p))
tab22$as <- c("(D4.1)", rep("",9),"(D4.2)", rep("",9),"(D4.3)", rep("",9))
tab23 <- tab22[,c(3,1,4,5,7,6)]
colnames(tab23) <- c("Condition","$p$","$K=1000$","$K=2000$","$K=5000$","$K=20000$")

kbl(tab23, 
    format = "latex",
    longtable = F, 
    booktabs = T, 
    linesep = "",
    digits = 3,
    align = "r",
    escape = FALSE,
    caption = "Example 1: the lower bounds of the SAUC by the simulation-based given different $K$.",
    label = "tab1",
    row.names = NA)  

## PLOT 3a: AS1 COMPARISON OF SETTING DIFFERENT K ----
## 
dfk2a <- dfk2[dfk2$as==1,]
p3a <- .plot.kest(
  dfk=dfk2a, lgd=lgd3, ylim=seq(0.25,0.9,0.05),
  title="A. Lower bounds of the SAUC under Condition (D4.1)")

## PLOT 3b: AS2 COMPARISON OF SETTING DIFFERENT K ----
## 
dfk2b <- dfk2[dfk2$as==2,]
p3b <- .plot.kest(
  dfk=dfk2b, lgd=lgd3, ylim=seq(0.25,0.9,0.05),
  title="B. Lower bounds of the SAUC under Condition (D4.2)")

## PLOT 3c: AS2 COMPARISON OF SETTING DIFFERENT K ----
## 
dfk2c <- dfk2[dfk2$as==3,]
p3c <- .plot.kest(
  dfk=dfk2c, lgd=lgd3, ylim=seq(0.25,0.9,0.05),
  title="C. Lower bounds of the SAUC under Condition (D4.3)")

## COMPARISON OF SETTING DIFFERENT K AND SAMPLES ----
##
b.MC.rep1 <- read.csv("SAS/ivd/ivdresult2-K1000-R10.csv")
b.MC.rep2 <- read.csv("SAS/ivd/ivdresult2-K2000-R10.csv")
b.MC.rep3 <- read.csv("SAS/ivd/ivdresult2-K5000-R10.csv")
b.MC.rep4 <- read.csv("SAS/ivd/ivdresult2-K20000-R10.csv")
dfr <- rbind.data.frame(b.MC.rep1,b.MC.rep2,b.MC.rep3,b.MC.rep4)
dfr.sauc <- sapply(1:nrow(dfr), function(i) .sauc.value.pb(par= unlist(par), bias=dfr[i,"minb"]) )

dfr2 <- data.frame(
  p = dfr$p,
  pp = as.factor(dfr$p),
  est = dfr.sauc,
  as = dfr$as,
  grp = rep(lgd3,each = 300)
)
dfr2$grp2 <- factor(dfr2$grp, levels = lgd3, ordered = TRUE)

## PLOT 3d: AS1 COMPARISON OF SETTING DIFFERENT K ----
## 
dfr2a <- dfr2[dfr2$as==1,]
p3d <- .plot.kbox(
  dfk=dfr2a, lgd=lgd3, ylim=seq(0.25,0.9,0.05),
  title="D. 10-time lower bounds of the SAUC under Condition (D4.1)")

## PLOT 3e: AS2 COMPARISON OF SETTING DIFFERENT K ----
## 
dfr2b <- dfr2[dfr2$as==2,]
p3e <- .plot.kbox(
  dfk=dfr2b, lgd=lgd3, ylim=seq(0.25,0.9,0.05),
  title="E. 10-time lower bounds of the SAUC under Condition (D4.2)")

## PLOT 3f: AS2 COMPARISON OF SETTING DIFFERENT K ----
## 
dfr2c <- dfr2[dfr2$as==3,]
p3f <- .plot.kbox(
  dfk=dfr2c, lgd=lgd3, ylim=seq(0.25,0.9,0.05),
  title="F. 10-time lower bounds of the SAUC under Condition (D4.3)")

p3 <- grid.arrange(p3a, p3b, p3c, p3d, p3e, p3f, ncol=3)
ggsave(filename = "fig-tab/ivd/eg2-k-ivd.eps", plot = p3, device = cairo_ps, width = 18, height = 12)

