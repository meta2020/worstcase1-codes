##
## EXAMPLE 2 TROPNIN DATA
##

library(kableExtra)
library(mada)
# devtools::install_github("meta2020/dtametasa")
library(dtametasa)
library(gridExtra)
library(ggplot2)
library(tidyr)
source("rfunc/data.pre.R")
source("rfunc/plotfunc.R")
source("rfunc/sauc.R")
source("rfunc/sauc.ci.R")
source("rfunc/sroc.R")
source("rfunc/Piao2019.R")
source("rfunc/Li2021.R")

## GENERATE ANALYTICAL DATA FOR SAS PROGRAMMING---------------------------------
## DATA PREPROCESS
example2 <- read.csv("example-data/eg2-IVD.csv")


## CONTINUITY CORRECTION BY ADDING 0.5 TO THE STUDIES WITH 0 CELLS
example2.1 <- correction(data = example2, type = "single") 
# write.csv(example2.1, "anadata/eg2-IVD-cc.csv", row.names = F)

## DATA AFTER LOGIT-TRANSFORMATION
example2.2 <- logit.data(example2.1)
example2.2$fpr <- 1-example2.2$spec

## META-ANALYSIS USING THE REITSMA MODEL WITHOUT CONSIDERING PUBLICATION BIAS

ma2 <- reitsma(data=example2, method = "ml", correction.control = "single")

## THE ESTIMATES
par <- data.frame(
  mu1 = ma2$coefficients[1], mu2 = -ma2$coefficients[2],
  tau11 = ma2$Psi[1], tau22= ma2$Psi[4], tau12= -ma2$Psi[2])
# write.csv(par, "anadata/ml-par-IVD.csv", row.names = F)


##
## THE SAUC WITHOUT PUBLICATION BIAS--------------------------------------------

fit <- dtametasa.fc(example2.1, p = 1)
fit$sauc.ci
sprintf("SAUC (95.CI) = %.3f (%.3f to %.3f)", fit$sauc.ci[1], fit$sauc.ci[2], fit$sauc.ci[3])


## WORST-CASE BOUNDS UNDER CONDITION D4.1-D4.3, PRODUCED BY SAS
df <- read.csv("SAS/IVDresult2-K2000-R1.csv")

## PLOT 1: LOWER BOUND ONLY ----
##
df.as1 <- df[df$as==1,]
df.as2 <- df[df$as==2,]
df.as3 <- df[df$as==3,]
srocd1 <- plot.sorc.lb(df = example2.2, dfas = df.as1, title = "A. Lower bounds of the SROC curves under Condition (D4.1)")
srocd2 <- plot.sorc.lb(df = example2.2, dfas = df.as2, title = "B. Lower bounds of the SROC curves under Condition (D4.2)")
srocd3 <- plot.sorc.lb(df = example2.2, dfas = df.as3, title = "C. Lower bounds of the SROC curves under Condition (D4.3)")

## THE BIAS-ADJUSTED SAUC GIVEN MARGINAL SELECTION PROBABILITY
df$sauc.lb <- sapply(1:30, function(i) sauc.value.pb(par= unlist(par), bias=df[i,"minb"]) )
# sauc.ci <- sapply(1:30, function(i) saucci(par, fit$var.ml, df$sauc.lb [i]))
# df$sauc.clb <- sauc.ci[1,]

lgd1 <- c("Estimate withou PB","Condition (D4.1)","Condition (D4.2)","Condition (D4.3)")
dfl <- data.frame(
  p = c(df$p[1:10],df$p),
  est = c(rep(fit$sauc.ci[1],10),df$sauc.lb),
  grp = rep(lgd1, each=10)
)

saucd <- ggplot(dfl, aes(x = p, y = est, color = grp, linetype = grp)) +
  geom_line(size=1) +
  scale_y_continuous(limits = c(0,1),n.breaks = 10, name = "SAUC") +
  scale_x_reverse(n.breaks = 10, name="Marginal selection probability") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.25, 0.25), 
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank())+
  scale_colour_manual("",
                      breaks = lgd1,
                      values = c("grey50","#4daf4a", "#e41a1c", "#377eb8"))+
  scale_linetype_manual("",
                        breaks = lgd1,
                        values = c(2,1,2,3))+
  ggtitle("D. Lower bounds of the SAUC under 3 conditions")
p1 <- grid.arrange(srocd1, srocd2, srocd3, saucd, ncol=2)
ggsave(filename = "fig-tab/ivd/eg2-ivd.eps", plot = p1, device = cairo_ps, width = 12, height = 12)

## TABLE 1: LOWER BOUND OF THE SAUC ----

tab12 <- dfl %>% spread(key = grp, value = est)
tab1 <- cbind.data.frame(p =　seq(0.1, 1, 0.1),
                         M = round(nrow(example2)*(1/seq(0.1, 1, 0.1)-1)),
                         A = sprintf("%.3f", tab12$`Condition (D4.1)`),
                         B = sprintf("%.3f", tab12$`Condition (D4.2)`),
                         C = sprintf("%.3f", tab12$`Condition (D4.3)`))
colnames(tab1) <- c("$p$", "$M$","Lower bound", "Lower bound", "Lower bound")
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
    caption = "Example 2: the lower bound of the SAUC by the simulation-based bound.",
    label = "tab2-ivd",
    row.names = NA)%>%
  add_header_above(c("","", "Condition (D4.1)"=1, "Condition (D4.2)"=1, "Condition (D4.3)"=1), escape = FALSE) %>% 
  footnote(general = "
           $M$ indicates the number of unpublished studies.", 
           escape = FALSE, threeparttable = TRUE,  general_title = "")



## THE BIAS-ADJUSTED SAUC GIVEN MARGINAL SELECTION PROBABILITY
# saucd1 <- sapply(1:10, function(i) sauc.value.pb(par= unlist(par), bias=df.as1[i,"minb"]) )
# sauc.ci1 <- sapply(1:10, function(i) saucci(par, fit$var.ml, saucd1[i]))
# saucd1.lb <- sauc.ci1[1,]
# saucd1.ub <- sauc.ci1[2,]
# saucd2 <- sapply(1:10, function(i) sauc.value.pb(par= unlist(par), bias=df.as2[i,"minb"]) )
# sauc.ci2 <- sapply(1:10, function(i) saucci(par, fit$var.ml, saucd2[i]))
# saucd2.lb <- sauc.ci2[1,]
# saucd2.ub <- sauc.ci2[2,]
# saucd3 <- sapply(1:10, function(i) sauc.value.pb(par= unlist(par), bias=df.as3[i,"minb"]) )
# sauc.ci3 <- sapply(1:10, function(i) saucci(par, fit$var.ml, saucd3[i]))
# saucd3.lb <- sauc.ci3[1,]
# saucd3.ub <- sauc.ci3[2,]


## COMPARISION WITH OTHER METHODS ----------------------------------------------
load("rfunc/comp-ivd-dt.RData")

## PLOT 2a: c1=1, SENSITIVITY ----
## 
ldg2a <- c("Estimate withou PB",
           "Simulation-based bound under Condition (D4.1)",
           "Zhou et al.'s method when c1=1 and c2=0",
           "95% lower confidence interval of Zhou et al.'s method")
dfc1 <- data.frame(
  p = rep(df$p[1:10],4),
  est = c(dfl[dfl$grp %in% c("Estimate withou PB","Condition (D4.1)"),"est"],
          zhou.mod3[1,],zhou.mod3[2,]),
  grp = rep(ldg2a, each=10)
)
p2a <- ggplot(dfc1, aes(x = p, y = est, color = grp, linetype = grp)) +
  geom_line(size=1) +
  scale_y_continuous(limits = c(0,1),n.breaks = 10, name = "SAUC") +
  scale_x_reverse(n.breaks = 10, name="Marginal selection probability") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.35, 0.25), 
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank())+
  scale_colour_manual("",
                      breaks = ldg2a,
                      values = c("grey50","#e41a1c","#377eb8","#377eb8"))+
  scale_linetype_manual("",
                        breaks = ldg2a,
                        values = c(2,1,1,3))+
  ggtitle("A. Comparison scenario 1")

## PLOT 2b: c1=0, SPECIFICITY  ----
## 
ldg2b <- c("Estimate withou PB",
           "Simulation-based bound under Condition (D4.2)",
           "Zhou et al.'s method when c1=1 and c2=0",
           "95% lower confidence interval of Zhou et al.'s method")
dfc2 <- data.frame(
  p = rep(df$p[1:10],4),
  est = c(dfl[dfl$grp %in% c("Estimate withou PB","Condition (D4.2)"),"est"],
          zhou.mod4[1,],zhou.mod4[2,]),
  grp = rep(ldg2b, each=10)
)
p2b <- ggplot(dfc2, aes(x = p, y = est, color = grp, linetype = grp)) +
  geom_line(size=1) +
  scale_y_continuous(limits = c(0,1),n.breaks = 10, name = "SAUC") +
  scale_x_reverse(n.breaks = 10, name="Marginal selection probability") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.35, 0.25), 
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank())+
  scale_colour_manual("",
                      breaks = ldg2b,
                      values = c("grey50","#e41a1c","#377eb8","#377eb8"))+
  scale_linetype_manual("",
                        breaks = ldg2b,
                        values = c(2,1,1,3))+
  ggtitle("B. Comparison scenario 2")

## PLOT 2c: c1=c2, lnDOR ----
## 
ldg2c <- c("Estimate withou PB",
           "Simulation-based bound under Condition (D4.3)",
           "Zhou et al.'s method when c1=1 and c2=1",
           "95% lower confidence interval of Zhou et al.'s method")
dfc3 <- data.frame(
  p = rep(df$p[1:10],4),
  est = c(dfl[dfl$grp %in% c("Estimate withou PB","Condition (D4.3)"),"est"],
          zhou.mod2[1,],zhou.mod2[2,]),
  grp = rep(ldg2c, each=10)
)
p2c <- ggplot(dfc3, aes(x = p, y = est, color = grp, linetype = grp)) +
  geom_line(size=1) +
  scale_y_continuous(limits = c(0,1),n.breaks = 10, name = "SAUC") +
  scale_x_reverse(n.breaks = 10, name="Marginal selection probability") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.35, 0.25), 
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank())+
  scale_colour_manual("",
                      breaks = ldg2c,
                      values = c("grey50","#e41a1c","#377eb8","#377eb8"))+
  scale_linetype_manual("",
                        breaks = ldg2c,
                        values = c(2,1,1,3))+
  ggtitle("C. Comparison scenario 3")

## PLOT 2d: ESTIMATED ONE (SUPP) ----
##
lgd2d <- c("Estimate withou PB",
           "Simulation-based bound under 3 conditions",
           "Zhou et al.'s method when c1 c2 are estimated",
           "95% lower confidence interval of Zhou et al.'s method")
dfl1 <- dfl[dfl$grp %in% c("Estimate withou PB","Condition (D4.1)"),]
dfl2 <- dfl[dfl$grp %in% c("Condition (D4.2)"),]
dfl3 <- dfl[dfl$grp %in% c("Condition (D4.3)"),]
dfl4 <- cbind.data.frame(
  p = rep(dfl1$p[1:10],4),
  est1 = c(dfl1$est,zhou.mod1[1,],zhou.mod1[2,]),
  est2 = c(rep(NA,10),dfl2$est,rep(NA,20)),
  est3 = c(rep(NA,10),dfl3$est,rep(NA,20)),
  grp2 = rep(lgd2d,each=10))

p2d <- ggplot(dfl4, aes(x = p, y = est1, color = grp2, linetype = grp2)) +
  geom_line(size=1) +
  geom_line(aes(x = p, y = est2, color = grp2, linetype=grp2), size=1)+
  geom_line(aes(x = p, y = est3, color = grp2, linetype=grp2), size=1)+
  geom_point(aes(df.piao$p, df.piao$sauc, shape = "Piao et al.'s method"), color="black", size = 2.5)+
  geom_point(aes(df.li$p, df.li$sauc, shape = "Li et al.'s method"), color="black", size = 2.5)+
  scale_y_continuous(limits = c(0,1),n.breaks = 10, name = "SAUC") +
  scale_x_reverse(n.breaks = 10, name="Marginal selection probability") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.35, 0.25), 
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank())+
  scale_colour_manual("",
                      breaks = c(lgd2d, "Piao et al.","Li et al."),
                      values = c("grey50","#e41a1c","#377eb8","#377eb8"),
                      guide = guide_legend(override.aes = 
                                           list(pch = c(rep(NA,4)),
                                                lty = c(2,1,1,3))))+
  scale_linetype_manual("",
                        breaks = c(lgd2d, "Piao et al.'s method","Li et al.'s method"),
                        values = c(2,1,1,3,0,0)) +
  scale_shape_manual("",
                     breaks = c(lgd2d, "Piao et al.'s method","Li et al.'s method"),
                     values = c(rep(0,4),16,18))+
  ggtitle("D. Comparison scenario 4")

p24 <- grid.arrange(p2a, p2b, p2c, p2d, ncol=2)
ggsave(filename = "fig-tab/ivd/eg2-comp-ivd.eps", plot = p24, device = cairo_ps, width = 12, height = 12)


## COMPARISON OF SETTING DIFFERENT K ----
## 
b.MC.result1 <- read.csv("SAS/IVDresult2-K1000-R1.csv")
b.MC.result2 <- read.csv("SAS/IVDresult2-K2000-R1.csv")
b.MC.result3 <- read.csv("SAS/IVDresult2-K5000-R1.csv")
b.MC.result4 <- read.csv("SAS/IVDresult2-K20000-R1.csv")
dfk <- rbind.data.frame(b.MC.result1,b.MC.result2,b.MC.result3,b.MC.result4)
dfk.sauc <- sapply(1:nrow(dfk), function(i) sauc.value.pb(par= unlist(par), bias=dfk[i,"minb"]) )

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
    caption = "Example 2: the lower bounds of the SAUC by the simulation-based given different $K$.",
    label = "tab1",
    row.names = NA)  
## PLOT 3a: AS1 COMPARISON OF SETTING DIFFERENT K ----
## 
dfk2a <- dfk2[dfk2$as==1,]
p3a <- ggplot(dfk2a, aes(x = p, y = est, color = grp, linetype = grp)) +
  geom_line(size=1) +
  scale_y_continuous(limits = c(0.25,0.9),n.breaks = 13, name = "SAUC") +
  scale_x_reverse(n.breaks = 10, name="Marginal selection probability") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.3, 0.25), 
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank())+
  scale_colour_manual("",
                      breaks = lgd3,
                      values = c("#4daf4a","#e41a1c","#984ea3","#ff7f00"))+
  scale_linetype_manual("",
                        breaks = lgd3,
                        values = c(1,1,2,2))+
  ggtitle("A. Lower bounds of the SAUC under Condition (D4.1)")

## PLOT 3b: AS2 COMPARISON OF SETTING DIFFERENT K ----
## 
dfk2b <- dfk2[dfk2$as==2,]
p3b <- ggplot(dfk2b, aes(x = p, y = est, color = grp, linetype = grp)) +
  geom_line(size=1) +
  scale_y_continuous(limits = c(0.25,0.9),n.breaks = 13, name = "SAUC") +
  scale_x_reverse(n.breaks = 10, name="Marginal selection probability") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.3, 0.25), 
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank())+
  scale_colour_manual("",
                      breaks = lgd3,
                      values = c("#4daf4a","#e41a1c","#984ea3","#ff7f00"))+
  scale_linetype_manual("",
                        breaks = lgd3,
                        values = c(1,1,2,2))+
  ggtitle("B. Lower bounds of the SAUC under Condition (D4.2)")

## PLOT 3c: AS2 COMPARISON OF SETTING DIFFERENT K ----
## 
dfk2c <- dfk2[dfk2$as==3,]
p3c <- ggplot(dfk2c, aes(x = p, y = est, color = grp, linetype = grp)) +
  geom_line(size=1) +
  scale_y_continuous(limits = c(0.25,0.9),n.breaks = 13, name = "SAUC") +
  scale_x_reverse(n.breaks = 10, name="Marginal selection probability") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.3, 0.25), 
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank())+
  scale_colour_manual("",
                      breaks = lgd3,
                      values = c("#4daf4a","#e41a1c","#984ea3","#ff7f00"))+
  scale_linetype_manual("",
                        breaks = lgd3,
                        values = c(1,1,2,2))+
  ggtitle("C. Lower bounds of the SAUC under Condition (D4.3)")

## COMPARISON OF SETTING DIFFERENT K AND SAMPLES ----
##
b.MC.rep1 <- read.csv("SAS/IVDresult2-K1000-R10.csv")
b.MC.rep2 <- read.csv("SAS/IVDresult2-K2000-R10.csv")
b.MC.rep3 <- read.csv("SAS/IVDresult2-K5000-R10.csv")
b.MC.rep4 <- read.csv("SAS/IVDresult2-K20000-R10.csv")
dfr <- rbind.data.frame(b.MC.rep1,b.MC.rep2,b.MC.rep3,b.MC.rep4)
dfr.sauc <- sapply(1:nrow(dfr), function(i) sauc.value.pb(par= unlist(par), bias=dfr[i,"minb"]) )

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
p3d <- ggplot(dfr2a, aes(x=pp, y=est, fill=grp2)) + 
  geom_boxplot(alpha = 0.5)+
  scale_x_discrete(limits = rev(levels(dfr2a$pp)),name="Marginal selection probability") +
  scale_y_continuous(limits = c(0.25,0.9),n.breaks = 13, name = "SAUC") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.3, 0.25), legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank())+
  scale_fill_manual("",
                    breaks = lgd3,
                    values = c("#4daf4a","#e41a1c","#984ea3","#ff7f00")) +
  ggtitle("D. 10-time lower bounds of the SAUC under Condition (D4.1)")

## PLOT 3e: AS2 COMPARISON OF SETTING DIFFERENT K ----
## 
dfr2b <- dfr2[dfr2$as==2,]
p3e <- ggplot(dfr2b, aes(x=pp, y=est, fill=grp2)) + 
  geom_boxplot(alpha = 0.5)+
  scale_x_discrete(limits = rev(levels(dfr2b$pp)),name="Marginal selection probability") +
  scale_y_continuous(limits = c(0.25,0.9),n.breaks = 13, name = "SAUC") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.3, 0.25), legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank())+
  scale_fill_manual("",
                    breaks = lgd3,
                    values = c("#4daf4a","#e41a1c","#984ea3","#ff7f00")) +
  ggtitle("E. 10-time lower bounds of the SAUC under Condition (D4.2)")

## PLOT 3f: AS2 COMPARISON OF SETTING DIFFERENT K ----
## 
dfr2c <- dfr2[dfr2$as==3,]
p3f <- ggplot(dfr2c, aes(x=pp, y=est, fill=grp2)) + 
  geom_boxplot(alpha = 0.5)+
  scale_x_discrete(limits = rev(levels(dfr2c$pp)),name="Marginal selection probability") +
  scale_y_continuous(limits = c(0.25,0.9),n.breaks = 13, name = "SAUC") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"),
        legend.key = element_rect (fill = "white"),
        legend.position = c(0.3, 0.25), legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank())+
  scale_fill_manual("",
                    breaks = lgd3,
                    values = c("#4daf4a","#e41a1c","#984ea3","#ff7f00")) +
  ggtitle("E. 10-time lower bounds of the SAUC under Condition (D4.3)")
p3 <- grid.arrange(p3a, p3b, p3c, p3d, p3e, p3f, ncol=3)
ggsave(filename = "fig-tab/ivd/eg2-k-ivd.eps", plot = p3, device = cairo_ps, width = 18, height = 12)

