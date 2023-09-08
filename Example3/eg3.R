##
## LOAD DATA
##
library(readxl)
library(kableExtra)
files.sources <- list.files(path = "R/")
x <- sapply(paste0("R/", files.sources), source)

## MEDIAN FU 
med.data <- read_excel("Ki67.xlsx", sheet = "MCT")
med.data$mty <- med.data$mct_mo/12

## OBTAIN ETA
eta <- cens.eta(data = med.data, med.year = mty,  n1 = n1, n0 = n0, s1_med = s1_mct, s0_med = s0_mct)$par


## OS DATA

os.data  <- read_excel("Ki67.xlsx", sheet = "OS")
os.data$ty <- os.data$t/12

## HR DATA

dataHR <-read_excel("Ki67.xlsx", sheet = "HR")
dataHR$u_lnHR <- log(dataHR$HR)
ln_ci.low   <- log(dataHR$ci.low)
ln_ci.up    <- log(dataHR$ci.up)
se          <- (ln_ci.up - ln_ci.low)/(2 * qnorm(0.975))
dataHR$v_lnHR <- se^2

## MERGED DATA (OS+HR)
os_hr <- merge(os.data, dataHR, all.x = TRUE)

## GENERATE 3RD YEAR DATA
data3 <- convert.dt(
  data = os_hr, tK = 3, study = study, ty = ty, 
  n1 = n1, n0 = n0, s1 = s1, s0 = s0, u_lnHR = u_lnHR, v_lnHR = v_lnHR, 
  eta = eta)

## REMOVE OBSERVATIONS WITH MISSING VALUES
data3 <- na.omit(data3)

## GENERATE 5RD YEAR DATA
data5 <- convert.dt(
  data = os_hr, tK = 5, study = study, ty = ty, 
  n1 = n1, n0 = n0, s1 = s1, s0 = s0, u_lnHR = u_lnHR, v_lnHR = v_lnHR, 
  eta = eta)

## REMOVE OBSERVATIONS WITH MISSING VALUES
data5 <- na.omit(data5)

data <- rbind.data.frame(data3,data5)

write.csv(data, "data-35y.csv", row.names = F)

## META-ANALYSIS USING THE HZ MODEL WITHOUT CONSIDERING PUBLICATION BIAS

## YEAR 3
data <- data3
y1  <- data$u_sen
y2  <- data$u_spe 
v1  <- data$v_sen
v2  <- data$v_spe 
v12 <- data$v_senspe 

ma3.llk  <- function(par) llk.BNM.ml(par, y1, y2, v1, v2, v12)
ma3 <- nlminb(c(rep(0.1,4), -0.1), ma3.llk)
par3 <- ma3$par
num.hessian <- hessian(ma3.llk, ma3$par)
rownames(num.hessian) <- colnames(num.hessian) <- c("u1", "u2", "t1", "t2", "r1")
var.ml3 <- solve(num.hessian)
sauc3 <- SAUC.ci(par3, var.ml3, sauc.type = c("sroc"), ci.level = 0.95)
l31 <- sauc3[2]-sauc3[1]
l32 <- sauc3[3]-sauc3[1]

## THE ESTIMATES
par3 <- data.frame(
  t = 3,
  theta1 = par3[1], theta2 = par3[2],
  tau11 = par3[3]^2, tau22= par3[4]^2, tau12= prod(par3[3:5]))
write.csv(par3, "ml-par3.csv", row.names = F)

## YEAR 5
data <- data5
y1  <- data$u_sen
y2  <- data$u_spe 
v1  <- data$v_sen
v2  <- data$v_spe 
v12 <- data$v_senspe 

ma3.llk  <- function(par) llk.BNM.ml(par, y1, y2, v1, v2, v12)
ma3 <- nlminb(c(rep(0.1,4), -0.1), ma3.llk)
par3 <- ma3$par
num.hessian <- hessian(ma3.llk, ma3$par)
rownames(num.hessian) <- colnames(num.hessian) <- c("u1", "u2", "t1", "t2", "r1")
var.ml3 <- solve(num.hessian)
sauc5 <- SAUC.ci(par3, var.ml3, sauc.type = c("sroc"), ci.level = 0.95)
l51 <- sauc5[2]-sauc5[1]
l52 <- sauc5[3]-sauc5[1]

## THE ESTIMATES
par3 <- data.frame(
  t = 5,
  theta1 = par3[1], theta2 = par3[2],
  tau11 = par3[3]^2, tau22= par3[4]^2, tau12= prod(par3[3:5]))
write.csv(par3, "ml-par5.csv", row.names = F)

## LOAD ESTIMATED BOUNDS
bound1 <- read.csv("MCbound3-E1.csv")
bound2 <- read.csv("MCbound3-E2.csv")
bound3 <- read.csv("MCbound3-E3.csv")
bound4 <- read.csv("MCbound5-E1.csv")
bound5 <- read.csv("MCbound5-E2.csv")
bound6 <- read.csv("MCbound5-E3.csv")

p <- seq(1,0.1,-0.1)

## Under CONDITION D1-D3
## PLOT 2

df1 <- data.frame(p=p, bound = rev(bound1$SAUC_LB), lb = rev(bound1$SAUC_LB+l31), ub = rev(bound1$SAUC_LB+l32))
df2 <- data.frame(p=p, bound = rev(bound2$SAUC_LB), lb = rev(bound2$SAUC_LB+l31), ub = rev(bound2$SAUC_LB+l32))
df3 <- data.frame(p=p, bound = rev(bound3$SAUC_LB), lb = rev(bound3$SAUC_LB+l31), ub = rev(bound3$SAUC_LB+l32))
df4 <- data.frame(p=p, bound = rev(bound4$SAUC_LB), lb = rev(bound4$SAUC_LB+l51), ub = rev(bound4$SAUC_LB+l52))
df5 <- data.frame(p=p, bound = rev(bound5$SAUC_LB), lb = rev(bound5$SAUC_LB+l51), ub = rev(bound5$SAUC_LB+l52))
df6 <- data.frame(p=p, bound = rev(bound6$SAUC_LB), lb = rev(bound6$SAUC_LB+l51), ub = rev(bound6$SAUC_LB+l52))

plot1 <- ggplot(df1, aes(x = p)) + 
  geom_ribbon(mapping=aes(ymin=lb,ymax=ub), fill="black", colour="white", alpha=0.1)+
  geom_line(aes(y = bound), colour="black") +
  geom_point(aes(y = bound), colour ="black", pch=19, size=2) + 
  scale_x_reverse(n.breaks = 10, name="Overall selection probability") + 
  scale_y_continuous(limits = c(0.4,0.7), n.breaks = 5, name = "Worst-case lower bound") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"))+
  ggtitle("(A) Year 3: Condition (E1)") 

plot2 <- ggplot(df2, aes(x = p)) + 
  geom_ribbon(mapping=aes(ymin=lb,ymax=ub), fill="black", colour="white", alpha=0.1)+
  geom_line(aes(y = bound), colour="black") +
  geom_point(aes(y = bound), colour ="black", pch=19, size=2) + 
  scale_x_reverse(n.breaks = 10, name="Overall selection probability") + 
  scale_y_continuous(limits = c(0.4,0.7), n.breaks = 5, name = "Worst-case lower bound") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"))+
  ggtitle("(B) Year 3: Condition (E2)") 

# plot3 <- ggplot(df3, aes(x = p)) + 
#   geom_ribbon(mapping=aes(ymin=lb,ymax=ub), fill="black", colour="white", alpha=0.1)+
#   geom_line(aes(y = bound), colour="black") +
#   geom_point(aes(y = bound), colour ="black", pch=19, size=2) + 
#   scale_x_reverse(n.breaks = 10, name="Overall selection probability") + 
#   scale_y_continuous(limits = c(0.4,0.7), n.breaks = 5, name = "Worst-case lower bound") +
#   theme(panel.background = element_rect(fill = "white", colour = "grey50"),
#         panel.grid.major = element_line(colour = "grey87"))+
#   ggtitle("(C) Condition (E3)") 

plot4 <- ggplot(df4, aes(x = p)) + 
  geom_ribbon(mapping=aes(ymin=lb,ymax=ub), fill="black", colour="white", alpha=0.1)+
  geom_line(aes(y = bound), colour="black") +
  geom_point(aes(y = bound), colour ="black", pch=19, size=2) + 
  scale_x_reverse(n.breaks = 10, name="Overall selection probability") + 
  scale_y_continuous(limits = c(0.4,0.7), n.breaks = 5, name = "Worst-case lower bound") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"))+
  ggtitle("(C) Year 5: Condition (E1)") 

plot5 <- ggplot(df5, aes(x = p)) + 
  geom_ribbon(mapping=aes(ymin=lb,ymax=ub), fill="black", colour="white", alpha=0.1)+
  geom_line(aes(y = bound), colour="black") +
  geom_point(aes(y = bound), colour ="black", pch=19, size=2) + 
  scale_x_reverse(n.breaks = 10, name="Overall selection probability") + 
  scale_y_continuous(limits = c(0.4,0.7), n.breaks = 5, name = "Worst-case lower bound") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"))+
  ggtitle("(D) Year 5: Condition (E2)") 

# plot6 <- ggplot(df6, aes(x = p)) + 
#   geom_ribbon(mapping=aes(ymin=lb,ymax=ub), fill="black", colour="white", alpha=0.1)+
#   geom_line(aes(y = bound), colour="black") +
#   geom_point(aes(y = bound), colour ="black", pch=19, size=2) + 
#   scale_x_reverse(n.breaks = 10, name="Overall selection probability") + 
#   scale_y_continuous(limits = c(0.4,0.7), n.breaks = 5, name = "Worst-case lower bound") +
#   theme(panel.background = element_rect(fill = "white", colour = "grey50"),
#         panel.grid.major = element_line(colour = "grey87"))+
#   ggtitle("(F) Condition (E3)") 

plot <- grid.arrange(plot1, plot2, plot4, plot5, ncol=2)
ggsave(filename = "eg3.eps", plot = plot, device = cairo_ps, 
       width = 10, height = 10) 
## COMPARISON PLOT
# setEPS(width = 12, height = 8); postscript("eg3.eps")
# par(mfrow = c(2,3))
# matplot(bound1[order(-bound1$P),2:3], ylab = "Worst-case bounds of the SAUC(3)", xlab = "Overall selction probability (p)", ylim = c(0.4,0.9), xaxt = "n", type = "b", pch = 1, col = 1, lty = 1)
# abline(h = 0.5, col = "grey", lty = 2)
# axis(1, at = 1:10, labels = p)
# title("(A) Constraint (5.1)", adj = 0, font.main = 1, cex.main = 1.5)
# matplot(bound2[order(-bound2$P),2:3], ylab = "Worst-case bounds of the SAUC(3)", xlab = "Overall selction probability (p)", ylim = c(0.4,0.9), xaxt = "n", type = "b", pch = 2, col = 1, lty = 1)
# abline(h = 0.5, col = "grey", lty = 2)
# axis(1, at = 1:10, labels = p)
# title("(B) Constraint (5.2)", adj = 0, font.main = 1, cex.main = 1.5)
# matplot(bound3[order(-bound3$P),2:3], ylab = "Worst-case bounds of the SAUC(3)", xlab = "Overall selction probability (p)", ylim = c(0.4,0.9), xaxt = "n", type = "b", pch = 5, col = 1, lty = 1)
# abline(h = 0.5, col = "grey", lty = 2)
# axis(1, at = 1:10, labels = p)
# title("(C) Constraint (5.3)", adj = 0, font.main = 1, cex.main = 1.5)
# matplot(bound4[order(-bound1$P),2:3], ylab = "Worst-case bounds of the SAUC(5)", xlab = "Overall selction probability (p)", ylim = c(0.4,0.9), xaxt = "n", type = "b", pch = 1, col = 1, lty = 1)
# abline(h = 0.5, col = "grey", lty = 2)
# axis(1, at = 1:10, labels = p)
# title("(D) Constraint (5.1)", adj = 0, font.main = 1, cex.main = 1.5)
# matplot(bound5[order(-bound2$P),2:3], ylab = "Worst-case bounds of the SAUC(5)", xlab = "Overall selction probability (p)", ylim = c(0.4,0.9), xaxt = "n", type = "b", pch = 2, col = 1, lty = 1)
# abline(h = 0.5, col = "grey", lty = 2)
# axis(1, at = 1:10, labels = p)
# title("(E) Constraint (5.2)", adj = 0, font.main = 1, cex.main = 1.5)
# matplot(bound6[order(-bound3$P),2:3], ylab = "Worst-case bounds of the SAUC(5)", xlab = "Overall selction probability (p)", ylim = c(0.4,0.9), xaxt = "n", type = "b", pch = 5, col = 1, lty = 1)
# abline(h = 0.5, col = "grey", lty = 2)
# axis(1, at = 1:10, labels = p)
# title("(F) Constraint (5.3)", adj = 0, font.main = 1, cex.main = 1.5)
# par(mfrow = c(1,1))
# dev.off()

## TABLE2

p <- seq(1,0.1,-0.1)
tab1 <- cbind.data.frame(p = rep(p,2), 
                         t = c(rep(3, 10), rep(5,10)),
                         A = c(sprintf("%.3f [%.3f, %.3f]", df1[,2], df1[,3], df1[,4]), sprintf("%.3f [%.3f, %.3f]", df4[,2], df4[,3], df4[,4])),
                         B = c(sprintf("%.3f [%.3f, %.3f]", df2[,2], df2[,3], df2[,4]), sprintf("%.3f [%.3f, %.3f]", df5[,2], df5[,3], df5[,4]))
                         # C = c(sprintf("%.3f [%.3f, %.3f]", df3[,2], df3[,3], df3[,4]), sprintf("%.3f [%.3f, %.3f]", df6[,2], df6[,3], df6[,4]))
                         )
colnames(tab1) <- c("p", "Year", "Constraint (D1)", "Constraint (D2)")

sink("tab3.tex")
kbl(tab1, 
    format = "latex",
    longtable = F, 
    booktabs = T, 
    linesep = "",
    digits = 3,
    align = "r",
    escape = FALSE,
    caption = "Example 3: the lower worst-case bound of the SAUC under the assumption of Condition (E1)-(E2) with 95% CI",
    label = "tab3",
    row.names = NA)
sink()
