## NOTES
## DATA IS SAVED IN "eg2-IVD.csv"
## 

library(kableExtra)
library(mada)
# devtools::install_github("meta2020/dtametasa")
library(dtametasa)
require(gridExtra)

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

## SAUC
p <- seq(1,0.1,-0.1)
## Zhou et al. METHOD 1
sauc1 <- sapply(p, function(p) dtametasa.fc(data=example2, p=p, correct.type = "all", c1.square = 1)$sauc.ci)
sauc1[2,2] <- (sauc1[2,1]+sauc1[2,3])/2
sauc1[3,2] <- (sauc1[3,1]+sauc1[3,3])/2
## Zhou et al. METHOD 2
sauc2 <- sapply(p, function(p) dtametasa.fc(data=example2, p=p, correct.type = "all", c1.square = 0)$sauc.ci)
## Zhou et al. METHOD 1
sauc3 <- sapply(p, function(p) dtametasa.fc(data=example2, p=p, correct.type = "all", c1.square = 0.5)$sauc.ci)

## BOUNDS UNDER CONDITION D1-D3
bound1 <- read.csv("MCbound-saucD1.csv")
bound2 <- read.csv("MCbound-saucD2.csv")
bound3 <- read.csv("MCbound-saucD3.csv")


ma21<- dtametasa.fc(data=example2, p=1, correct.type = "all")
l1 <- ma21$sauc.ci[2]-ma21$sauc.ci[1]
l2 <- ma21$sauc.ci[3]-ma21$sauc.ci[1]

## Under CONDITION D1-D3
## PLOT 2

df1 <- data.frame(p=p, bound = rev(bound1$SAUC_LB), lb = rev(bound1$SAUC_LB+l1), ub = rev(bound1$SAUC_LB+l2))
df2 <- data.frame(p=p, bound = rev(bound2$SAUC_LB), lb = rev(bound2$SAUC_LB+l1), ub = rev(bound2$SAUC_LB+l2))
df3 <- data.frame(p=p, bound = rev(bound3$SAUC_LB), lb = rev(bound3$SAUC_LB+l1), ub = rev(bound3$SAUC_LB+l2))

plot1 <- ggplot(df1, aes(x = p)) + 
  geom_ribbon(mapping=aes(ymin=lb,ymax=ub), fill="black", colour="white", alpha=0.1)+
  geom_line(aes(y = bound), colour="black") +
  geom_point(aes(y = bound), colour ="black", pch=19, size=2) + 
  scale_x_reverse(n.breaks = 10, name="Overall selection probability") + 
  scale_y_continuous(limits = c(0.45,1), n.breaks = 7, name = "Worst-case lower bound") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"))+
  ggtitle("(A) Condition (D1)") 

plot2 <- ggplot(df2, aes(x = p)) + 
  geom_ribbon(mapping=aes(ymin=lb,ymax=ub), fill="black", colour="white", alpha=0.1)+
  geom_line(aes(y = bound), colour="black") +
  geom_point(aes(y = bound), colour ="black", pch=19, size=2) + 
  scale_x_reverse(n.breaks = 10, name="Overall selection probability") + 
  scale_y_continuous(limits = c(0.45,1), n.breaks = 7, name = "Worst-case lower bound") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"))+
  ggtitle("(B) Condition (D2)") 

plot3 <- ggplot(df3, aes(x = p)) + 
  geom_ribbon(mapping=aes(ymin=lb,ymax=ub), fill="black", colour="white", alpha=0.1)+
  geom_line(aes(y = bound), colour="black") +
  geom_point(aes(y = bound), colour ="black", pch=19, size=2) + 
  scale_x_reverse(n.breaks = 10, name="Overall selection probability") + 
  scale_y_continuous(limits = c(0.45,1), n.breaks = 7, name = "Worst-case lower bound") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "grey87"))+
  ggtitle("(C) Condition (D3)") 
plot <- grid.arrange(plot1, plot2, plot3, ncol=3)
ggsave(filename = "eg2.eps", plot = plot, device = cairo_ps, 
       width = 12, height = 4) 

## COMPARISON PLOT
# setEPS(width = 12, height = 4); postscript("eg2.eps")
# par(mfrow = c(1,3))
# matplot(bound1[order(-bound1$P),2:3], ylab = "Worst-case bounds of the SAUC", xlab = "Overall selction probability (p)", ylim = c(0.5,1), xaxt = "n", type = "b", pch = 1, col = 1, lty = 1)
# abline(h = 0.5, col = "grey", lty = 2)
# axis(1, at = 1:10, labels = p)
# title("(A) Constraint (4.1)", adj = 0, font.main = 1, cex.main = 1.5)
# matplot(bound2[order(-bound2$P),2:3], ylab = "Worst-case bounds of the SAUC", xlab = "Overall selction probability (p)", ylim = c(0.5,1), xaxt = "n", type = "b", pch = 2, col = 1, lty = 1)
# abline(h = 0.5, col = "grey", lty = 2)
# axis(1, at = 1:10, labels = p)
# title("(B) Constraint (4.2)", adj = 0, font.main = 1, cex.main = 1.5)
# matplot(bound3[order(-bound3$P),2:3], ylab = "Worst-case bounds of the SAUC", xlab = "Overall selction probability (p)", ylim = c(0.5,1), xaxt = "n", type = "b", pch = 5, col = 1, lty = 1)
# abline(h = 0.5, col = "grey", lty = 2)
# axis(1, at = 1:10, labels = p)
# title("(C) Constraint (4.3)", adj = 0, font.main = 1, cex.main = 1.5)
# par(mfrow = c(1,1))
# dev.off()

## TABLE2


tab1 <- cbind.data.frame(p = p, 
                         A = sprintf("%.3f [%.3f, %.3f]", df1[,2], df1[,3], df1[,4]), 
                         B = sprintf("%.3f [%.3f, %.3f]", df2[,2], df2[,3], df2[,4]),
                         C = sprintf("%.3f [%.3f, %.3f]", df3[,2], df3[,3], df3[,4]))
colnames(tab1) <- c("p", "Condition (D1)", "Condition (D2)", "Condition (D3)")

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
