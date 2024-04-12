## 
## PLOT FUNCTION
## 

sroc.func.pb <- function(x, par, bias){
  
  u1  <- par[1]
  u2  <- par[2]
  t12 <- par[5]
  t22 <- par[4]
  
  plogis(u1 - (t12/t22) * (qlogis(x) + u2) + bias)
  
}


sauc.value.pb <- function(par, bias){
  
  sauc <- integrate({function (x) sroc.func.pb(x, par, bias)}, 0, 1)
  sauc$value
  
}


plot.sorc.lb<- function(
    df = ldata, dfas = df.as1,
    title = "A. SROC curves under Condition (D4.1)"){
  ggplot(df, aes(fpr, sens)) + geom_point(shape = "x") + 
    scale_x_continuous(limits = c(0,1), name = "1 - Specificity") +
    scale_y_continuous(limits = c(0,1), name = "Sensitivity") +
    geom_function(fun = function(x) sroc.func.pb(x, par = unlist(par), bias = dfas[dfas$p==1,"minb"]), aes(color ="p=1"))+
    geom_function(fun = function(x) sroc.func.pb(x, par = unlist(par), bias = dfas[dfas$p==0.8,"minb"]), aes(color ="p=0.8"))+
    geom_function(fun = function(x) sroc.func.pb(x, par = unlist(par), bias = dfas[dfas$p==0.6,"minb"]), aes(color ="p=0.6"))+
    geom_function(fun = function(x) sroc.func.pb(x, par = unlist(par), bias = dfas[dfas$p==0.4,"minb"]), aes(color ="p=0.4"))+
    geom_function(fun = function(x) sroc.func.pb(x, par = unlist(par), bias = dfas[dfas$p==0.2,"minb"]), aes(color ="p=0.2"))+
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
          panel.grid.major = element_line(colour = "grey87"),
          legend.key = element_rect (fill = "white"),
          legend.position = c(0.75, 0.25), legend.background = element_rect(fill = "white", color = "black"))+
    ggtitle(title) +
    scale_colour_manual("Marginal selection probability",
                        breaks = paste0("p=",c(1, 0.8, 0.6, 0.4, 0.2)),
                        values = c("black", "#e41a1c", "#377eb8", "#4daf4a", "#984ea3"))
}

plot.sauc.comp <- function(
    as.num = 1, 
    title = "A. Lower bounds of the SAUC under Condition (D4.1)"){
  df.k1.as1 <- b.MC.result1[b.MC.result1$as==as.num,]
  k1.saucd1 <- sapply(1:10, function(i) sauc.value.pb(par= unlist(par), bias=df.k1.as1[i,"minb"]) )
  
  df.k2.as1 <- b.MC.result2[b.MC.result2$as==as.num,]
  k2.saucd1 <- sapply(1:10, function(i) sauc.value.pb(par= unlist(par), bias=df.k2.as1[i,"minb"]) )
  
  df.k3.as1 <- b.MC.result3[b.MC.result3$as==as.num,]
  k3.saucd1 <- sapply(1:10, function(i) sauc.value.pb(par= unlist(par), bias=df.k3.as1[i,"minb"]) )
  
  df <- data.frame(p = seq(1, 0.1, -0.1),
                   sauc = rep(fit$sauc.ci[1], 10),
                   k1.saucd1,k2.saucd1,k3.saucd1)
  
  ggplot(df, aes(x=p)) + 
    geom_line(aes(y = sauc, color = "Estimate without PB"), lty=2, size=1) +
    geom_line(aes(y = k1.saucd1, color = "Simulation-based worst-case bound (K=1000)"), lty=1, size=1) +
    geom_line(aes(y = k2.saucd1, color = "Simulation-based worst-case bound (K=2000)"), lty=2, size=1) +
    geom_line(aes(y = k3.saucd1, color = "Simulation-based worst-case bound (K=20000)"),lty=3, size=1) +
    scale_x_reverse(n.breaks = 10, name="Marginal selection probability") +
    scale_y_continuous(limits = c(0,1), name = "SAUC") +
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
          panel.grid.major = element_line(colour = "grey87"),
          legend.key = element_rect (fill = "white"),
          legend.position = c(0.4, 0.25), 
          legend.background = element_rect(fill = "white", color = "black"),
          legend.title = element_blank())+
    scale_colour_manual("",
                        breaks = c("Simulation-based worst-case bound (K=1000)",
                                   "Simulation-based worst-case bound (K=2000)",
                                   "Simulation-based worst-case bound (K=20000)",
                                   "Estimate without PB"),
                        values = c("#808800","#dc143c","#60100b","grey50"),
                        guide = guide_legend(override.aes = list(lty = c(1,2,3,2))))+
    ggtitle(title)
}
