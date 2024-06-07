## 
## PLOT FUNCTION
## 

.sroc.func.pb <- function(x, par, bias){
  
  u1  <- par[1]
  u2  <- par[2]
  t12 <- par[5]
  t22 <- par[4]
  
  plogis(u1 - (t12/t22) * (qlogis(x) + u2) + bias)
  
}

.sroc.func.pb.par <- function(x, par){
  
  u1  <- par[1]
  u2  <- par[2]
  t12 <- par[5]
  t22 <- par[4]
  
  plogis(u1 - (t12/t22) * (qlogis(x) + u2))
  
}


.sauc.value.pb <- function(par, bias){
  
  if(is.na(bias)) NA else {
  sauc <- integrate({function (x) .sroc.func.pb(x, par, bias)}, 0, 1)
  sauc$value}
  
}

.plot.sroc.combine <- function(
    ldata = example2.2, dfas = df[df$as==1,], zhou.par = zhou.par1,
    title = "A. SROC curves under Condition (D4.1)"){
  
  ggplot(ldata, aes(fpr, sens)) + geom_point(shape = "x") + 
    scale_x_continuous(limits = c(0,1), name = "1 - Specificity") +
    scale_y_continuous(limits = c(0,1), name = "Sensitivity") +
    geom_function(fun = function(x) .sroc.func.pb(x, par = unlist(par), bias = dfas[dfas$p==1,"minb"]), aes(color ="p=1"))+
    geom_function(fun = function(x) .sroc.func.pb(x, par = unlist(par), bias = dfas[dfas$p==0.8,"minb"]), aes(color ="p=0.8"))+
    geom_function(fun = function(x) .sroc.func.pb(x, par = unlist(par), bias = dfas[dfas$p==0.6,"minb"]), aes(color ="p=0.6"))+
    geom_function(fun = function(x) .sroc.func.pb(x, par = unlist(par), bias = dfas[dfas$p==0.4,"minb"]), aes(color ="p=0.4"))+
    geom_function(fun = function(x) .sroc.func.pb(x, par = unlist(par), bias = dfas[dfas$p==0.2,"minb"]), aes(color ="p=0.2"))+
    geom_function(fun = function(x) .sroc.func.pb.par(x, par = zhou.par[,1]), aes(color ="p=1"), linetype=2)+
    geom_function(fun = function(x) .sroc.func.pb.par(x, par = zhou.par[,3]), aes(color ="p=0.8"), linetype=2)+
    geom_function(fun = function(x) .sroc.func.pb.par(x, par = zhou.par[,5]), aes(color ="p=0.6"), linetype=2)+
    geom_function(fun = function(x) .sroc.func.pb.par(x, par = zhou.par[,7]), aes(color ="p=0.4"), linetype=2)+
    geom_function(fun = function(x) .sroc.func.pb.par(x, par = zhou.par[,9]), aes(color ="p=0.2"), linetype=2)+
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
          panel.grid.major = element_line(colour = "grey87"),
          legend.key = element_rect (fill = "white"),
          legend.position = c(0.75, 0.25), legend.background = element_rect(fill = "white", color = "black"))+
    ggtitle(title) +
    scale_colour_manual("Marginal selection probability",
                        breaks = paste0("p=",c(1, 0.8, 0.6, 0.4, 0.2)),
                        values = c("black", "#e41a1c", "#377eb8", "#4daf4a", "#984ea3"))
  
}

.plot.sorc.npar <- function(
    ldata = example2.2, mle = par, dfas = df[df$as==1,],
    title = "A. SROC curves under Condition (D4.1)"){
  
  ggplot(ldata, aes(fpr, sens)) + geom_point(shape = "x") + 
    scale_x_continuous(limits = c(0,1), name = "1 - Specificity") +
    scale_y_continuous(limits = c(0,1), name = "Sensitivity") +
    geom_function(fun = function(x) .sroc.func.pb(x, par = unlist(mle), bias = dfas[dfas$p==1,"minb"]), aes(color ="p=1"))+
    geom_function(fun = function(x) .sroc.func.pb(x, par = unlist(mle), bias = dfas[dfas$p==0.8,"minb"]), aes(color ="p=0.8"))+
    geom_function(fun = function(x) .sroc.func.pb(x, par = unlist(mle), bias = dfas[dfas$p==0.6,"minb"]), aes(color ="p=0.6"))+
    geom_function(fun = function(x) .sroc.func.pb(x, par = unlist(mle), bias = dfas[dfas$p==0.4,"minb"]), aes(color ="p=0.4"))+
    geom_function(fun = function(x) .sroc.func.pb(x, par = unlist(mle), bias = dfas[dfas$p==0.2,"minb"]), aes(color ="p=0.2"))+
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
          panel.grid.major = element_line(colour = "grey87"),
          legend.key = element_rect (fill = "white"),
          legend.position = c(0.75, 0.25), legend.background = element_rect(fill = "white", color = "black"))+
    ggtitle(title) +
    scale_colour_manual("Marginal selection probability",
                        breaks = paste0("p=",c(1, 0.8, 0.6, 0.4, 0.2)),
                        values = c("black", "#e41a1c", "#377eb8", "#4daf4a", "#984ea3"))
  
}

.plot.sorc.par <- function(
    ldata = example2.2, zhou.par = zhou.par1,
    title = "A. SROC curves by Zhou's method"){
  
  ggplot(ldata, aes(fpr, sens)) + geom_point(shape = "x") + 
    scale_x_continuous(limits = c(0,1), name = "1 - Specificity") +
    scale_y_continuous(limits = c(0,1), name = "Sensitivity") +
    geom_function(fun = function(x) .sroc.func.pb.par(x, par = zhou.par[,1]), aes(color ="p=1"), linetype=1)+
    geom_function(fun = function(x) .sroc.func.pb.par(x, par = zhou.par[,3]), aes(color ="p=0.8"), linetype=1)+
    geom_function(fun = function(x) .sroc.func.pb.par(x, par = zhou.par[,5]), aes(color ="p=0.6"), linetype=1)+
    geom_function(fun = function(x) .sroc.func.pb.par(x, par = zhou.par[,7]), aes(color ="p=0.4"), linetype=1)+
    geom_function(fun = function(x) .sroc.func.pb.par(x, par = zhou.par[,9]), aes(color ="p=0.2"), linetype=1)+
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
          panel.grid.major = element_line(colour = "grey87"),
          legend.key = element_rect (fill = "white"),
          legend.position = c(0.75, 0.25), legend.background = element_rect(fill = "white", color = "black"))+
    ggtitle(title) +
    scale_colour_manual("Marginal selection probability",
                        breaks = paste0("p=",c(1, 0.8, 0.6, 0.4, 0.2)),
                        values = c("black", "#e41a1c", "#377eb8", "#4daf4a", "#984ea3"))
  
}

