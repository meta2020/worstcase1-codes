## 
## PLOT FUNCTION
## 

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
  
  plogis(u1+b1 - (t12/t22) * (qlogis(x) + b2 + u2) )
  
}

sauc.value <- function(par, bias){
  
  sauc <- integrate({function (x) sroc.func(x, par, bias)}, 0, 1)
  sauc$value
  
}

## PLOT ALL THE 10 SIMULATION-BASD BOUNDS
## 
plot10 <- function(y3d1_wide = y3d1_wide, 
                   title = "A. SROC curves under Condition (D1)",
                   spd1 = spd1, sed1 = sed1){
  
  ggplot(df, aes(x, y)) + geom_point(shape = "x") + 
    scale_x_continuous(limits = c(0,1), name = "1 - Specificity") +
    scale_y_continuous(limits = c(0,1), name = "Sensitivity") +
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==1,"sim1"]), aes(color ="p=1"), alpha=0.2)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==1,"sim2"]), aes(color ="p=1"), alpha=0.2)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==1,"sim3"]), aes(color ="p=1"), alpha=0.2)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==1,"sim4"]), aes(color ="p=1"), alpha=0.2)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==1,"sim5"]), aes(color ="p=1"), alpha=0.2)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==1,"sim6"]), aes(color ="p=1"), alpha=0.2)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==1,"sim7"]), aes(color ="p=1"), alpha=0.2)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==1,"sim8"]), aes(color ="p=1"), alpha=0.2)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==1,"sim9"]), aes(color ="p=1"), alpha=0.2)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==1,"sim10"]), aes(color ="p=1"), alpha=0.2)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.8,"sim1"]), aes(color ="p=0.8"), alpha=0.2)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.8,"sim2"]), aes(color ="p=0.8"), alpha=0.2)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.8,"sim3"]), aes(color ="p=0.8"), alpha=0.2)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.8,"sim4"]), aes(color ="p=0.8"), alpha=0.2)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.8,"sim5"]), aes(color ="p=0.8"), alpha=0.2)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.8,"sim6"]), aes(color ="p=0.8"), alpha=0.2)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.8,"sim7"]), aes(color ="p=0.8"), alpha=0.2)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.8,"sim8"]), aes(color ="p=0.8"), alpha=0.2)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.8,"sim9"]), aes(color ="p=0.8"), alpha=0.2)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.8,"sim10"]), aes(color ="p=0.8"), alpha=0.2)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.6,"sim1"]), aes(color ="p=0.6"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.6,"sim2"]), aes(color ="p=0.6"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.6,"sim3"]), aes(color ="p=0.6"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.6,"sim4"]), aes(color ="p=0.6"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.6,"sim5"]), aes(color ="p=0.6"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.6,"sim6"]), aes(color ="p=0.6"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.6,"sim7"]), aes(color ="p=0.6"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.6,"sim8"]), aes(color ="p=0.6"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.6,"sim9"]), aes(color ="p=0.6"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.6,"sim10"]), aes(color ="p=0.6"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.4,"sim1"]), aes(color ="p=0.4"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.4,"sim2"]), aes(color ="p=0.4"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.4,"sim3"]), aes(color ="p=0.4"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.4,"sim4"]), aes(color ="p=0.4"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.4,"sim5"]), aes(color ="p=0.4"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.4,"sim6"]), aes(color ="p=0.4"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.4,"sim7"]), aes(color ="p=0.4"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.4,"sim8"]), aes(color ="p=0.4"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.4,"sim9"]), aes(color ="p=0.4"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.4,"sim10"]), aes(color ="p=0.4"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.2,"sim1"]), aes(color ="p=0.2"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.2,"sim2"]), aes(color ="p=0.2"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.2,"sim3"]), aes(color ="p=0.2"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.2,"sim4"]), aes(color ="p=0.2"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.2,"sim5"]), aes(color ="p=0.2"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.2,"sim6"]), aes(color ="p=0.2"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.2,"sim7"]), aes(color ="p=0.2"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.2,"sim8"]), aes(color ="p=0.2"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.2,"sim9"]), aes(color ="p=0.2"), alpha=0.3)+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1_wide[y3d1_wide$p==0.2,"sim10"]), aes(color ="p=0.2"), alpha=0.3)+
    # geom_point(data = data.frame(x = 1-spd1[c(10,9,7,5,3)], y = sed1[c(10,9,7,5,3)]), aes(color = paste0("p=",c(1, 0.9, 0.7, 0.5, 0.3))))+
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
          panel.grid.major = element_line(colour = "grey87"),
          legend.key = element_rect (fill = "white"),
          legend.position = c(0.75, 0.2), legend.background = element_rect(fill = "white", color = "black"))+
    ggtitle(title) +
    scale_colour_manual(TeX("Marginal selection probability"),
                        breaks = c(paste0("p=",c(1, 0.8, 0.6, 0.4, 0.2))),
                        values = c("black", "#e41a1c", "#377eb8", "#4daf4a", "#984ea3"))
}



plot.ave <- function(y3d1 = y3d1, 
                     title = "A. SROC curves under Condition (D1)",
                     spd1 = spd1, sed1 = sed1){
  
  ggplot(df, aes(x, y)) + geom_point(shape = "x") + 
    scale_x_continuous(limits = c(0,1), name = "1 - Specificity") +
    scale_y_continuous(limits = c(0,1), name = "Sensitivity") +
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1[y3d1$p==1,"x"]), aes(color ="p=1"))+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1[y3d1$p==0.8,"x"]), aes(color ="p=0.8"))+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1[y3d1$p==0.6,"x"]), aes(color ="p=0.6"))+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1[y3d1$p==0.4,"x"]), aes(color ="p=0.4"))+
    geom_function(fun = function(x) sroc.func(x, par = unlist(par), bias = y3d1[y3d1$p==0.2,"x"]), aes(color ="p=0.2"))+
    # geom_point(data = data.frame(x = 1-spd1[c(10,8,6,4,2)], y = sed1[c(10,8,6,4,2)]), aes(color = paste0("p=",c(1, 0.8, 0.6, 0.4, 0.2))))+
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
          panel.grid.major = element_line(colour = "grey87"),
          legend.key = element_rect (fill = "white"),
          legend.position = c(0.75, 0.25), legend.background = element_rect(fill = "white", color = "black"))+
    ggtitle(title) +
    scale_colour_manual("Marginal selection probability",
                        breaks = paste0("p=",c(1, 0.8, 0.6, 0.4, 0.2)),
                        values = c("black", "#e41a1c", "#377eb8", "#4daf4a", "#984ea3"))
}
