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

# .plot.sroc.combine <- function(
#     ldata = example2.2, dfas = df[df$as==1,], zhou.par = zhou.par1,
#     title = "A. SROC curves under Condition (D4.1)"){
#   
#   ggplot(ldata, aes(fpr, sens)) + geom_point(shape = "x") + 
#     scale_x_continuous(limits = c(0,1), name = "1 - Specificity") +
#     scale_y_continuous(limits = c(0,1), name = "Sensitivity") +
#     geom_function(fun = function(x) .sroc.func.pb(x, par = unlist(par), bias = dfas[dfas$p==1,"minb"]), aes(color ="p=1"))+
#     geom_function(fun = function(x) .sroc.func.pb(x, par = unlist(par), bias = dfas[dfas$p==0.8,"minb"]), aes(color ="p=0.8"))+
#     geom_function(fun = function(x) .sroc.func.pb(x, par = unlist(par), bias = dfas[dfas$p==0.6,"minb"]), aes(color ="p=0.6"))+
#     geom_function(fun = function(x) .sroc.func.pb(x, par = unlist(par), bias = dfas[dfas$p==0.4,"minb"]), aes(color ="p=0.4"))+
#     geom_function(fun = function(x) .sroc.func.pb(x, par = unlist(par), bias = dfas[dfas$p==0.2,"minb"]), aes(color ="p=0.2"))+
#     geom_function(fun = function(x) .sroc.func.pb.par(x, par = zhou.par[,1]), aes(color ="p=1"), linetype=2)+
#     geom_function(fun = function(x) .sroc.func.pb.par(x, par = zhou.par[,3]), aes(color ="p=0.8"), linetype=2)+
#     geom_function(fun = function(x) .sroc.func.pb.par(x, par = zhou.par[,5]), aes(color ="p=0.6"), linetype=2)+
#     geom_function(fun = function(x) .sroc.func.pb.par(x, par = zhou.par[,7]), aes(color ="p=0.4"), linetype=2)+
#     geom_function(fun = function(x) .sroc.func.pb.par(x, par = zhou.par[,9]), aes(color ="p=0.2"), linetype=2)+
#     theme(panel.background = element_rect(fill = "white", colour = "grey50"),
#           panel.grid.major = element_line(colour = "grey87"),
#           legend.key = element_rect (fill = "white"),
#           legend.position = c(0.75, 0.15), legend.background = element_rect(fill = "white", color = "black"))+
#     ggtitle(title) +
#     scale_colour_manual("Marginal selection probability",
#                         breaks = paste0("p=",c(1, 0.8, 0.6, 0.4, 0.2)),
#                         values = c("black", "#e41a1c", "#377eb8", "#4daf4a", "#984ea3"))
#   
# }

## FUNCTION TO PLOT THE SROC OF THE NP METHOD
.plot.sorc.npar <- function(
    ldata = example2.2, mle = par, dfas = df[df$as==1,],
    title = "A. SROC curves under Condition (D4.1)"){
  
  ggplot(ldata, aes(fpr, sens)) + geom_point(shape = "x") + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    scale_x_continuous(limits = c(0,1), name = "1 - Specificity") +
    scale_y_continuous(limits = c(0,1), name = "Sensitivity") +
    geom_function(fun = function(x) .sroc.func.pb(x, par = unlist(mle), bias = dfas[dfas$p==1,"minb"]), aes(color ="p=1"), linewidth = 1)+
    geom_function(fun = function(x) .sroc.func.pb(x, par = unlist(mle), bias = dfas[dfas$p==0.8,"minb"]), aes(color ="p=0.8"), linewidth = 1)+
    geom_function(fun = function(x) .sroc.func.pb(x, par = unlist(mle), bias = dfas[dfas$p==0.6,"minb"]), aes(color ="p=0.6"), linewidth = 1)+
    geom_function(fun = function(x) .sroc.func.pb(x, par = unlist(mle), bias = dfas[dfas$p==0.4,"minb"]), aes(color ="p=0.4"), linewidth = 1)+
    geom_function(fun = function(x) .sroc.func.pb(x, par = unlist(mle), bias = dfas[dfas$p==0.2,"minb"]), aes(color ="p=0.2"), linewidth = 1)+
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
          panel.grid.major = element_line(colour = "grey87"),
          legend.key = element_rect (fill = "white"),
          legend.position = c(0.7, 0.2), 
          legend.background = element_rect(fill = "white", color = "black"),
          legend.text=element_text(size=15),
          legend.title = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          plot.title = element_text(size = 15))+
    ggtitle(title) +
    scale_colour_manual("Marginal selection probability",
                        breaks = paste0("p=",c(1, 0.8, 0.6, 0.4, 0.2)),
                        values = c("black", "#e41a1c", "#377eb8", "#4daf4a", "#984ea3"))
  
}

## FUNCTION TO PLOT THE SROC OF THE P METHOD
.plot.sorc.par <- function(
    ldata = example2.2, zhou.par = zhou.par1,
    title = "A. SROC curves by Zhou's method"){
  
  ggplot(ldata, aes(fpr, sens)) + geom_point(shape = "x") + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    scale_x_continuous(limits = c(0,1), name = "1 - Specificity") +
    scale_y_continuous(limits = c(0,1), name = "Sensitivity") +
    geom_function(fun = function(x) .sroc.func.pb.par(x, par = zhou.par[,1]), aes(color ="p=1"), linewidth=1)+
    geom_function(fun = function(x) .sroc.func.pb.par(x, par = zhou.par[,3]), aes(color ="p=0.8"), linewidth=1)+
    geom_function(fun = function(x) .sroc.func.pb.par(x, par = zhou.par[,5]), aes(color ="p=0.6"), linewidth=1)+
    geom_function(fun = function(x) .sroc.func.pb.par(x, par = zhou.par[,7]), aes(color ="p=0.4"), linewidth=1)+
    geom_function(fun = function(x) .sroc.func.pb.par(x, par = zhou.par[,9]), aes(color ="p=0.2"), linewidth=1)+
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
          panel.grid.major = element_line(colour = "grey87"),
          legend.key = element_rect (fill = "white"),
          legend.position = c(0.7, 0.2), 
          legend.background = element_rect(fill = "white", color = "black"),
          legend.text=element_text(size=15),
          legend.title = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          plot.title = element_text(size = 15))+
    ggtitle(title) +
    scale_colour_manual("Marginal selection probability",
                        breaks = paste0("p=",c(1, 0.8, 0.6, 0.4, 0.2)),
                        values = c("black", "#e41a1c", "#377eb8", "#4daf4a", "#984ea3"))
  
}

## FUNCTION OT PLOT THE SAUC OF ALL METHODS
.plot.sauc.all <- function(
    dfc = dfc1, ldg=ldg2a, 
    title = "G. Bounds for the SAUC given PB from sensitivity"){
  
  ggplot(dfc, aes(x = p, y = est, color = grp, linetype = grp)) +
    geom_line(linewidth=1.2) +
    scale_y_continuous(limits = c(0,1),n.breaks = 10, name = "SAUC") +
    scale_x_reverse(n.breaks = 10, name="Marginal selection probability") +
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
          panel.grid.major = element_line(colour = "grey87"),
          legend.key = element_rect (fill = "white"),
          legend.position = c(0.5, 0.15),
          legend.background = element_rect(fill = "white", color = "black"),
          legend.title = element_blank(),
          legend.text=element_text(size=15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          plot.title = element_text(size = 15))+
    scale_colour_manual("",
                        breaks = ldg,
                        values = c("grey50","#e41a1c","#377eb8","#377eb8"))+
    scale_linetype_manual("",
                          breaks = ldg,
                          values = c(2,1,1,3))+
    ggtitle(title)
  
}

.plot.sauc.single <- function(
    dfc=dfc1, ldg=ldg2a,
    title="D. Simulation-based lower bounds under Condition (D4.1)"){
  
  ggplot(dfc, aes(x = p, y = est, color = grp, linetype = grp)) +
    geom_line(linewidth=1.2) +
    scale_y_continuous(limits = c(0,1),n.breaks = 10, name = "SAUC") +
    scale_x_reverse(n.breaks = 10, name="Marginal selection probability") + 
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
          panel.grid.major = element_line(colour = "grey87"),
          legend.key = element_rect (fill = "white"),
          legend.position = c(0.5, 0.15), 
          legend.background = element_rect(fill = "white", color = "black"),
          legend.title = element_blank(),
          legend.text=element_text(size=15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          plot.title = element_text(size = 15))+
    scale_colour_manual("",
                        breaks = ldg,
                        values = c("grey50","#e41a1c"))+
    scale_linetype_manual("",
                          breaks = ldg,
                          values = c(2,1))+
    ggtitle(title)
  
}


.plot.kest <- function(
    dfk=dfk2a, lgd=lgd3, ylim = c(0.4,0.75),
    title="A. Lower bounds of the SAUC under Condition (D4.1)"){
  
  ggplot(dfk, aes(x = p, y = est, color = grp, linetype = grp)) +
    geom_line(size=1) +
    scale_y_continuous(limits = range(ylim), breaks = ylim, name = "SAUC") +
    scale_x_reverse(n.breaks = 10, name="Marginal selection probability") + 
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
          panel.grid.major = element_line(colour = "grey87"),
          legend.key = element_rect (fill = "white"),
          legend.position = c(0.3, 0.25), 
          legend.background = element_rect(fill = "white", color = "black"),
          legend.title = element_blank(),
          legend.text=element_text(size=15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          plot.title = element_text(size = 15))+
    scale_colour_manual("",
                        breaks = lgd,
                        values = c("#4daf4a","#e41a1c","#984ea3","#ff7f00"))+
    scale_linetype_manual("",
                          breaks = lgd,
                          values = c(1,1,2,2))+
    ggtitle(title)
  
}

.plot.kbox <- function(
    dfk=dfr2a, lgd=lgd3, ylim = seq(0.4,0.75, 0.05),
    title="D. 10-time lower bounds of the SAUC under Condition (D4.1)"){
  
  ggplot(dfk, aes(x=pp, y=est, fill=grp2)) + 
    geom_boxplot(alpha = 0.5)+
    scale_x_discrete(limits = rev(levels(dfk$pp)),name="Marginal selection probability") +
    scale_y_continuous(limits = range(ylim), breaks = ylim, name = "SAUC") +
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
          panel.grid.major = element_line(colour = "grey87"),
          legend.key = element_rect (fill = "white"),
          legend.position = c(0.3, 0.25), legend.background = element_rect(fill = "white", color = "black"),
          legend.title = element_blank(),
          legend.text=element_text(size=15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          plot.title = element_text(size = 14))+
    scale_fill_manual("",
                      breaks = lgd,
                      values = c("#4daf4a","#e41a1c","#984ea3","#ff7f00")) +
    ggtitle(title)
}
