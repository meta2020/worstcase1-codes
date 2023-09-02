##******************************************************************************
##
## PARAMETER IN G DISTRIBUTON
##
##******************************************************************************


cens.eta <- function(data, med.year, n1, n0, s1_med, s0_med, eta.range = c(0,1), init.eta = 0.01){
	
	dt.mct <- data.frame(n1 = eval(substitute(n1), data), 
																						n0 = eval(substitute(n0), data),
																						med.year = eval(substitute(med.year), data),
																						s1_med = eval(substitute(s1_med), data),
																						s0_med = eval(substitute(s0_med), data))
	
	S_tf     <-  with(dt.mct, (s0_med * n0 + s1_med * n1) / (n0+n1) )
	
	U <- function(eta) sum( (0.5 - S_tf * exp(- dt.mct$med.year * eta))^2 )  ## G = P(C>t) =  exp(- eta t)
	
	int.eta <- c(eta = init.eta)
	
	nlminb(int.eta, U, lower = eta.range[1], upper = eta.range[2])
	
	
}
