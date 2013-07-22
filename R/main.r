#' Generate random numbers from Tukey's g-and-h
#' 
#' @param n n of samples to generate
#' @param g parameter
#' @param h parameter
#' @param A parameter
#' @param B parameter
rgh <- function(n, g=0, h=0, A=0, B=1){
	if(!B > 0) stop("B must be positive")
	
	z <- rnorm(n)
	res <- A + B * trans.gh(z, g, h)
	res
}

#' Quantile function for Tukey's g-and-h
#' Currently working for h >= 0 only
#'
#' @param p probs in (0, 1)
#' @param g parameter
#' @param h parameter
#' @param A parameter
#' @param B parameter
qgh <- function(p, g=0, h=0, A=0, B=1){
	if(!B > 0) stop("B must be positive")
	
	res <- A + B * std.qgh(p, g, h)
	res
}

#' Distribution function for Tukey's g-and-h
#' Currently working for h=0 only
#'
#' @param x numeric
#' @param g parameter
#' @param h parameter
#' @param A parameter
#' @param B parameter
pgh <- function(x, g=0, A=0, B=1){
	if(!B > 0) stop("B must be positive")
	
	z <- (x - A) / B
	res <- std.pgh(z, g)
	res
}

#' Density function for Tukey's g-and-h
#' Currently working for h=0 only
#'
#' @param x numeric
#' @param g parameter
#' @param h parameter
#' @param A parameter
#' @param B parameter
dgh <- function(x, g=0, h=0, A=0, B=1){
	if(!B > 0) stop("B must be positive")
	
	z <- (x - A) / B
	res <- (1 / B) * std.dgh(z, g)
	res
}

#' Fit g-and-h via quatile method
#'
#' @param smpl sample vector
#' @param verbose show intermediate tables
qFit <- function(smpl, verbose=FALSE){
	# compute quantiles for fitting
	q.tabl <- quantile.table(smpl)	
	if (verbose) print(q.tabl)

	# fit A and g 
	g.tabl <- g.table(q.tabl)
	if (verbose) print(g.tabl)

	# fit B and h
	h.tabl <- h.table(g.tabl)
	if (verbose) print(h.tabl)

	res <- c(g=g.tabl$g, h=h.tabl$h, A=g.tabl$A, B=h.tabl$B)
	res	
}


#' Calculate goodness of fit
#'
#' @param smpl sample vector
#' @param fit.pars results from qFit()
#' @param plot.fit Do you want to show fit plots?
goodnessFit <- function(smpl, fit.pars, plot.fit=FALSE){
  p <- c( .01, .025, .05, .1, .25, .4)
  p <- sort(c(1 - p, p))
  
  A <- fit.pars["A"]
  B <- fit.pars["B"]
  g <- fit.pars["g"]
  h <- fit.pars["h"]
  
  x <- quantiles(smpl, p)
  xfit <- qgh(p, g, h, A, B)
  
  if (plot.fit){
  plot(x, xfit, pch=19, cex=0.3)
  abline(0, 1, col="red")
  
  plot(p, ((x - xfit) / x)*100 , pch=19, cex=0.8, ylab="est. error %")
  abline(h=0, col="red")
  }
  
  res <- abs((f(p) - q_gh(p, g, h, A, B)) / f(p)) * 100 
  
  list(tab = data.frame(p=p, target=x, fitted=xfit, err.Pcent=round(res, 3)),
       overall = mean(res))
  
}