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