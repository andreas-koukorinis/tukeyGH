#' Tukey's g-and-h transformation
#'
#' @param x numeric vector
#' @param g g parameter
#' @param h h parameter
trans.gh <- function(x, g, h){
	
	if(g==0){		
		res <- x * exp((h * x^2) / 2) 			
	}else{
		res <- ((exp(g * x) - 1) / g) * exp((h * x^2 ) / 2)
	}
	
	res
}


#' Derivative of Tukey's g-and-h transformation
#'
#' @param x numeric vector
#' @param g g parameter
#' @param h h parameter
deriv.gh <- function(x, g, h){
	
	if(g==0){
		res <- (h * x^2 + 1) * exp((h * x^2) / 2)
	}else{
		res <- (((exp(g * x) - 1) / g) * (h * x) + exp(g * x)) * exp((h * x^2) / 2)
	}
	
	res
}

#' Inverse of Tukey's g transformation (When h=0 and g is not 0)
#' When g=0 there is no need for an inverse since T(x)=x
#'
#' @param x numeric vector
#' @param g g parameter
inv.trans.g <- function(x, g){
	
	res <- log(x * g + 1) / g
	res			
}
   	
#' Standard density: h = 0 (Shifted logNormal case)
#'
#' @param x numeric vector
#' @param g g parameter
std.dgh <- function(x, g){

    if(g==0){
    	res <- dnorm(x)
    }else{
    	support <- (x * g + 1) > 0
    	res <- rep(0, length(support))  
    	x.inv <- inv.trans.g(x[support], g)
    	res[support] <- abs(1 / deriv.gh(x.inv, g, h=0)) * dnorm(x.inv)
    }
   
    res    	
}

#' Standard distribution: h = 0 (Shifted logNormal case)
#'
#' @param x numeric vector
#' @param g g parameter
std.pgh <- function(x, g){

    if(g==0){
    	res <- pnorm(x)
    }else{
    	if(g > 0){
    		support <- (x * g + 1) > 0
    		res <- rep(0, length(support))  
    		x.inv <- inv.trans.g(x[support], g)
    		res[support] <- pnorm(x.inv)
    	}else{
    		support <- (x * g + 1) > 0
    		res <- rep(1, length(support))  
    		x.inv <- inv.trans.g(x[support], g)
    		res[support] <- pnorm(x.inv)
        }
    }
    
    res    	
}

#' Standard quantiles: Not defined for h < 0
#'
#' @param p probs in (0, 1)
#' @param g parameter
#' @param h parameter (h must be non-negative)
std.qgh <- function(p, g, h){
	if(h < 0) stop("h must be non-negative")
	
	z <- qnorm(p)
	res <- trans.gh(z, g, h)
	res
}