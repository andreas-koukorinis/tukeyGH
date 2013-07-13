# Tukey's g-and-h transformation
trans.gh <- function(z, g, h){
	
	if(g==0){		
		res <- z * exp((h * z^2) / 2) 			
	}else{
		res <- ((exp(g * z) - 1) / g) * exp((h * z^2 ) / 2)
	}
	
	res
}

# Derivative of Tukey's g-and-h transformation
deriv.gh <- function(z, g, h){
	
	if(g==0){
		res <- (h * z^2 + 1) * exp((h * z^2) / 2)
	}else{
		res <- (((exp(g * z) - 1) / g) * (h * z) + exp(g * z)) * exp((h * z^2) / 2)
	}
	
	res
}

# Inverse of Tukey's g transformation (h=0, !g==0)
inv.trans.g <- function(x, g){
	
	res <- log(x * g + 1) / g
	res			
}
   	
# Standard density: h = 0
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

# Standard distribution: h = 0
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

# Standard quantiles: Not defined for h < 0
std.qgh <- function(p, g, h){
	if(h < 0) stop("h must be non-negative")
	
	z <- qnorm(p)
	res <- trans.gh(z, g, h)
	res
}