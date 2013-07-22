#' Depth function
#'
#' Internal function
#' Given a sample size n compute the depths for Tukey's letter values
#'
#' @param n sample size
get.depth <- function(n){
  res <- c()
  depth <- n
  
  while (depth > 1){
    new.depth <- (floor(depth) + 1) / 2
    res <- c(res, new.depth)
    depth <- new.depth
  }
  
  res <- sort(c(res, n + 1 - res))
  res
}


#' Compute Tukey's letter values
#'
#' Internal function
#' Depends on get.depth()
#'
#' @param smpl sample vector
tukey.quantiles <- function(smpl){
  smpl <- sort(smpl)
  
  n <- length(smpl)  
  depth <- get.depth(n)
  floor.depth <- floor(depth)
  
  # grab letter values
  res <- ifelse(depth==floor.depth, 
                smpl[depth], 
                (smpl[floor.depth] + smpl[floor.depth + 1]) / 2) 
  res
}


#' Summarize quantiles for a given sample
#'
#' Internal function
#' Depends on tukey.quantiles()
#' Depends on get.depth()
#'
#' @param smpl sample vector
quantile.table <- function(smpl){
  n <- length(smpl)
  depth <- get.depth(n)
  N <- length(depth)
  
  lower.depth <- depth[1:(N / 2)]
  lower.depth <- rev(lower.depth)
  upper.depth <- depth[(N / 2 + 1):N]	
  
  quants <- tukey.quantiles(smpl)
  lower.quantiles <- quants[ 1:(N / 2) ]
  lower.quantiles <- rev(lower.quantiles)
  upper.quantiles <- quants[ ((N / 2) + 1):N ]
  
  # compute corresponding z quantiles
  k <- N / 2 - 1
  lower.p <- (1 / 2^(0:k)) / 2 
  upper.p <- 1 - lower.p 
  
  z.lowerQ <- qnorm(lower.p)
  z.upperQ <- qnorm(upper.p)
  
  res <- data.frame(lDepth=lower.depth, uDepth=upper.depth,
  			 		lTail=lower.p, uTail=upper.p, 
  			 		lZQ=z.lowerQ, uZQ=z.upperQ, 
  			 		lLetters=lower.quantiles, uLetters=upper.quantiles)
  
  res  
}


#' Show g for each quantile and estimate g
#'
#' Internal function
#' Depends on quantile.table()
#'
#' @param tabl result from calling quantile.table()
g.table <- function(tabl){
  
  half.area <- tabl$lTail[-1]
  z <- tabl$lZQ[-1]
  
  A <- m <- tabl$lLetters[1]
  LHS <- m - tabl$lLetters[-1]
  UHS <- tabl$uLetters[-1] - m
  z.factor <- -1 * (1 / z)
  g.p <- z.factor * log(UHS / LHS)
  
  g <- median(g.p)
  gMAD <- mad(g.p)
  
  gtab <- data.frame(halfArea=half.area,
  					z = z,
   					zFact=z.factor,
    				LHS=LHS, UHS=UHS, 
    				ratio=UHS/LHS, g.p=g.p)
    				
  res <- list(gtab=gtab, g=g, gMAD=gMAD, A=A)
  res
}

#' Estimate h: given A and g
#'
#' Internal function
#' Depends on g.table()
#'
#' @param tabl result from calling g.table()
h.table <- function(tabl){ 
  g <- tabl$g

  z <- tabl$gtab$z
  halfArea <- tabl$gtab$halfArea
  UHS <- tabl$gtab$UHS
  LHS <- tabl$gtab$LHS
  
  z.sqrd.5 <- z^2 / 2
  
  UHS.star <- g * UHS / (exp(-g * z) - 1)
  LHS.star <- g * LHS / (1 - exp(g * z))
  Spread.star <- g * (UHS + LHS) / (exp(- g * z) - exp(g * z))
  
  # estimate h and B via least-squares
  X <- cbind(1, z.sqrd.5)
  
  #y <- log(UHS.star)
  #y <- log(LHS.star)
  y <- log(Spread.star)
  
  est <- solve(t(X) %*% X) %*% (t(X) %*% y)
  est <- as.numeric(est)
  
  B <- exp(est[1])
  h <- est[2]
    
  htab <- data.frame(halfArea=halfArea, z.sqrd.5=z.sqrd.5,
                     UHS.star=UHS.star, LHS.star=LHS.star, 
                     Spread.star=Spread.star)  
    
  res  <- list(htab=htab, h=h, B=B)
  res  
}