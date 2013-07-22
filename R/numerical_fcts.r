#' Numerical inverse when h > 0
#' Not complete 
#'
#' @export
num.inv.A <- function(){
	if(!h > 0) stop("h must be positive")
	
	if(g==0){
		f <- function(y) trans.gh(y, g=0, h) - x
		# set limits to (-50, 50) for now.
		res <- uniroot(f, c(-50, 50))$root
		res	
	}
 
 	if(g > 0){}
 	
 	if(g < 0){}
 	
 	res

}