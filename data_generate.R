generate_data = function(x, graphonname,fixed_u){
  n <- length(x)
  XX <- x %*% t(rep(1, n))
  YY <- rep(1, n) %*% t(x)
  switch(graphonname,
         "f1" = (XX + YY) / 2 - 0.15,
	     "f2" = 0.5*(cos(0.1/((XX-0.5)^3+(YY-0.5)^3+0.01)))*(pmax(XX^(2/3),YY^(2/3))) + 0.4,
         "f3" = ((XX^2 + YY^2) / 3 * cos(1 / (XX^4 + YY^4)) + 0.15),
         stop("Unknown distribution")
  )
}
