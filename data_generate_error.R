data_add_error = function(data, errorscale){
    n = dim(data)[1]
    error = matrix(runif(n*n, -errorscale, errorscale),nrow = n)
    return(data+error)
}