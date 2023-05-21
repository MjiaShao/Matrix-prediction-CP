setwd("C:/phd4/network yuan/conformal_pred/github_CP")
source("./data_generate.R")
source("./data_generate_error.R")
library(softImpute)

set.seed(1)
alpha = 0.1
totalboot = 1000
scaledata = 5
Nall = c(50,100,200,400)
xi_n <- seq(0.1,0.9,0.1)
graphonname_all <-c("f1","f2", "f3")
for(n in Nall){
for(fixed_u in xi_n){
for(graphonname in graphonname_all){
cover = rep(0,totalboot)
CIlen = rep(0,totalboot)
timerecord = rep(0,totalboot)
for(k in 1:totalboot){
  x <- runif(n)
  x[n] <- fixed_u
  data <- generate_data(x, graphonname,fixed_u)*scaledata
  data <- data_add_error(data, 0.1)
  data[lower.tri(data)] <- t(data)[lower.tri(data)]
  diag(data) <- 0
  row_index <- n
  col_index <- n-1
  pred <- data[row_index,col_index]
  data[row_index,col_index] = NA
  data[col_index,row_index] = NA
  timestart = Sys.time()
  imputed_SVD <- softImpute(data,rank=2,lambda=1)
  imputed_matnew <- data
  imputed_mat <- imputed_SVD$u %*% diag(imputed_SVD$d)  %*% t(imputed_SVD$v)
  errors <- data[!is.na(data)] - imputed_mat[!is.na(data)]
  nboot =100
  predall = c()
  for (i in 1:nboot){
    errorsnew <- sample(errors, size = length(errors), replace = FALSE)
    imputed_matnew[!is.na(data)] <- imputed_mat[!is.na(data)] + errorsnew
    imputed_SVDnew <- softImpute(imputed_matnew,rank=2,lambda=1)
    completenew <- complete(imputed_matnew, imputed_SVDnew)
    predall[i] <- completenew[row_index,col_index]+sample(errors,1)
  }
  lower_bound <- quantile(predall, alpha/2)
  upper_bound <- quantile(predall, 1-alpha/2)
  if (lower_bound <pred &upper_bound >pred){
    cover[k] = 1
  }
  CIlen[k] = upper_bound- lower_bound
  timeend = Sys.time()
  timerecord[k] = as.numeric(difftime(timeend, timestart, units = "secs"))
}
df = data.frame(cover = cover, CIlen = CIlen, timerecord=timerecord)
write.csv(df, sprintf("./result/conf_net_soft_%d_%d_%s_con_%d.csv",n,floor(fixed_u*10),graphonname,floor(scaledata)))
}}}

