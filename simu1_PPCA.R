
setwd("C:/phd4/network yuan/conformal_pred/github_CP")

source("./Function_main_algebraic.R")
source("./Functions_supp_algebraic.R")
source("./OtherFunctions.R")
source("./data_generate.R")
source("./data_generate_error.R")
library(FactoMineR)
library(softImpute)
library(missMDA)

set.seed(200)
alpha = 0.1
totalboot = 300
r =2
sigma <- 0.1
scaledata = 5

n = 50;
xi_n <- seq(0.1,0.9,0.1)
graphonname_all <-c("f1","f2", "f3")
for(fixed_u in xi_n){
for(graphonname in graphonname_all){
indMissVar <- (n-1):n
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
  
  # Get the row and column index of the selected entry
  row_index <- n
  col_index <- n-1
  pred <- data[row_index,col_index]
  Y = data
  data[row_index,col_index] = NA
  data[col_index,row_index] = NA
  YNA <- data
  obs <- c() #indexes of the complete rows
  for (i in 1:nrow(YNA)){
      if(sum(is.na(YNA[i,]))==0){
        obs <- c(obs,i)
      }
    }
  timestart = Sys.time()
  YNAimp <- imputePCA(YNA,ncp=r) 
  YNAimp <- YNAimp$completeObs
  M = 1 - is.na(YNA)
  indRegVarAll <- setdiff(1:n,indMissVar)
  indRegVar <- NULL
  if (is.null(indRegVar)){
    indRegVar=setdiff(1:n,indMissVar)
  }
  res_estim_MNAR <- Mean_covariances_estimations_algebraic(YNA,indMissVar,indRegVar,r,opt="agg", YNAimp = YNAimp, opt_data = "real")
  MeanMNAR <- res_estim_MNAR$mean
  CovMNAR <- res_estim_MNAR$cov
  B <- matrix(rnorm(n * r), nrow = r, ncol = n)
  res_PPCA_MNAR <- Results_PPCA_imputation(CovMNAR,MeanMNAR,YNA,Y,B,indMissVar,M,r,sigma)
  pointpred <-  res_PPCA_MNAR$YNA[n,n-1]
  lower_bound <- pointpred - qnorm(1-alpha/2)*abs(CovMNAR[n-1,n-1])
  upper_bound <- pointpred + qnorm(1-alpha/2)*abs(CovMNAR[n-1,n-1])
  
  if (lower_bound <pred &upper_bound >pred){
    cover[k] = 1
  }
  CIlen[k] = upper_bound- lower_bound
  timeend = Sys.time()
  timerecord[k] = as.numeric(difftime(timeend, timestart, units = "secs"))
}
df = data.frame(cover = cover, CIlen = CIlen, timerecord = timerecord)
write.csv(df, sprintf("./result/conf_net_ppca_%d_%d_%s_con_%d.csv",n,floor(fixed_u*10),graphonname,floor(scaledata)))
}}

