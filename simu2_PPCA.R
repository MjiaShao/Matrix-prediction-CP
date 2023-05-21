source("./Function_main_algebraic.R")
source("./Functions_supp_algebraic.R")
source("./OtherFunctions.R")
source("./data_generate.R")
source("./data_generate_error.R")
library(FactoMineR)
library(missMDA)


set.seed(200)
r =2
sigma <- 0.1
alpha = 0.1
totalboot = 300
typemissing = "maxmissing"
scaledata = 5
n = 50
graphonname_all <-c("f1","f2", "f3")
misperall <- seq(0.05,0.20,0.05)
for(graphonname in graphonname_all){
	if (graphonname == "f1"){
		fixed_u = 0.9
	} else if (graphonname == "f2"){
		fixed_u = 0.7
	} else if (graphonname == "f3"){
		fixed_u = 0.6
	}
for(misper in misperall){
cover = rep(0,totalboot)
CIlen = rep(0,totalboot)
timerecord = rep(0,totalboot)
misn = ceiling(n*misper)

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
  Y = data
  if(typemissing == 'minmissing'){
     idx <- order(data[upper.tri(data)], decreasing = FALSE)[misn]
     data[data<=data[upper.tri(data)][idx]] <- NA
  } else if (typemissing == 'maxmissing'){
     idx <- order(data[upper.tri(data)], decreasing = TRUE)[misn]
     data[data>=data[upper.tri(data)][idx]] <- NA
  }
  diag(data) = 0  
  data[row_index,col_index] = NA
  data[col_index,row_index] = NA
  YNA <- data
  timestart = Sys.time()
   YNAimp <- imputePCA(YNA,ncp=r)
  YNAimp <- YNAimp$completeObs
  M = 1 - is.na(YNA)
  indMissVar = rowMeans(M)!=1
  indMissVar = seq(n)[indMissVar]
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
  print(k)
}
df = data.frame(cover = cover, CIlen = CIlen, timerecord = timerecord)
write.csv(df, sprintf("./result/conf_net_ppca_%d_%d_%s_%d_con_new_rand_%s_%d.csv",n,floor(fixed_u*10),graphonname,floor(misper*100),typemissing, floor(scaledata)))
print(misper)
}
  print(graphonname)
}
