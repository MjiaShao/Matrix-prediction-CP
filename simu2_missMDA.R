
library(missMDA)
library(Matrix)
source("./data_generate.R")
source("./data_generate_error.R")

args <- commandArgs(trailingOnly = TRUE)
set.seed(200)
typemissing = "maxmissing"
alpha = 0.1
totalboot = 300
scaledata = 5
Nall = c(50,100,200)
graphonname_all <-c("f1","f2", "f3")
misperall <- seq(0.05,0.20,0.05)
for(n in Nall){
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
  
  # Get the row and column index of the selected entry
  row_index <- n
  col_index <- n-1
  pred <- data[row_index,col_index]
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

  timestart = Sys.time()
  imputed_mat <- MIPCA(data, ncp = 2,nboot = 100)
  everyimpute = imputed_mat$res.MI
  nboot =100
  predall = c()
  for (i in 1:nboot){
    predma = everyimpute[[i]]
    predall[i] = predma[row_index,col_index]
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
df = data.frame(cover = cover, CIlen = CIlen, timerecord = timerecord)
write.csv(df, sprintf("./result/conf_net_missMDA_%d_%d_%s_con_%d_rand_%s_%d.csv",n,floor(fixed_u*10),graphonname,floor(misper*100),typemissing, floor(scaledata)))

}}}