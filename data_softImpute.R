library(R.matlab)
library(softImpute)
dataori = readMat('./data/Final_data.mat')
alpha = 0.1
typemissing <-   "maxmissing"
typeall <- c("NC","SZ")
nummissenall <- c(300, 1000)
seedall <- seq(8)
for(type in typeall){
for(nummissen in nummissenall){
for(seed in seedall){
if (type == "SZ"){
  dataall = dataori$cor.SZ.fisherZ.correct
} else {
  dataall = dataori$cor.NC.fisherZ.correct
}
totalboot = dim(dataall)[3]
n = dim(dataall)[1]
covermean = c()
CIlenmean = c()
timemean = c()
coverall = rep(0, n*totalboot)
dim(coverall) = c(n, totalboot)
CIlenall  = rep(0, n*totalboot)
dim(CIlenall) =  c(n, totalboot)
timeall =rep(0, n*totalboot)
dim(timeall) = c(n, totalboot)
set.seed(seed)
for(rowinx in 1:n){
  cover = rep(0,totalboot)
  CIlen = rep(0,totalboot)
  timerecord = rep(0,totalboot)
  colindices = seq(n)
  colindices = colindices[-rowinx]
  for(k in 1:totalboot){
    data <- as.matrix(dataall[,,k])
    row_index <- rowinx
    col_index <- sample(colindices,1)
    pred <- data[row_index,col_index]
    if(typemissing == 'minmissing'){
      idx <- order(data[upper.tri(data)], decreasing = FALSE)[nummissen]
      data[data<=data[upper.tri(data)][idx]] <- NA
    } else if (typemissing == 'maxmissing'){
      idx <- order(data[upper.tri(data)], decreasing = TRUE)[nummissen]
      data[data>=data[upper.tri(data)][idx]] <- NA
    }
    diag(data) = 0
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
      completenew <- softImpute::complete(imputed_matnew, imputed_SVDnew)
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
  covermean[rowinx] = mean(cover)
  CIlenmean[rowinx] = mean(CIlen)
  timemean[rowinx] = mean(timerecord)
  coverall[rowinx,] = cover
  CIlenall[rowinx,] = CIlen
  timeall[rowinx,] = timerecord
}    

df = data.frame(covermean = covermean, CIlenmean = CIlenmean, timemean = timemean)
write.csv(df, sprintf("./result/data_softImpute_%s_random_mislarge_%d_%s_%d.csv",type,floor(nummissen),typemissing,floor(seed)))
write.table(coverall, sprintf("./result/data_softImpute_%s_random_mislarge_%d_%s_%d_cover.txt",type,floor(nummissen),typemissing,floor(seed)), append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE)
write.table(CIlenall, sprintf("./result/data_softImpute_%s_random_mislarge_%d_%s_%d_CIlen.txt",type,floor(nummissen),typemissing,floor(seed)), append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE)
write.table(timeall, sprintf("./result/data_softImpute_%s_random_mislarge_%d_%s_%d_timecost.txt",type,floor(nummissen),typemissing,floor(seed)), append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE)

}}}



