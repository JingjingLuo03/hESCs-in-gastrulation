library(RColorBrewer)
library(ggplot2)
library(idopNetwork)
library(pracma)
library(stats)
library(reshape2)
list_tf<- c()
for (i in 1:length(toptf_list)) {
  list_tf<- union(list_tf, toptf_list[[i]])
}

tfMatr<-exp_new[netgene,]
times=seq(1,dim(tfMatr)[2],1)
fit_result1<-c()
for(j in 1:dim(tfMatr)[1]){
  data<-tfMatr[j,]
  result<-get_Legendre_par(initial_f_par=rep(0.01,r),data,r,times)
  fit_result1<-cbind(fit_result1,Legendre(result$par,r,times))
}
colnames(fit_result1)<-rownames(tfMatr)
############################functions##############################################
LgdP <- expression( 1,tt,
                    ( 3* tt^2 - 1 )/2 , 
                    ( 5 * tt^3 - 3* tt )/2, 
                    ( 35 * tt^4 - 30 * tt^2 + 3)/8,
                    ( 63 * tt^5 - 70 * tt^3 + 15 * tt )/8,
                    ( 231 * tt^6 - 315 * tt^4 + 105 * tt^2 - 5)/16  )

Legendre<-function(par,r,times){
  tnum = length(times) 
  f<-c()
  for(t in 1:tnum ){
    tt <- -1 + 2*(times[t] - times[1])/(times[tnum] - times[1])
    ff<-0
    for(i in 1:r){
      #eval执行表达式输出结果
      ff<-ff+par[i]*eval(LgdP[i])
    }
    f<-c(f,ff)
  }
  return(f)
}

s.mle<-function(par,data,r,times){
  y <- Legendre(par,r,times)
  yi <- data
  res <- sum((yi-y)^2)
  return(res)
}

get_Legendre_par<-function(initial_f_par,data,r,times){
  #最优化函数,下山单纯形法
  a <- optim(initial_f_par,s.mle,data=data,r=r,times=times,method = "Nelder-Mead",control=list(maxit=2000,trace=FALSE))
  curve_par_i<-a$par
  return(list(par=a$par,res=a$value))
}