data<-read.table("rh200100100064.dat")
B=10000
n=length(data)
boot=apply(matrix(sample(data,size=n*B,replace=T),nrow=B,ncol=n),1,mean)
quantile(boot, c(0.025, 0.975))

