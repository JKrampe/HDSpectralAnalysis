library(doParallel)
require(dqrng)
require(doRNG)
cl<-makeCluster(detectCores()-1)
print(detectCores())
registerDoParallel(cl)
registerDoRNG(seed = 1,once=TRUE)

require(abind)
abind2=function(...) abind(...,along=3)

require(network)
require(R.matlab)

result=NULL
resultb=NULL

Numbers=c(19,20)
State=c("Open","Closed")



for(ab in State)
  for(ba in Numbers)
  {
    
    Frequencyband=c(1/25*2*pi,1/14*2*pi)
  
    tmp=R.matlab::readMat(sprintf('Data_Subset/EEG_Cat_Study4_Resting_Data_S%i_Eyes%s_Gaussianized.mat',ba,ab))

    if(ab=="Open")
    {
      Data_array=tmp$DATA2
    }
    if(ab=='Closed')
    {
      Data_array=tmp$DATA1
    }
      
    
    
    Channel_loc=read.csv("Data_Subset/Channel_locations.csv")
    
    erg1a=foreach(a=icount(dim(Data_array)[1]),.combine = cbind,.inorder = FALSE,.packages = c("dqrng","glmnet","abind","foreach","iterators","HDSpectralAnalysis","glassoFast"),.errorhandling = "stop") %dopar% {
      
      X_t=t(Data_array[,,a])
      p=dim(X_t)[2]
      n=dim(X_t)[1]
      PC_Network(X_t,lambda=sqrt(log(p)/n)*exp(seq(from = -5,to = 4,length.out =  50)),Frequencyband =  Frequencyband,FDR_alpha = 0.1,Delta_0 = 0)  
    }
    
    save(erg1a,file=sprintf("Data_Subset/EEG_Cat_Study4_Resting_Eyes%s_FinalEpochs_S%i_Result.RData",ab,ba))
    erg2_mean=apply(sapply(erg1a[2,],function(x) x,simplify = "array"),1:2,mean)
    erg3_mean=apply(sapply(erg1a[3,],function(x) x,simplify = "array"),1:2,mean)
    
    
    
    result=abind(result,erg2_mean,along=3)
    resultb=abind(resultb,erg3_mean,along=3)
    
    tikzDevice::tikz(file = sprintf("PC_Network_%s_%i.tex",ab,ba),width = 16,height = 9)
    plot(network(erg2_mean>0.5,directed = F,loops = F),coord=Channel_loc[,2:3],label = Channel_loc$V5,vertex.cex = 1,arrowhead.cex = 1,jitter=FALSE)
    dev.off()
    
  }


result

tikzDevice::tikz(file = sprintf("PC_Network_Pat_19.tex",ab,ba),width = 16,height = 9)
par(mfrow=c(1,2))
plot(network(result[,,1]>0.5,directed = F,loops = F),coord=Channel_loc[,2:3],label = Channel_loc$V5,vertex.cex = 1,arrowhead.cex = 1,jitter=FALSE,main="Eyes Open")
plot(network(result[,,3]>0.5,directed = F,loops = F),coord=Channel_loc[,2:3],label = Channel_loc$V5,vertex.cex = 1,arrowhead.cex = 1,jitter=FALSE,main="Eyes Closed")
dev.off()

tikzDevice::tikz(file = sprintf("PC_Network_Pat_20.tex",ab,ba),width = 16,height = 9)
par(mfrow=c(1,2))
plot(network(result[,,2]>0.5,directed = F,loops = F),coord=Channel_loc[,2:3],label = Channel_loc$V5,vertex.cex = 1,arrowhead.cex = 1,jitter=FALSE,main="Eyes Open")
plot(network(result[,,4]>0.5,directed = F,loops = F),coord=Channel_loc[,2:3],label = Channel_loc$V5,vertex.cex = 1,arrowhead.cex = 1,jitter=FALSE,main="Eyes Closed")
dev.off()



tikzDevice::tikz(file = sprintf("PC_Network_Pat_Coh_19.tex",ab,ba),width = 16,height = 9)
par(mfrow=c(1,2))
plot(network(resultb[,,1]>0.5,directed = F,loops = F),coord=Channel_loc[,2:3],label = Channel_loc$V5,vertex.cex = 1,arrowhead.cex = 1,jitter=FALSE,main="Coherence")
plot(network(result[,,1]>0.5,directed = F,loops = F),coord=Channel_loc[,2:3],label = Channel_loc$V5,vertex.cex = 1,arrowhead.cex = 1,jitter=FALSE,main="Partial Coherence")
dev.off()



pdf(file = sprintf("PC_Network_Pat_19.pdf"),width = 16,height = 9)
par(mfrow=c(1,2))
plot(network(result[,,1]>0.5,directed = F,loops = F),coord=Channel_loc[,2:3],label = Channel_loc$V5,vertex.cex = 1,arrowhead.cex = 1,jitter=FALSE,main="Eyes Open")
plot(network(result[,,3]>0.5,directed = F,loops = F),coord=Channel_loc[,2:3],label = Channel_loc$V5,vertex.cex = 1,arrowhead.cex = 1,jitter=FALSE,main="Eyes Closed")
dev.off()

pdf(file = sprintf("PC_Network_Pat_20.pdf"),width = 16,height = 9)
par(mfrow=c(1,2))
plot(network(result[,,2]>0.5,directed = F,loops = F),coord=Channel_loc[,2:3],label = Channel_loc$V5,vertex.cex = 1,arrowhead.cex = 1,jitter=FALSE,main="Eyes Open")
plot(network(result[,,4]>0.5,directed = F,loops = F),coord=Channel_loc[,2:3],label = Channel_loc$V5,vertex.cex = 1,arrowhead.cex = 1,jitter=FALSE,main="Eyes Closed")
dev.off()



pdf(file = sprintf("PC_Network_Pat_Coh_19.pdf"),width = 16,height = 9)
par(mfrow=c(1,2))
plot(network(resultb[,,1]>0.5,directed = F,loops = F),coord=Channel_loc[,2:3],label = Channel_loc$V5,vertex.cex = 1,arrowhead.cex = 1,jitter=FALSE,main="Coherence")
plot(network(result[,,1]>0.5,directed = F,loops = F),coord=Channel_loc[,2:3],label = Channel_loc$V5,vertex.cex = 1,arrowhead.cex = 1,jitter=FALSE,main="Partial Coherence")
dev.off()

p=dim(Data_array)[1]
q=(p^2-p)/2
q_ind=sapply(1:p,function(i) sapply(1:p,function(j) i>j),simplify = "matrix")

write.csv(data.frame(Student_state=as.vector(outer(Numbers,State,FUN = paste)),Coherence_connectivy=round(apply(resultb>0.5,3,function(x) mean(x[q_ind]))*100,1),partial_coherence_connectivity=round(apply(result>0.5,3,function(x) mean(x[q_ind]))*100,1)),file = "Connectivity.csv")

