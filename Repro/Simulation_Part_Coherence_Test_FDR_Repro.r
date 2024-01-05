
args=(commandArgs(TRUE))

library(HDSpec)

if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  n = 50
  B=1000
  B2=1000
  file='A_Sigma.RData'
  H=20
  d=1
  
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}


library(doParallel)
require(dqrng)
require(doRNG)
cl<-makeCluster(detectCores()-1)
print(detectCores())
registerDoParallel(cl)
registerDoRNG(seed = 1,once=TRUE)


require(abind)

load(file)


q=(p^2-p)/2
q_ind=sapply(1:p,function(i) sapply(1:p,function(j) i>j),simplify = "matrix")
FourierFrequencies=0:(n-1)*2*pi/n

lambda=sqrt(log(p)/n)*exp(seq(from = -8,to = 5,length.out =  100))

M3=16

FourierFrequencies_h=1:(M3-1)*pi/M3

ttilde=2*log(q*M3)


True_Spec_inv=HDSpec::f_spec_inv_ARMA(A = A,B = B_m,Sigma = Sigma,omega_list = FourierFrequencies_h)
True_part_coh=f_inv_to_part_rho(True_Spec_inv)
True_R=apply(abs(True_part_coh),1:2,max)



sim_Res1=foreach(a=icount(B),.combine = cbind,.inorder = FALSE,.packages = c("dqrng","glmnet","abind","foreach","iterators","HDSpectralAnalysis","glassoFast"),.errorhandling = "stop",.export = "q") %dopar% {
  X_t=simulateVARMAs(A = A,B=B_m,Sigma = Sigma,n = n)
  
  
  X_t=t(t(X_t)-apply(X_t,2,mean))
  
  Z=X_to_Z(X_t)
  
  
  M0=adaptive_bandwidth_TLM(X_t,Kern = Kernel,even = T)
  M1=max(4,M0[1])
  N1=floor(M1/(log(M1)^(2/r)))
  
  FourierFrequencies_h1=1:(N1-1)*pi/N1
  
  tmp2=VARLASSO(X_t,VAR_order = VAR_order,lambda = lambda,method = "LASSO",crit="bic_m",mbic2_c =log(p))
  
  M2=max(4,adaptive_bandwidth_TLM(tmp2$eps,Kern = Kernel,even = T)[1])
  N2=floor(M2/(log(M2)^(2/r)))
  
  FourierFrequencies_h2=1:(N2-1)*pi/N2
  
  Spec1=HDSpec::spec_inv_freq(Z = Z,Frequencies = FourierFrequencies_h1,lambda = lambda,Kernel = Kernel,M = M1,M_corr = M0[2])
  Spec1a=Spec1$spec*as.numeric(abs(Spec1$spec)>sapply(Spec1$lambdas,function(i) matrix(i,p,p),simplify = "array"))
  
  Spec_VAR_esta=spec_inv_preWhiteningVAR_freq(X_t = X_t,lambda=lambda,Kernel = Kernel,M = M2,VAR_order = VAR_order,Frequencies = FourierFrequencies_h2,method = "LASSO",crit="bic_m",Check_pos_def = T,M_corr = M0[2],mbic2_c =log(p),gamma=log(log(p)))
  
  erg_t=sapply(c(0,0.1,0.2),function(Delta_0)
  {
    TX_erg=function(TX,Gtx=Gt,FDR_alpha=FDR_alpha,t_hat=NULL,d=NULL,qx=q)
    {
      
      
      if(is.null(t_hat))
      {
        t_hat=t_inf(TX[q_ind],Gt=function(t) Gtx(t,d),q = qx,ttilde = 2*log(qx*(d)),lengt.out = 200,alpha = FDR_alpha)
      }
      return(list(FDR=FDR_TX(TX[q_ind],t_hat,True_R[q_ind],Delta_0),Power=Emp_power(TX[q_ind],t_hat,True_R[q_ind],Delta_0),
                  Test=TX>=t_hat,
                  t_hat=t_hat))
    }
    
    
    Test_TX=max_de_part_coh_H_0(Z=Z,spec = Spec1$spec,Kernel = Kernel,M = M1,DEMEAN = T,Delta_0 = Delta_0,FourierFrequencies_h = FourierFrequencies_h1)
    Test_TX4=max_de_part_coh_H_0_pre_Filt(Z=Z,spec = Spec_VAR_esta$f_inv,Kernel = Kernel,M = M2,DEMEAN = T,Delta_0 = Delta_0,FourierFrequencies_h =  FourierFrequencies_h2,Filter = Array_to_A(Spec_VAR_esta$VAR_A))
    
    test_pen=apply(abs(Spec1a)>Delta_0,1:2,median)
    test_pen2=apply(abs(Spec_VAR_esta$f_inv)>Delta_0,1:2,median)
    test_pen3=apply(abs(Spec1$spec)>Delta_0,1:2,median)
    
    sapply(seq(0.05,0.3,by = 0.05),function(x)
      cbind(TX_erg(Test_TX4$TX_chi*n/M2,Gtx = Gt2t,FDR_alpha = x,d=length(FourierFrequencies_h2)),
            TX_erg(test_pen,Gtx = Gt2t,FDR_alpha = x,t_hat = 0.5),
            TX_erg(test_pen3,Gtx = Gt2t,FDR_alpha = x,t_hat = 0.5),
        TX_erg(Test_TX$TX_chi*n/M1,Gtx = Gt2t,FDR_alpha = x,d=length(FourierFrequencies_h1)),
            TX_erg(test_pen2,Gtx = Gt2t,FDR_alpha = x,t_hat = 0.5)),simplify = "array")
  },simplify = "array"
  )
  
  
  
  tmp=list(FDR=erg_t[1,,,],
           POWER=erg_t[2,,,],
           Test=erg_t[3,,,],
           M1=M1,
           M2=M2)
}


erg_FDR=sapply(sim_Res1[1,],function(x) x,simplify = "array")
erg_POWER=sapply(sim_Res1[2,],function(x) x,simplify = "array")
erg_M2=unlist(sapply(sim_Res1[4,],function(x) x,simplify = "array"))
erg_M2a=unlist(sapply(sim_Res1[5,],function(x) x,simplify = "array"))

erg_POWER_m=apply(erg_POWER,1:3,function(x) mean(unlist(x)))
erg_FDR_m=apply(erg_FDR,1:3,function(x) mean(unlist(x)))

erg_POWER_sd=apply(erg_POWER,1:3,function(x) sd(unlist(x)))
erg_FDR_sd=apply(erg_FDR,1:3,function(x) sd(unlist(x)))


save(A,B_m,Sigma,erg_M2,erg_M2a,True_R,erg_POWER_m,erg_FDR_m,erg_FDR,erg_POWER,erg_FDR_sd,erg_POWER_sd,file=sprintf("Spec_FDR_n=%i_M1=ADAP_VAR_p=%i_d1=%0.2f_d2=%0.2f_%s",n,VAR_order,Delta1,Delta2,file))

#Following operation requires large amount of RAM. if p is large. 
lose_list=function(y) return(apply(y,1:3,function(x) x[[1]]))
erg_TEST=apply(sapply(sim_Res1[3,],function(x) lose_list(x),simplify = "array"),1:4,mean)
erg_TEST_m=array(erg_TEST,c(p,p,dim(erg_FDR_m)[1],dim(erg_FDR_m)[2],dim(erg_FDR_m)[3]))




if(is.null(M1))
{
  save(A,B_m,Sigma,erg_M2,erg_M2a,True_R,erg_TEST_m,erg_POWER_m,erg_FDR_m,erg_FDR,erg_POWER,erg_FDR_sd,erg_POWER_sd,file=sprintf("Spec_FDR_n=%i_M1=ADAP_VAR_p=%i_d1=%0.2f_d2=%0.2f_%s",n,VAR_order,Delta1,Delta2,file))
}else{
  save(A,B_m,Sigma,erg_M2,erg_M2a,True_R,erg_TEST_m,erg_POWER_m,erg_FDR_m,erg_FDR,erg_POWER,erg_FDR_sd,erg_POWER_sd,file=sprintf("Spec_FDR_n=%i_M1=%i_VAR_p=%i_d1=%0.2f_d2=%0.2f_%s",n,round(M1),VAR_order,Delta1,Delta2,file))  
}
