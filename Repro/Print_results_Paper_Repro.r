# Print Table -------------------------------------------------------------

require(abind)

files=c("DGP_VARMA(1,1)_p50.RData","DGP_VARMA(1,1)_p100.RData"
        ,"DGP_VARMA(1,1)_p200.RData"
        ,"DGP_VARMA(0,5)_p50.RData","DGP_VARMA(0,5)_p100.RData"
        ,"DGP_VARMA(0,5)_p200.RData"
)


ns=c(190,2^9,2^11,2^12)


FDRs=NULL
Powers=NULL
FDRs_sd=NULL
Powers_sd=NULL
names=NULL
for(file in files)
  for(n in ns)
  {
    load(sprintf("Spec_FDR_n=%i_M1=%s_VAR_p=%i_d1=%0.2f_d2=%0.2f_%s",n,"ADAP",round(log(n,base = 10)),0,0.2,file))
    names=c(names,sprintf("%s_%i",substr(file,5,nchar(file)-6),n))
    if(is.null(FDRs))
    {
      FDRs=erg_FDR_m
      Powers=erg_POWER_m
      FDRs_sd=erg_FDR_sd
      Powers_sd=erg_FDR_sd
    }else{
      FDRs=abind(FDRs,erg_FDR_m,along = 4)
      Powers=abind(Powers,erg_POWER_m,along = 4)
      FDRs_sd=abind(FDRs_sd,erg_FDR_sd,along = 4)
      Powers_sd=abind(Powers_sd,erg_POWER_sd,along = 4)  
    }
  }

FDRs=round(FDRs,2)
Powers=round(Powers,2)
FDRs_sd=round(FDRs_sd,2)
Powers_sd=round(Powers_sd,2)

mean_sd=function(m,sd)
{
  sapply(1:length(m),function(i) return(sprintf("%0.2f(%0.2f)",m[i],sd[i])))
}


require(xtable)

tmp=cbind(mean_sd(FDRs[1,1,1,],FDRs_sd[1,1,1,]),mean_sd(Powers[1,1,1,],Powers_sd[1,1,1,]),
          mean_sd(FDRs[1,2,1,],FDRs_sd[1,2,1,]),mean_sd(Powers[1,2,1,],Powers_sd[1,2,1,]),
          mean_sd(FDRs[2,1,1,],FDRs_sd[2,1,1,]),mean_sd(Powers[2,1,1,],Powers_sd[2,1,1,]))

dimnames(tmp)=list(names,c('FDR_0.05','Power_0.05','FDR_0.1','Power_0.1','Reg_FDR','Reg_Power'))
print.xtable(xtable(tmp),
             hline.after = (1:6)*4,file="Result_table_out.tex")


tmp=cbind(mean_sd(FDRs[1,1,3,],FDRs_sd[1,1,3,]),mean_sd(Powers[1,1,3,],Powers_sd[1,1,3,]),
          mean_sd(FDRs[1,2,3,],FDRs_sd[1,2,3,]),mean_sd(Powers[1,2,3,],Powers_sd[1,2,3,]),
          mean_sd(FDRs[2,1,3,],FDRs_sd[2,1,3,]),mean_sd(Powers[2,1,3,],Powers_sd[2,1,3,]))

dimnames(tmp)=list(names,c('FDR_0.05','Power_0.05','FDR_0.1','Power_0.1','Reg_FDR','Reg_Power'))
print.xtable(xtable(tmp),
             hline.after = (1:6)*4,file="Result_table_out_Delta_0.2.tex")






# Figures

# False Positve/False Negative --------------------------------------------


require(HDSpectralAnalysis)
files=c("DGP_VARMA(1,1)_p50.RData","DGP_VARMA(0,5)_p50.RData")

ns=c(2^9,2^12)


FDRs=NULL
Powers=NULL
names=NULL
Delta1=0
for(file in files)
  
{
  load(sprintf("Spec_FDR_n=%i_M1=%s_VAR_p=%i_d1=%0.2f_d2=%0.2f_%s",ns[1],"ADAP",ceiling(log(ns[1],base = 10)),0,0.2,file))
  p=dim(True_R)[1]
  q=(p^2-p)/2
  q_ind=sapply(1:p,function(i) sapply(1:p,function(j) i>j),simplify = "matrix")
  name_proc=substr(file,5,nchar(file)-6)
  names=c(names,sprintf("%i_%s",p,substr(file,5,nchar(file)-6)))
  
  
  
  erg_TEST_m2=erg_TEST_m
  erg_POWER_m2=erg_POWER_m
  erg_FDR_m2=erg_FDR_m
  load(sprintf("Spec_FDR_n=%i_M1=%s_VAR_p=%i_d1=%0.2f_d2=%0.2f_%s",ns[2],"ADAP",ceiling(log(ns[2],base=10)),0,0.2,file))
  
  
  tmp1=(1-erg_TEST_m[,,1,1,1])
  tmp1[q_ind]=(1-erg_TEST_m2[,,1,1,1])[q_ind]
  
  tmp2=(1-erg_TEST_m[,,2,1,1])
  tmp2[q_ind]=(1-erg_TEST_m2[,,2,1,1])[q_ind]
  
  tikzDevice::tikz(sprintf("n=512_4096_%s_d=%0.2f.tex",name_proc,Delta1),width = 16/100*p,height = 9/100*p,onefile = T)
  par(mfrow=c(1,2))
  # par(bg="black")
  plotmat(((tmp1)*(True_R>Delta1 & True_R<1)-(1-tmp1)*(True_R<=Delta1)),range_A = c(-1.001,1.001),color1 = "blue","red",main=sprintf("Testing"))
  lines(p:1,p:1,col="black")
  plotmat(((tmp2)*(True_R>Delta1 & True_R<1)-(1-tmp2)*(True_R<=Delta1)),range_A = c(-1.001,1.001),color1 = "blue","red",main="Regularizing")
  lines(p:1,p:1,col="black")
  dev.off()
}


# FDR/Power Curve Lambda --------------------------------------------------
load(file="Spec_FDR_Lambda_Seq_n=512_M1=ADAP_VAR_p=3_d1=0.00_d2=0.20_DGP_VARMA(1,1)_p50.RData")
p=50
n=512
lambda=sqrt(log(p)/n)*exp(seq(from = -8,to = 5,length.out =  100))



tikzDevice::tikz("FDR_Power_lambda_005.tex",width = 16/2,height = 7/2,onefile = T)
i=1
plot(y = erg_FDR_m[1,i,1,],x = lambda,ylim=c(0,1),type="l",log="x",ylab="FDR(solid) Power(dashed)",xlim=range(sqrt(log(p)/n)*exp(seq(from = -4,to = 3,length.out =  100))),xlab="$\\lambda$")
abline(h=0.05*i,lty=4)
lines(y = erg_FDR_m[2,i,1,],x = lambda,col="red")
lines(y = erg_POWER_m[1,i,1,],x = lambda,lty=2,lwd=2)
lines(y = erg_POWER_m[2,i,1,],x = lambda,col="red",lty=2,lwd=2)
legend("topright",legend = c("Testing","Regularizing"),col=c("black","red"),lty=c(1,1))
dev.off()

tikzDevice::tikz("FDR_Power_lambda_020.tex",width = 16/2,height = 7/2,onefile = T)
i=4
plot(y = erg_FDR_m[1,i,1,],x = lambda,ylim=c(0,1),type="l",log="x",ylab="FDR(solid) Power(dashed)",xlim=range(sqrt(log(p)/n)*exp(seq(from = -4,to = 3,length.out =  100))),xlab="$\\lambda$")
abline(h=0.05*i,lty=4)
lines(y = erg_FDR_m[2,i,1,],x = lambda,col="red")
lines(y = erg_POWER_m[1,i,1,],x = lambda,lty=2,lwd=2)
lines(y = erg_POWER_m[2,i,1,],x = lambda,col="red",lty=2,lwd=2)
legend("topright",legend = c("Testing","Regularizing"),col=c("black","red"),lty=c(1,1))
dev.off()
