#################################################
#   DATOS DE PESQUERÍA
#################################################
rm(list=ls())
año=1955:1987
cap=c(115.40,118.20,126.40,130.70,146,159.0,148.70,147.6,169.5,162.3,203,195,176.70,143.6,165.1,142.5,202,243.93,157.78,123,90.27,144.55,102.38,102.14,92.9,101.68,101.06,86.75,74.06,83.15,73.84,104.94,99.65)
esf=c(6.67,7.56,7.67,8.04,8.96,9.24,12.3,10.41,12.13,11.12,18.73,18.34,17.65,14.35,19.15,19.71,28.49,49.78,31.75,26.45,19.37,27.02,21.15,17.31,15.15,18.55,17.39,14.78,11.41,12.47,10.13,15.14,15.69)


############MODELO DE SCHAEFER####################################
u=cap/esf
MS=lm(u~esf)
beta=abs(coefficients(MS))
RMS=(beta[1]^2)/(4*beta[2]) 
fRMS=beta[1]/(2*beta[2])

############MODELO DE FOX#########################################
MF=lm(log(u)~esf)
beta1=abs(coefficients(MF))
RMS1=(exp(beta1[1]))/(exp(1)*beta1[2])

#########MODELO GENERALIZADO#########################################
nla=function(n,a){ 
  G=(n^(n/(n-1)))/(n-1)
  MNL=nls(u~b0*(1-(b0/(b1*G))*esf)^(1/(n-1)),start = c(b0= 14.92822,b1=160 )) #Parámetros iniciales obtenidos con el ajuste del método de Schaefer
  par(mfrow=c(1,2))
  b0=coef(MNL)[1]
  b1=coef(MNL)[2]	
  fRMS=((n-1)/n)*G*b1/b0
  plot(esf,cap,ylab = expression(C[e]),xlab ="f",main="Mod. Logístico Generalizado",sub=paste("n=",eval(n)))
  abline(v=fRMS,lty=2,col="purple", lwd=3)
  curve(x*b0*(1-(b0/(G*b1))*x)^(1/(n-1)),add=T,min(esf),max(esf),col=a,sub=paste("n=",eval(n)), lwd=3)
  text(fRMS+2,75, labels = expression(paste(f["RMS"])))
  text(fRMS+2,140, labels = c("RMS"))
  plot(esf,u,ylab = expression(U[e]),xlab ="f", ylim=c(0,20), main = "Ajuste de datos")#,sub=paste("n=",eval(n))
  curve(b0*(1-(b0/(G*b1))*x)^(1/(n-1)),add=T,min(esf),max(esf),col=a,sub=paste("n=",eval(n)), lwd=3)
  print(fRMS)
  return(b1) #Regresa el RMS
}
nla(2, 7)
legend(x=40, y=20, legend=c("n=0.05", "n=0.5", "n=0.75", "n=0.99", "n=1.01", "n=1.5", "n=2"), 
            col=c(4,1,3,2,5,6,7), pch=16, bty="n", y.intersp = .5)

