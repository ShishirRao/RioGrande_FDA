library(stringr)
library(ggplot2)
library(ggplot2)


######################################################
FLM_Fourier_create<- function(y,name,ylimit){

  daybasis65 = create.fourier.basis(c(0,365),65)
  
  tempSmooth_f=smooth.basis(day.5,Albu_flow,daybasis65)
  tempfd_f =tempSmooth_f$fd
  
  
  #store the functional data covariates in a list
  #fourier
  templist_f = vector("list",2)
  templist_f[[1]] = rep(1,16)
  templist_f[[2]] = tempfd_f
  
  
  conbasis = create.constant.basis(c(0,365))
  betabasis = create.fourier.basis(c(0,365),5)
  betalist = vector("list",2)
  betalist[[1]] = conbasis
  betalist[[2]] = betabasis
  
  
  fRegressList = fRegress(y,templist_f,betalist)
  
  
  betaestlist = fRegressList$betaestlist
  tempbetafd = betaestlist[[2]]$fd
  #plot(tempbetafd, xlab="Day",
  #     ylab="Beta for geoslope")
  
  y.hat1 = fRegressList$yhatfdobj       # Extract the fitted values
  resid.y1 <- y - y.hat1                # Residuals
  
  
  #calcualate squared multiple correlation and F ratio 
  SSE1.1 = sum(resid.y1^2)
  SSE0 = sum((y - mean(y))^2)
  
  RSQ1 = (SSE0-SSE1.1)/SSE0
  Fratio1 = ((SSE0-SSE1.1)/(fRegressList$df-1))/(SSE1.1/(length(y)-fRegressList$df))
  
  F.res = Fperm.fd(y, templist_f,betalist, plotres=F)
  
  #Confidence intervals
  resid = y - fRegressList$yhatfdobj
  SigmaE. = sum(resid^2)/(length(y)-fRegressList$df)
  SigmaE = SigmaE. *diag(rep(1,length(y)))
  y2cMap = tempSmooth_f$y2cMap
  stderrList = fRegress.stderr(fRegressList, y2cMap,SigmaE)
  
  
  betafdPar = fRegressList$betaestlist[[2]]
  betafd = betafdPar$fd
  betastderrList = stderrList$betastderrlist
  betastderrfd = betastderrList[[2]]
  
  maxfd = betafd+2*betastderrfd
  minfd = betafd-2*betastderrfd
  
  plot(betafd, xlab="Day",
       ylab="Flow Reg. Coeff.",ylim=c(-ylimit,ylimit),
       lwd=2,cex.lab=1.5, cex.axis=3, cex.main=1.5, cex.sub=1.5)
  
  title(main = name,cex.lab=1.5, cex.axis=3, cex.main=1.5, cex.sub=1.5)
  title(sub = paste("RSQ = ",round(RSQ1,3),",",
                    "Fratio = ",round(Fratio1,3),"\n",
                    "P-value = ",round(F.res$pval,3),",",
                    "df =", round(fRegressList$df,1)),
                    line = -0.9,cex.lab=1.5, cex.axis=3, cex.main=1.5, cex.sub=1.5)
  #legend("topleft",paste("RSQ = ",round(RSQ1,3),"\n","Fratio = ",round(Fratio1,3)),"\n",
  #        "P-value = ",F.res$pval)
  lines(betafd+2*betastderrfd,lwd=1,lty=2)
  lines(betafd-2*betastderrfd,lwd=1,lty=2)
  
}


#########################FLM_Fourier_create_smoothing############
FLM_Fourier_create_smoothing<- function(y,name,lambda_input,ylimit){
  
  daybasis65 = create.fourier.basis(c(0,365),65)
  
  tempSmooth_f=smooth.basis(day.5,Albu_flow,daybasis65)
  tempfd_f =tempSmooth_f$fd
  
  
  #store the functional data covariates in a list
  #fourier
  templist_f = vector("list",2)
  templist_f[[1]] = rep(1,16)
  templist_f[[2]] = tempfd_f
  
  Lcoef = c(0,(2*pi/365)^2,0)
  
  # Harmonic accelerator
  harmaccelLfd = vec2Lfd(Lcoef, c(0,365))
  
  conbasis = create.constant.basis(c(0,365))
  betalist = vector("list",2)
  betalist[[1]] = conbasis
  
  #use a saturated beta
  betabasis = create.fourier.basis(c(0, 365), 16)
  lambda = 10^lambda_input
  betafdPar = fdPar(betabasis, accelLfd, lambda)
  betalist[[2]] = betafdPar
  
  fRegressList = fRegress(y,templist_f,betalist)
  
  betaestlist = fRegressList$betaestlist
  tempbetafd = betaestlist[[2]]$fd
  
  y.hat1 = fRegressList$yhatfdobj          # Extract the fitted values
  resid.y1 <- y - y.hat1                 # Residuals
  
  #calcualate squared multiple correlation and F ratio 
  SSE1.1 = sum(resid.y1^2)
  SSE0 = sum((y - mean(y))^2)
  
  RSQ1 = (SSE0-SSE1.1)/SSE0
  Fratio1 = ((SSE0-SSE1.1)/(fRegressList$df-1))/(SSE1.1/(length(y)-fRegressList$df))
  
  F.res = Fperm.fd(y, templist_f,betalist, plotres=F)
  
  #Confidence intervals
  resid = y - fRegressList$yhatfdobj
  SigmaE. = sum(resid^2)/(length(y)-fRegressList$df)
  SigmaE = SigmaE. *diag(rep(1,length(y)))
  y2cMap = tempSmooth_f$y2cMap
  stderrList = fRegress.stderr(fRegressList, y2cMap,SigmaE)
  
  
  betafdPar = fRegressList$betaestlist[[2]]
  betafd = betafdPar$fd
  betastderrList = stderrList$betastderrlist
  betastderrfd = betastderrList[[2]]
  
  maxfd = betafd+2*betastderrfd
  minfd = betafd-2*betastderrfd
  
  plot(betafd, xlab="Day",
       ylab="Flow Reg. Coeff.",ylim=c(-ylimit,ylimit),
       lwd=2)
  title(main = name)
  subtitle = paste("RSQ = ",round(RSQ1,3),",",
                    "Fratio = ",round(Fratio1,3),",",
                    "P-value = ",round(F.res$pval,3),"\n",
                    "df =", round(fRegressList$df,1),",",
                    "CV =",round(fRegressList$OCV,3))
  
  mtext(side=1, line=-1, at=10, adj=0, cex=0.6, subtitle)
  #legend("topleft",paste("RSQ = ",round(RSQ1,3),"\n","Fratio = ",round(Fratio1,3)),"\n",
  #        "P-value = ",F.res$pval)
  lines(betafd+2*betastderrfd,lwd=1,lty=2)
  lines(betafd-2*betastderrfd,lwd=1,lty=2)
}

###################################################
par(mfrow=c(2,2))

#print FLM output with reduced dimension beta
FLM_Fourier_create(log(Fish$geoslope+1),"Log(Geoslope+1)", 15e-6)
FLM_Fourier_create(log(Fish$YOY_y+1),"Log(YOY+1)",70e-6)
FLM_Fourier_create(log(Fish$recruit_slope+1),"Log(Recruitment Slope+1)",1e-6)
FLM_Fourier_create(log(Fish$Oct_Index+1),"Log(Oct Index+1)",50e-6)

par(mfrow=c(3,3))
for(lambda in 10:16){
  FLM_Fourier_create_smoothing(log(Fish$geoslope+1),
                               paste("log(Geoslope+1) with lambda = ",lambda),
                               lambda,15e-6)
}
