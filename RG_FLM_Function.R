library(stringr)

############## FLM_Fourier_create ####################


######################################################
FLM_Fourier_create<- function(y,name){

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
  
  y.hat1 = fRegressList$yhatfdobj          # Extract the fitted values
  resid.y1 <- geoslope - y.hat1                 # Residuals
  
  
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
       ylab="Flow Reg. Coeff.",ylim=c(-15e-6,15e-6),
       lwd=2)
  title(main = name)
  title(sub = paste("RSQ = ",round(RSQ1,3),"\n",
                    "Fratio = ",round(Fratio1,3),"\n",
                    "P-value = ",F.res$pval,"\n"),
                    line = -0.5,cex.sub = 0.8)
  #legend("topleft",paste("RSQ = ",round(RSQ1,3),"\n","Fratio = ",round(Fratio1,3)),"\n",
  #        "P-value = ",F.res$pval)
  lines(betafd+2*betastderrfd,lwd=1,lty=2)
  lines(betafd-2*betastderrfd,lwd=1,lty=2)
}

############## FLM_Fourier_create_smoothing ####################


################################################################




FLM_Fourier_create_smoothing<- function(y,name,lambda_input){
  
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
  resid.y1 <- geoslope - y.hat1                 # Residuals
  
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
       ylab="Flow Reg. Coeff.",ylim=c(-15e-6,15e-6),
       lwd=2)
  title(main = name)
  subtitle = paste("RSQ = ",round(RSQ1,3),",",
                    "Fratio = ",round(Fratio1,3),"\n",
                    "P-value = ",F.res$pval,",",
                    "lambda = ",lambda_input)
  mtext(side=1, line=-1, at=10, adj=0, cex=0.5, subtitle)
  #legend("topleft",paste("RSQ = ",round(RSQ1,3),"\n","Fratio = ",round(Fratio1,3)),"\n",
  #        "P-value = ",F.res$pval)
  lines(betafd+2*betastderrfd,lwd=1,lty=2)
  lines(betafd-2*betastderrfd,lwd=1,lty=2)
}

?mtext


FLM_Fourier_create(Fish_logResponse$geoslope,"Geoslope")

par(mfrow=c(3,3))
for(lambda in 10:16){
  FLM_Fourier_create_smoothing(Fish_logResponse$geoslope,"Geoslope",lambda)
}


FLM_Fourier_create(Fish_logResponse$YOY_y,"YOY")


FLM_Fourier_create(Fish_logResponse$recruit_slope,"Recruitment Slope")
FLM_Fourier_create(Fish_logResponse$Oct_Index,"October Index")
