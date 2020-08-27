library(fda)


########### b-spline example ############
(mean_age = mean(growth$age))
?create.bspline.basis

hgtbasis = create.bspline.basis(c(1,18), 35, 6, growth$age)
hgtbasis = create.bspline.basis(c(1,18), 12, 6, growth$age)
plot(hgtbasis)
?Data2fd
hgtffd = Data2fd(growth$hgtf,growth$age,hgtbasis)

plotfit.fd(growth$hgtf,age,hgtffd)


############ fourier basis ###########
?create.fourier.basis()

# Create a minimal Fourier basis for annual data
#  using 3 basis functions
yearbasis3 <- create.fourier.basis(c(0,365),
                                   axes=list("axesIntervals") )
#  plot the basis
plot(yearbasis3)


# fourier basis with period = 365
daybasis65 = create.fourier.basis(c(0,365),65)
plot(daybasis65)
summary(daybasis65)


CanadianWeather
daily

###### setting up a functional regression model ##########

# response variable is a scalar -- annual rainfall from 35 Canadian weather stations
weatherData = daily
annualprec = log10(apply(daily$precav,2,sum))
length(annualprec)

# creating the basis function 
tempbasis =create.fourier.basis(c(0,365),65)
plot(tempbasis)

# creating the functional response variable
tempSmooth=smooth.basis(day.5,daily$tempav,tempbasis)
plot(tempSmooth)
tempfd =tempSmooth$fd

#capture the fd coefficients in a list for regresssion
templist = vector("list",2)
templist[[1]] = rep(1,35)
templist[[2]] = tempfd

#define the type of beta
conbasis = create.constant.basis(c(0,365))
betabasis = create.fourier.basis(c(0,365),5)
betalist = vector("list",2)
betalist[[1]] = conbasis
betalist[[2]] = betabasis

#perform regression
fRegressList = fRegress(annualprec,templist,betalist)

betaestlist = fRegressList$betaestlist
tempbetafd = betaestlist[[2]]$fd
plot(tempbetafd, xlab="Day",
     ylab="Beta for temperature")

#plot y vs. fitted y
plot(annualprec~fRegressList$yhatfdobj)
# look at residuals vs y
residuals = (annualprec-fRegressList$yhatfdobj)
plot(residuals~annualprec) 

#calcualate squared multiple correlation and F ratio 
annualprechat1 = fRegressList$yhatfdobj
annualprecres1 = annualprec - annualprechat1
SSE1.1 = sum(annualprecres1^2)
SSE0 = sum((annualprec - mean(annualprec))^2)

RSQ1 = (SSE0-SSE1.1)/SSE0
Fratio1 = ((SSE0-SSE1.1)/(fRegressList$df-1))/(SSE1.1/(length(annualprec)-fRegressList$df))


############Using a smoothing parameter########
Lcoef = c(0,(2*pi/365)^2,0)
plot(Lcoef)

# Harmonic accelerator
harmaccelLfd = vec2Lfd(Lcoef, c(0,365))
plot(harmaccelLfd[[1]])

# defining the beta estimate by a functional parameter object that
# incorporates both this roughness penalty and a level of
# smoothing
betabasis = create.fourier.basis(c(0, 365), 35)
lambda = 10^12.5
betafdPar = fdPar(betabasis, harmaccelLfd, lambda)
betalist[[2]] = betafdPar

# regression using lambda and harmonic accelerator
annPrecTemp = fRegress(annualprec, templist,
                       betalist)
betaestlist2 = annPrecTemp$betaestlist
annualprechat2 = annPrecTemp$yhatfdobj

#Is the fit better than without lambda?
plot(annualprec~fRegressList$yhatfdobj)
points(annualprec~annPrecTemp$yhatfdobj,lwd=2)

residuals2 = (annualprec-annPrecTemp$yhatfdobj)
plot(residuals~annualprec) 
points(residuals2~annualprec,lwd =2)

#calcualate squared multiple correlation and F ratio 
SSE1.2 = sum((annualprec-annualprechat2)^2)
(RSQ2 = (SSE0 - SSE1.2)/SSE0)
(Fratio2 = ((SSE0-SSE1.2)/(annPrecTemp$df-1))/(SSE1.2/(length(annualprec)-annPrecTemp$df)))

print(c(RSQ1,Fratio1))
print(c(RSQ2,Fratio2))

#Confidence intervals
resid = annualprec - annPrecTemp$yhatfdobj
SigmaE.= sum(resid^2)/(35-annPrecTemp$df)
SigmaE = SigmaE. *diag(rep(1,35))
y2cMap = tempSmooth$y2cMap
stderrList = fRegress.stderr(annPrecTemp, y2cMap,
                             SigmaE)


betafdPar = annPrecTemp$betaestlist[[2]]
betafd = betafdPar$fd
betastderrList = stderrList$betastderrlist
betastderrfd = betastderrList[[2]]
plot(betafd, xlab="Day",
     ylab="Temperature Reg. Coeff.",
     ylim=c(-6e-4,1.2e-03), lwd=2)
lines(betafd+2*betastderrfd, lty=2, lwd=1)
lines(betafd-2*betastderrfd, lty=2, lwd=1)


###### are these models any better than a constant value of Beta ##########

betalist[[2]] = fdPar(conbasis)
fRegressList = fRegress(annualprec, templist,
                        betalist)
betaestlist = fRegressList$betaestlist

annualprechat = fRegressList$yhatfdobj
SSE1 = sum((annualprec-annualprechat)^2)
RSQ = (SSE0 - SSE1)/SSE0
Fratio = ((SSE0-SSE1)/1)/(SSE1/33)

print(c(RSQ1,Fratio1))
print(c(RSQ2,Fratio2))
print(c(RSQ,Fratio))


#chosing a smoothing parameter
# initialise smoothing parameters
betabasis = create.fourier.basis(c(0, 365), 35)
lambda = 10^12.5
betafdPar = fdPar(betabasis, harmaccelLfd, lambda)
betalist[[2]] = betafdPar


loglam = seq(5,15,0.5)
nlam = length(loglam)
SSE.CV = matrix(0,nlam,1)
for (ilam in 1:nlam) {
  lambda = 10^loglam[ilam]
  betalisti = betalist
  betafdPar2 = betalisti[[2]]
  betafdPar2$lambda = lambda
  betalisti[[2]] = betafdPar2
  fRegi = fRegress.CV(annualprec, templist,
                      betalisti)
  SSE.CV[ilam] = fRegi$SSE.CV
}

plot(SSE.CV~loglam)


######### understanding the functional object creation #########
daytime = (1:365)-0.5

plot(daily, col =1, lty=1)
lines(daily$tempav[,2])
lines(daily$tempav[3,])


JJindex = c(182:365, 1:181)
tempmat = daily$tempav[JJindex,]

tempbasis = create.fourier.basis(c(0,365),65)

tempfd = smooth.basis(daytime, tempmat, tempbasis)$fd

tempfd$fdnames = list("Day (July 2 to June 30)",
                      "Weather Station",
                      "Mean temperature (deg. C)")

plot(tempfd, col=1, lty=1)


####### sine wave example ##############
seq = seq(0,10,1)
length(seq)
(nbasis = length(seq)+4-2)
basis13 = create.bspline.basis(c(0,10),nbasis,4,seq)

tvec = seq(0,1,len=13)
sinecoef = sin(2*pi*tvec)
class(sinecoef)
plot(sinecoef)
sinefd = smooth.basis(c(0,13),sinecoef,basis13)

sinefd = fd(sinecoef, basis13, list("t","","f(t)"))
op = par(cex=1.2)
plot(sinefd, lwd=2)
points(tvec*10, sinecoef, lwd=2)
par(op)
?fd
?s


############# example code #################
##
######## Simulated data example 1: a simple regression smooth  ########
##
#  Warning:  In this and all simulated data examples, your results
#  probably won't be the same as we saw when we ran the example because
#  random numbers depend on the seed value in effect at the time of the
#  analysis.
#
#  Set up 51 observation points equally spaced between 0 and 1
n = 51
argvals = seq(0,1,len=n)
#  The true curve values are sine function values with period 1/2
(x = sin(4*pi*argvals))

plot(x)
plot(x~argvals)
?rnorm
#  Add independent Gaussian errors with std. dev. 0.2 to the true values
sigerr = 0.2
y = x + rnorm(x)*sigerr
lines(y~argvals)

#  When we ran this code, we got these values of y (rounded to two
#  decimals):
y = c(0.27,  0.05,  0.58,  0.91,  1.07,  0.98,  0.54,  0.94,  1.13,  0.64,
      0.64,  0.60,  0.24,  0.15, -0.20, -0.63, -0.40, -1.22, -1.11, -0.76,
      -1.11, -0.69, -0.54, -0.50, -0.35, -0.15,  0.27,  0.35,  0.65,  0.75,
      0.75,  0.91,  1.04,  1.04,  1.04,  0.46,  0.30, -0.01, -0.19, -0.42,
      -0.63, -0.78, -1.01, -1.08, -0.91, -0.92, -0.72, -0.84, -0.38, -0.23,
      0.02)
#  Set up a B-spline basis system of order 4 (piecewise cubic) and with
#  knots at 0, 0.1, ..., 0.9 and 1.0, and plot the basis functions
length(y)

length(seq(0,1,0.1))

nbasis = 13 # 11+4-2
basisobj = create.bspline.basis(c(0,1),nbasis)
plot(basisobj)
#  Smooth the data, outputting only the functional data object for the
#  fitted curve.  Note that in this simple case we can supply the basis
#  object as the "fdParobj" parameter
ys = smooth.basis(argvals=argvals, y=y, fdParobj=basisobj)
Ys = smooth.basis(argvals=argvals, y=y, fdParobj=basisobj,
                  returnMatrix=TRUE)
# Ys[[7]] = Ys$y2cMap is sparse;  everything else is the same
plot(ys$fd)


all.equal(ys[-7], Ys[-7])


(xfd = ys$fd)
Xfd = Ys$fd

#  Plot the curve along with the data
plotfit.fd(y, argvals, xfd)
#  Compute the root-mean-squared-error (RMSE) of the fit relative to the
#  truth
RMSE = sqrt(mean((eval.fd(argvals, xfd) - x)^2))
print(RMSE)  #  We obtained 0.069
#  RMSE = 0.069 seems good relative to the standard error of 0.2.
#  Range through numbers of basis functions from 4 to 12 to see if we
#  can do better.  We want the best RMSE, but we also want the smallest
#  number of basis functions, which in this case is the degrees of
#  freedom for error (df).  Small df implies a stable estimate.
#  Note: 4 basis functions is as small as we can use without changing the
#  order of the spline.  Also display the gcv statistic to see what it
#  likes.
for (nbasis in 4:12) {
  basisobj = create.bspline.basis(c(0,1),nbasis)
  ys = smooth.basis(argvals, y, basisobj)
  xfd = ys$fd
  gcv = ys$gcv
  RMSE = sqrt(mean((eval.fd(argvals, xfd) - x)^2))
  # progress report:
   cat(paste(nbasis,round(RMSE,3),round(gcv,3),"\n"))
}
#  We got RMSE = 0.062 for 10 basis functions as optimal, but gcv liked
#  almost the same thing, namely 9 basis functions.  Both RMSE and gcv
#  agreed emphatically that 7 or fewer basis functions was not enough.
#  Unlike RMSE, however, gcv does not depend on knowing the truth.
#  Plot the result for 10 basis functions along with "*" at the true
#  values
nbasis = 10
basisobj = create.bspline.basis(c(0,1),10)
xfd = smooth.basis(argvals, y, basisobj)$fd
plotfit.fd(y, argvals, xfd)
points(argvals,x, pch="*")
#  Homework:
#  Repeat all this with various values of sigerr and various values of n


?smooth.basis

