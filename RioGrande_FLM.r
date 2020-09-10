library("lubridate")
library("ggplot2")
library("dplyr")
library("fda")

######### Hydrology ###################
# Read discharge data, select relevant years and convert to a matrix form 

Albu = read.csv("F:/wild/PhD/ICON/Internship/USGS guage/RioGrande.Albuquerque.1965.2019.csv",header=T)

#slice the dataset for years 2002-2018. Remove 2009 since fish data is not available for this year
Albu = Albu[(Albu$Year>=2002 & Albu$Year <=2018 & Albu$Year != 2009),]

#plot hydrographs yearwise
RioGrande.Albuquerque.2002.2018 = ggplot(Albu)+geom_line(aes(x=Day,y=Q))+facet_wrap(~Year,scales="free_x")+
  ggtitle("RioGrande.Albuquerque.2002.2018")+theme_bw()

length(unique(Albu$Year)) # 16 years of data

#Convert discharge data into a matrix with 365 rows and 16 columns
Albu_flow = matrix(0, nrow=365, ncol= 16)

for(i in 1:16){
  Albu_flow[,i] = Albu[((i*365)-364):(i*365),"Q"]
}

################## Fish Ecology ######################
Fish = read.csv("file:///F:/wild/PhD/ICON/Internship/19_6_2020_new_analysis/Trends_flow_RGSM.06182020.csv",header=T)

# select the response variables: geoslope, YOY_Y, Recruit_slope and October CPUE
names(Fish)
Fish = select(Fish,c("geoslope","YOY_y","recruit_slope","Oct_Index"))
Fish_logResponse = log(Fish+1)
y = Fish_logResponse$geoslope

############################################

######### Functional linear model implementation ##############
# Step 1 : Create a basis function -- fourier or b-spline?
# Step 2 : Convert the hydrograph in to a functional data object. Assess the fit using
# roughness and penalty
# Step 3 : Setup the functional predictor or functional regression coefficient(Beta)
# Step 4 : Select the appropriate response variable and 
# perform regression of the form
# Scalar respone ~ fRegress (time series of explanatory variable)  + optional scalar explanatory variable
# Step 5: Interpret the functional predictor variable
# Step 6: Validate the model using generalized cross validation (GCV), goodness of fit and 
# regression coefficient
# Step 7: Create different e-flow scenarios as functional data objects and check ecological response
################ 


########### Step 1 : Create a basis function -- fourier or b-spline? #############
# According to Ramsay et al. fourier series is used for explanatory variables which are periodic
# and b-splines are preferred for non-periodic functions
# Stewart-Koster et al. use b-spline for periodic functions. So, trying both for now

# fourier basis function with the number of basis functions as 65.
# Chosing 65 as nbasis to get started with.
# Also, the canadian weather example from Ramsay et al., uses 65 basis function for precip data 
daybasis65 = create.fourier.basis(c(0,365),65)
plot(daybasis65)

########### Step 2 : Create a functional data object using basis function and explanatory variable  #############
# using fourier basis

# check for year 2002 to see the fit
#tempSmooth_f=smooth.basis(dayOfYear,log(Albu_flow[1:365]+1),daybasis65) #dayOfYear = c(1,365)
#plot(tempSmooth_f)
#summary(tempSmooth_f)
# add actual discharge values on the plot of functional data object
#plot(log(Albu_flow[1:365]+1)

# do it for all the 16 years
tempSmooth_f=smooth.basis(day.5,Albu_flow,daybasis65)
plot(tempSmooth_f)
#save the functional data object
tempfd_f =tempSmooth_f$fd


#store the functional data covariates in a list
#fourier
templist_f = vector("list",2)
templist_f[[1]] = rep(1,16)
templist_f[[2]] = tempfd_f


################## Step 3: Create a functional predictor or functional regression coefficient(Beta) #########
# from Ramsay et al. creating a regression coefficient (beta) with 5 fourier basis function
conbasis = create.constant.basis(c(0,365))
betabasis = create.fourier.basis(c(0,365),5)
betalist = vector("list",2)
betalist[[1]] = conbasis
betalist[[2]] = betabasis


######## Step 4 : Select the appropriate response variable and perform regression of the form ############
# Scalar respone ~ fRegress (time series of explanatory variable)  + optional scalar explanatory variable

# Using log(geoslope) as the response variable to begin with
fRegressList1 = fRegress(y,templist_f,betalist)


############# Step 5: Interpret the functional predictor variable #################
## The regression coefficients are stored in fRegressList$betaestlist
betaestlist1 = fRegressList1$betaestlist
tempbetafd1 = betaestlist1[[2]]$fd
plot(tempbetafd1, xlab="Day",
     ylab="Beta for geoslope")

y.hat1 = fRegressList1$yhatfdobj          # Extract the fitted values
resid.y1 <- y - y.hat1                 # Residuals

# permutation based approach to calculate p-value. See Ramsay et al. for details
F.res = Fperm.fd(y, templist_f,betalist, plotres=F)

# plot y vs. fitted y
plot(y~y.hat1,xlab="Fitted values",ylab="Observed geoslope",main="Geoslope fit")
abline(0,1)
plot(resid.y1 ~ y,xlab = "Fitted values",ylab = "Residuals",main = "Geoslope: Residual plot")

#calcualate squared multiple correlation and F ratio 
SSE1.1 = sum(resid.y1^2)
SSE0 = sum((y - mean(y))^2)
RSQ1 = (SSE0-SSE1.1)/SSE0
Fratio1 = ((SSE0-SSE1.1)/(fRegressList1$df-1))/(SSE1.1/(length(y)-fRegressList1$df))

#Confidence intervals based on Ramsay et. al 2009
resid = y - fRegressList1$yhatfdobj
SigmaE. = sum(resid^2)/(length(y)-fRegressList1$df)
SigmaE = SigmaE. *diag(rep(1,length(y)))
y2cMap = tempSmooth_f$y2cMap
stderrList = fRegress.stderr(fRegressList1, y2cMap,SigmaE)

betafdPar = fRegressList1$betaestlist[[2]]
betafd = betafdPar$fd
betastderrList = stderrList$betastderrlist
betastderrfd = betastderrList[[2]]

maxfd = betafd+2*betastderrfd
minfd = betafd-2*betastderrfd

#plot regression coefficient
plot.fd(betafd, xlab="Day", ylab="Flow Reg. Coeff.",cex.axis=3,
     ylim = c(-15e-6,15e-6),
     lwd=2,cex.lab=5, cex.main=1, cex.sub=1.5)

lines(betafd+2*betastderrfd,lwd=1,lty=2)
lines(betafd-2*betastderrfd,lwd=1,lty=2)

title(main = name)
title(sub = paste("RSQ = ",round(RSQ1,3),"\n",
                  "Fratio = ",round(Fratio1,3),"\n",
                  "P-value = ",F.res$pval,"\n"),
      line = -0.5,cex.sub = 0.8)

subtitle = paste("RSQ = ",round(RSQ1,3),"\n",
                 "Fratio = ",round(Fratio1,3),"\n",
                 "P-value = ",F.res$pval,"\n")


############Using a smoothing parameter########
Lcoef = c(0,(2*pi/365)^2,0)
plot(Lcoef)

# Harmonic accelerator
harmaccelLfd = vec2Lfd(Lcoef, c(0,365))

accelLfd <- int2Lfd(2)
# defining the beta estimate by a functional parameter object that
# incorporates both this roughness penalty and a level of
# smoothing
betabasis = create.fourier.basis(c(0, 365), 65)
lambda = 10^12
betafdPar = fdPar(betabasis, accelLfd, lambda)
betalist[[2]] = betafdPar

# regression using lambda and harmonic accelerator
fRegressList2 = fRegress(y, templist_f,betalist)
betaestlist2 = fRegressList2$betaestlist
tempbetafd2 = betaestlist2[[2]]$fd
plot(tempbetafd2, xlab="Day",
     ylab="Beta for geoslope")

y.hat2 = fRegressList2$yhatfdobj          # Extract the fitted values
resid.y2 <- y - y.hat2               # Residuals

# plot y vs. fitted y
plot(y~y.hat1)
points(y~y.hat2, col = "Red")
plot(resid.y1 ~ y)
points(resid.y2~y, col = "red")

#calcualate squared multiple correlation and F ratio 
SSE1.2 = sum(resid.y2^2)
SSE0 = sum((y - mean(y))^2)

RSQ2 = (SSE0-SSE1.2)/SSE0
Fratio2 = ((SSE0-SSE1.2)/(fRegressList2$df-1))/(SSE1.1/(length(y)-fRegressList2$df))

#chosing lambda using CV approach
loglam = seq(10,15,0.5)
nlam = length(loglam)
SSE.CV = matrix(0,nlam,1)
SSE.GCV = matrix(0,nlam,1)
SSE.OCV = matrix(0,nlam,1)
for (ilam in 1:nlam) {
  lambda = 10^loglam[ilam]
  betalisti = betalist
  betafdPar2 = betalisti[[2]]
  betafdPar2$lambda = lambda
  betalisti[[2]] = betafdPar2
  fRegi = fRegress.CV(y, templist_f,betalisti)
  fRegi2 = fRegress(y, templist_f,betalisti)
  SSE.CV[ilam] = fRegi$SSE.CV
  SSE.GCV[ilam] = fRegi2$gcv
  SSE.OCV[ilam] = fRegi2$OCV
}

plot(SSE.CV~loglam, xlab = "Log(lambda)", ylab = "Crossvalidation score", 
     main = "Smoothing and penalty estimation for Geoslope")

