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
plot(log(Fish$geoslope+1)~Fish$year)
Fish = select(Fish,c("geoslope","YOY_y","recruit_slope","Oct_Index"))
Fish_logResponse = log(Fish+1)
geoslope = Fish_logResponse$geoslope

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
tempSmooth_f=smooth.basis(dayOfYear,Albu_flow[1:365],daybasis65) #dayOfYear = c(1,365)
plot(tempSmooth_f)
summary(tempSmooth_f)
# add actual discharge values on the plot of functional data object
points(Albu_flow[1:365])

# do it for all the 16 years
tempSmooth_f=smooth.basis(dayOfYear,Albu_flow,daybasis65)
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

# not sure how to do the b-spline version of beta yet. TBD

######## Step 4 : Select the appropriate response variable and perform regression of the form ############
# Scalar respone ~ fRegress (time series of explanatory variable)  + optional scalar explanatory variable

# Using log(geoslope) as the response variable to begin with
fRegressList = fRegress(geoslope,templist_f,betalist)


############# Step 5: Interpret the functional predictor variable #################
## The regression coefficients are stored in fRegressList$betaestlist

betaestlist = fRegressList$betaestlist
view(betaestlist)
geoslopebetafd = betaestlist[[2]]$fd
plot(geoslopebetafd, xlab="Day",
     ylab="Beta for geoslope")

y.hat = fRegressList$yhatfdobj          # Extract the fitted values
resid.y <- geoslope - y.hat                 # Residuals

plot(geoslope~y.hat)

# Develop approximate pointwise confidence intervals
sigmae. <- sum(resid.y^2)/(length(geoslope)-fRegressList$df)   #    Using 2* SE of the reg coefficient at each time t
sigmae <- sigmae.*diag(rep(1,length(geoslope)))
y2cMap <- tempSmooth_f$y2cMap
stderrList <- fRegress.stderr(fRegressList, y2cMap, sigmae)

# Extract the functional  reg coef for plotting
betafdpar <- betaestlist[[2]]
betafd <- betafdpar$fd
betastderrList <- stderrList$betastderrlist
betastderrfd <- betastderrList[[2]]


plot(geoslopebetafd, xlab="Day", ylab = "Flow reg coeff")
lines(betafd+2*betastderrfd)
lines(betafd-2*betastderrfd)


#Assess the quality of this fit
geoslopehat1 = fRegressList$yhatfdobj
geosloperes1 = Fish_logResponse$geoslope - geoslopehat1
(SSE1.1 = sum(geosloperes1^2))
(SSE0 = sum((Fish_logResponse$geoslope - mean(Fish_logResponse$geoslope))^2))


(RSQ1 = (SSE0-SSE1.1)/SSE0)
Fratio1 = ((SSE0-SSE1)/5)/(SSE1/29)



fRegressList$df




