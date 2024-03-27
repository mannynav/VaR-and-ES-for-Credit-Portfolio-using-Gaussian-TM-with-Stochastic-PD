
#Packages.
library(fitdistrplus)
library(GCPM)
library(QRM)
library(pracma)
library(actuar)

#Size of portfolio under study. LGD data should be in global environment.
PortSize <- length(as.vector(na.omit(as.numeric(lgd$V3))))

#Set correlation.
beta <- 0.20

#Vector to store n-1 summed variables.
vectS<-c(0)

#Vector of homogenous correlations for each loan.
corrVect<- rep(beta, PortSize)

#Number of simulations.
R <- 1500

#Loss given default from data set.
LGD <- suppressWarnings(as.vector(na.omit(as.numeric(lgd$V3))))

#Deterministic Exposure at default.
EAD <- rep(1,PortSize)

#For loop to generate vector of S_n-1 = X_1 + X_2 + ... + X_n-1.
for (i in 1:R) {
  
  #One factor model with F ~ N(0,1).
  F <- rnorm(1)
  
  #Idiosyncratic variables.
  epsilon <- rnorm(PortSize)
  
  #The critical variable X. If X > threshold, no default. If X <= threshold, obligor defaulted.
  X <- sqrt(corrVect)*F + sqrt(1-corrVect)*epsilon
  
  #Stochastic default probability.
  PD <- rprobitnorm(PortSize,qnorm(0.05)/sqrt(1-0.20),0.20/(1-0.20))
  
  #For default only model, the threshold is di = qnorm(PDi).
  threshold <- qnorm(PD)
  
  #Apply threshold to critical variable X.
  defaultedLoans <- X < threshold
  
  #Get obligors that have defaulted.
  defaultedLoans[defaultedLoans == 'TRUE'] <- 1
  
  #Vector for obligors that have defaulted along with their LGD and EAD.
  lossVariable  <- EAD * LGD * defaultedLoans
  
  #Vector to store sum n-1 of Exposure at Default x Loss Given Default for each simulation i.
  vectS[i] <- sum(lossVariable[1:length(lossVariable)-1])
}

#Single sample of loss portfolio that we use fit distributions to.
defaultedLoans[defaultedLoans == 'TRUE'] <- 1
singleSampSPD<-EAD*LGD*defaultedLoans
hist(singleSampSPD)
summary(singleSampSPD)

#Adjust sample to get values above 0.001
AdjustedSample <- singleSampSPD[singleSampSPD>0.001]
hist(AdjustedSample)

#Fit adjusted sample to five distributions.
fit_betaSPD <- fitdist(AdjustedSample,"beta", start = list(shape1 = 0.01, shape2 = 0.01))
fit_paretoSPD <- fitdist(AdjustedSample, distr = "pareto", start = list(shape = 1, scale = 1))
fit_lognormSPD <- fitdist(AdjustedSample, distr = "lnorm")
fit_weibullSPD <- fitdist(AdjustedSample,"weibull",lower = c(0,0))
fit_gammaSPD<-fitdist(AdjustedSample, "gamma")

#Plots of fits.
plot(fit_betaSPD)
plot(fit_paretoSPD)
plot(fit_lognormSPD)
plot(fit_weibullSPD)
plot(fit_gammaSPD)

#Summaries of fits.
summary(fit_betaSPD)
summary(fit_paretoSPD)
summary(fit_lognormSPD)
summary(fit_weibullSPD)
summary(fit_gammaSPD)

summary(vectS)
ScaledTotalLossesSPD <- vectS

#Histogram of distribution of losses.
hist(vectS,xlab= c("Losses"),breaks = 50, main = c("Histogram of Loss Distribution"))

#This is the conditional MC estimate for each S_(n-1) in ScaledTotalLossesSPD.
FncondSPD <- function(x,vectS) {
  s<-c(0) 
  for (i in 1:R) {
    #This distribution function is the mixture distribution function. See MixtureFitExample 4.8.
    s[i] <- pmgp(x-vectS[i], fit_mgp$estimate[1], fit_mgp$estimate[2],fit_mgp$estimate[3], fit_mgp$estimate[4],fit_mgp$estimate[5],fit_mgp$estimate[6])
  }
  return(mean(s))
}

#Quantile.
alpha <- 0.95

#Solve FncondSPD(1/(1-q)-1,ScaledTotalLossesSPD) = alphaSPD for q, store as qu.
qu <- uniroot(f = function(q) alpha - FncondSPD(1/ (1 - q) - 1,vectS), interval = c(0,1))$root

#VaR estimate.
VaR <- 1 /(1 - qu) - 1

#Compute inverse of estimated Fncond.
FnInvsCond <- function(x) {
  sINV<-c(0) 
  for (i in 1:R) {
    sINV[i] <- pmgp(x-vectS[i], fit_mgp$estimate[1], fit_mgp$estimate[2],fit_mgp$estimate[3], fit_mgp$estimate[4],fit_mgp$estimate[5],fit_mgp$estimate[6])
  }
  return(1-mean(sINV))
}

#Vectorize the inverse distribution function.
FnInvsCondVect <- Vectorize(FnInvsCond, "x")

#Integrate the inverse distribution function from the VaRBMM estimate to infinity to get the estimate for m.
m <- integral(FnInvsCondVect,VaR,Inf)

#Compute Expected shortfall.
ES <- VaR + (m)/(1-alpha)
