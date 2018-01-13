#Packages needed
library(dplyr)
library(BEST)

getMode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#Set expected weights so we can standardize to zero
japWgt <- 170
britWgt <- 150
nzWgt <- 120

#Load Japanese data and British datum subtracting weight so expected mean is 0
japData <- read.csv("japaneseData.csv")$wgt - japWgt
britData <- read.csv("britishData.csv")$wgt - britWgt

priorData <- append(japData, britData)

#Load NZ data
nzData <- read.csv("nzData.csv")$wgt - nzWgt

#Run MCMC for previous data, only set a prior for mu because we expect it to be zero but gave low certainty
priorModel <- BESTmcmc(priorData, prior = list(muM=0, muSD=1000))
plot(priorModel)
summary(priorModel)

nzModel.naive <- BESTmcmc(nzData, prior = list(muM=0, muSD=1000))
plot(nzModel.naive)
summary(nzModel.naive)
summary(nzModel.naive)


nzModel.prior <- BESTmcmc(nzData, prior = list(muM=mean(priorModel$mu), 
                                               muSD=sd(priorModel$mu)*10, 
                                               sigmaMode=summary(nzModel.naive)[2,3], 
                                               sigmaSD=sd(priorModel$sigma*5)))

plot(nzModel.prior)
summary(nzModel.prior)

mu <- seq(-40, 40, 1)
plot(mu, dnorm(mu, mean(priorModel$mu), getMode(priorModel$sigma)), type='l')
lines(mu, dnorm(mu, mean(nzModel.prior$mu), getMode(nzModel.prior$sigma)), type='l')

plot(density(priorModel$mu), xlim=c(-40, 40), lty=2)
lines(density(nzModel.naive$mu), lty=2)
lines(density(nzModel.prior$mu), col="red")
abline(v=0, lty=2)

inverse_gamma_specification <- function(mu, sigma) {
  sigma2 = sigma^2
  mu2 = mu^2
  if(sigma^2 < Inf) {
    nu = sqrt(2*(2+mu^2/sigma^2))
    nu2 = 2*nu
    nu1 = 2
    err = 2*mu^2*gamma(nu/2)^2-(sigma^2+mu^2)*(nu-2)*gamma((nu-1)/2)^2
    while(abs(nu2-nu1) > 1e-12) {
      if(err > 0) {
        nu1 = nu
        if(nu < nu2) {
          nu = nu2
        } else {
          nu = 2*nu
          nu2 = nu
        }
      } else {
        nu2 = nu
      }
      nu =  (nu1+nu2)/2
      err = 2*mu^2*gamma(nu/2)^2-(sigma^2+mu^2)*(nu-2)*gamma((nu-1)/2)^2
    }
    s = (sigma^2+mu^2)*(nu-2)
  } else {
    nu = 2
    s = 2*mu^2/pi
  }
  c(nu=nu, s=s)
}   

inverse_gamma_specification(5,5)
