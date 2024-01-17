library(nimble);library(rstan)

head(age_vbgc)

ij<-age_vbgc%>%
    arrange(ID, age_month)%>%
    group_by(ID)%>%
    mutate(j = 1:n())%>%
    arrange(ID)%>%
    group_by(ID)%>%
    mutate(i = cur_group_id())%>%
    ungroup(i, j)%>%
    mutate(J = max(j))%>%
    dplyr::select(ID, i, J, j, age_month, length)

#individuals max num captures (J)    
iJ<-ij%>%
  distinct(iJ)
  
summary(ij)

head(ij)
tail(ij)
#number of individuals
n = max(iJ$i)

model<-nimble::nimbleCode({
  
  for (i in 1:n){ 
    
    for(j in 1:J){
    
    #likelihood
    
    # y = y[i](t[j]) eq. 1 observed length
    # x = g[i](t[j]) = mean distribution specified in terms of the VB growth curve
    
    #Schofield et al. 2013 eq. 1
    y[i,j] ~ T(dnorm(x[i,j], tau.y), 0, Inf)
    #heirarchical
    Linf[i] ~ T(dnorm(beta_l0 + beta_l1, tau.Linf), 0, Inf)
    logit(K[i]) ~ dnorm(beta_k0 + beta_k1, tau.K)      
      #Schofield et al. 2013 eq. 3
    x[i,j] <- Linf[i]*(1 - K[i]^(age[i,j] - t0))
    
    #Midway et al. 2015
    #y[i] ~ dnorm(y.hat[i], tau.y)
    #y.hat[i] <- Linf *(1 - exp(-K *(age[i] - t0)))
    }
  }  
  #prior
  
  beta_l0 ~ dnorm(0, 1000)
  beta_l1 ~ dnorm(0, 1000)
  beta_k0 ~ dlogis(0,1)
  beta_k1 ~ dlogis(0,1)
  
  #Midway et al. 2015
  tau.y <- pow(sigma_y, -2)
  sigma_y ~ dunif(0, 50)
  
  tau.Linf <- pow(sigma_Linf, -2)
  sigma_Linf ~ dunif(0, 50)
  
  tau.K <- pow(sigma_K, -2)
  sigma_K ~ dunif(0, 5)
  
  #remove Linf and K from priors for heirarchical model
  #Linf ~ dnorm(3, 0.4)
  #K ~ dunif(0, 1)
  t0 ~ dnorm(-2.5, 1)
  
})

inits_fun = function(){
  
  #K_init = 0.5
  t0_init = -2.5
  
  
  return(list(t0 = t0_init, beta_10 = 0, beta_l1 = 0, beta_k0 = 0, beta_k1 = 0))
  
}
#Midway
#params = c("Linf", "K", "t0", "sigma_y", "y.hat")
#Schofield
params = c("Linf", "K", "t0", "sigma_y", "x")

data = list(y = as.numeric(ij$length), 
            age = as.numeric(ij$age_month),
            J = as.numeric(iJ$J),
            i = as.numeric(ij$i))


Rmodel <- nimble::nimbleModel(code = model,
                              constants = list(n = n),
                              data = data,
                              inits = inits_fun())

cmodel<-nimble::compileNimble(Rmodel)
mcmc<-buildMCMC(Rmodel, monitors = params)
cmcmc   <- compileNimble(mcmc, project = Rmodel)
samples <- runMCMC(cmcmc, niter = 1000, nburnin = 100, nchains = 3, samplesAsCodaMCMC = T)

coda.samples<-coda::mcmc(samples)

summ<-summary(samples)
summ
summ$statistics[1]
summ$quantiles[c(1,5)]

bayesplot::mcmc_trace(coda.samples)
coda::gelman.diag(coda.samples)

