library(nimble);library(rstan)

head(age_vbgc)

ij<-age_vbgc%>%
    arrange(ID, age_month)%>%
  ungroup()%>%
    mutate(ind = 1:n())%>%
    group_by(ID)%>%
    mutate(j = 1:n(),
           J = max(ind),
           i = cur_group_id() 
           )%>%
    #arrange(ID)%>%
    #group_by(ID)%>%
    #mutate(i = cur_group_id())%>%
    #mutate(J = max(j))%>%
    dplyr::select(ID, i, ind, j, J, age_month, length)%>%
  ungroup()

head(ij)
tail(ij)
summary(ij)

age_ij<-ij%>%
  dplyr::select(ID, j, age_month)%>%
  group_by(ID)%>%
  tidyr::pivot_wider(names_from = j, values_from = age_month)%>%
  ungroup()%>%
  dplyr::select(-ID)

colnames(age_ij)<-NULL
age_ij
  
length_ij<-ij%>%
  dplyr::select(ID, j, length)%>%
  group_by(ID)%>%
  tidyr::pivot_wider(names_from = j, values_from = length)%>%
  ungroup()%>%
  dplyr::select(-ID)

colnames(length_ij)<-NULL
as.matrix(length_ij)

#number of individuals
#n = max(ij$i)
n = nrow(age_ij)
n
nrow(length_ij)

J = max(ij$j)
J

model<-nimble::nimbleCode({
  
  for (i in 1:n){ 

    for(j in 1:J){
    
    #likelihood
    
    # y = y[i](t[j]) eq. 1 observed length
    # x = g[i](t[j]) = mean distribution specified in terms of the VB growth curve
    
     #Schofield et al. 2013 eq. 1
    y[i,j] ~ T(dnorm(x[i,j], tau.y), 0, Inf)
    #y[i] ~ T(dnorm(x[i], tau.y), 0, Inf)
     #Schofield heirarchical
    Linf[i] ~ T(dnorm(beta_l, tau.Linf), 0, Inf)
    logit(K[i]) ~ dnorm(beta_k, tau.K)
    log(t0) ~ dnorm(beta_t0, tau.t0)
      #Schofield et al. 2013 eq. 3
    x[i,j] <- Linf[i]*(1 - K[i]^(age[i,j] - t0))
    #x[i] <- Linf[i]*(1 - K[i]^(age[i] - t0))
    
    #Midway et al. 2015
    #y[i] ~ dnorm(y.hat[i], tau.y)
    #y.hat[i] <- Linf *(1 - exp(-K *(age[i] - t0)))
    }
  }
  #prior
  
  beta_l ~ dnorm(0, 1000)
  #beta_l1 ~ dnorm(0, 1000)
  beta_k ~ dlogis(0, 1)
  #beta_k1 ~ dlogis(0,1)
  beta_t0 ~ dnorm(0, 5)
  
  #Midway et al. 2015
  tau.y <- pow(sigma_y, -2)
  sigma_y ~ dunif(0, 50)
  
  tau.Linf <- pow(sigma_Linf, -2)
  sigma_Linf ~ dunif(0, 50)
  
  tau.K <- pow(sigma_K, -2)
  sigma_K ~ dunif(0, 5)
  
  tau.t0 <- pow(sigma_t0, -2)
  sigma_t0 ~ dunif(0, 5)
  
  #remove Linf and K from priors for heirarchical model
  #Linf ~ dnorm(3, 0.4)
  #K ~ dunif(0, 1)
  #t0 ~ dnorm(-2.5, 1)
  
  })

inits_fun = function(){
  
  #K_init = 0.5
  #t0_init = -2.5
  
  
  return(list(beta_1 = 0, beta_k = 0, beta_t0 = 0))
  
}
#Midway
#params = c("Linf", "K", "t0", "sigma_y", "y.hat")
#Schofield
params = c("Linf", "K", "t0", "sigma_y")

list(y = matrix(seizures$seize, ncol = J, nrow = n, byrow = TRUE))

data = list(y = as.matrix(length_ij), 
            age = as.matrix(age_ij))


Rmodel <- nimble::nimbleModel(code = model,
                              constants = list(n = n, J = J),
                              data = data,
                              inits = inits_fun())

cmodel<-nimble::compileNimble(Rmodel)
mcmc<-buildMCMC(Rmodel, monitors = params)
cmcmc   <- compileNimble(mcmc, project = Rmodel)
samples <- runMCMC(cmcmc, niter = 1000, nburnin = 100, nchains = 3, samplesAsCodaMCMC = T)

coda.samples<-coda::mcmc(samples)

summ<-summary(samples)
summ
samples.df<-as.data.frame(as.matrix(samples))
summ$statistics[805]
summ$quantiles[c(1,5)]

bayesplot::mcmc_trace(coda.samples)
coda::gelman.diag(coda.samples)

