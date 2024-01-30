## this is where I worked some things out before putting it in the measurements_demo.Rmd file
# L. Crowe Jan 2024

library(nimble);library(rstan)

age_vbgc<-age_calc%>%
  filter(!is.na(length_use))%>%
  dplyr::rename("length" = "length_use")

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
#as.matrix(age_ij)

length_ij<-ij%>%
  dplyr::select(ID, j, length)%>%
  group_by(ID)%>%
  tidyr::pivot_wider(names_from = j, values_from = length)%>%
  ungroup()%>%
  dplyr::select(-ID)

colnames(length_ij)<-NULL
as.matrix(length_ij)[1,7]

max_ind_j<-ij%>%
  dplyr::select(ID, j)%>%
  group_by(ID)%>%
  filter(j == max(j))%>%
  ungroup()%>%
  dplyr::select(-ID)

max_ind_j
#number of individuals
#n = max(ij$i)
n = nrow(age_ij)
n
nrow(length_ij)

#capture occasion
#J = max(ij$j)

head(ij)

model<-nimble::nimbleCode({
  
  for (i in 1:n){ 
    
    for(j in 1:max_ind[i]){
      
      #likelihood
      
      #notes ---
      # y = y[i](t[j]) eq. 1 observed length
      # x = g[i](t[j]) = predicted length/growth curve
      # ---
      
      #Schofield et al. 2013 eq. 1
      y[i,j] ~ dnorm(x[i,j], tau_y) #constrain to positive
      
      #Schofield heirarchical
      Linf[i] ~ dnorm(beta_Linf, tau_Linf) #constrain to positive
      logit(K[i]) ~ dlogis(beta_K, tau_K) # K = e^-k, constrain to 0, 1
      log(t0p) ~ dlnorm(beta_t0p, tau_t0p) # throw a negative on here to get t0
      
      #Schofield et al. 2013 eq. 3
      x[i,j] <- Linf[i]*(1 - K[i]^(age[i,j] + t0p))
      
      #x[i,j] <- Linf[i] *(1 - exp(-K[i] *(age[i,j] + t0)))
      
    }
  }
  #prior
  
  beta_Linf ~ dnorm(0, 1000)
  beta_K ~ dlogis(0, 1)
  beta_t0p ~ dlnorm(0, 5)
  
  tau_y <- pow(sigma_y, -2)
  sigma_y ~ dnorm(0, 500)
  
  tau_Linf <- pow(sigma_Linf, -2)
  sigma_Linf ~ dnorm(0, 5)
  
  tau_K <- pow(sigma_K, -2)
  sigma_K ~ dnorm(0, 5)
  
  tau_t0p <- pow(sigma_t0p, -2)
  sigma_t0p ~ dnorm(0, 5)
  
  #remove Linf and K from priors for heirarchical model
  #Linf ~ dnorm(3, 0.4)
  #K ~ dunif(0, 1)
  #t0 ~ dnorm(-2.5, 1)
  
  })

inits_fun = function(){
  
  #K_init = 0.5
  #t0_init = -2.5
  
  
  return(list(beta_Linf = 0, beta_K = 0, beta_t0p = 0))
  
}


params = c("Linf", "K", "t0p", "sigma_y", "beta_Linf", "beta_K", "beta_t0p")
params = c("t0p", "beta_Linf", "beta_K", "beta_t0p")

data = list(y = as.matrix(length_ij), 
            age = as.matrix(age_ij))


Rmodel <- nimble::nimbleModel(code = model,
                              constants = list(n = n, max_ind = max_ind_j$j),
                              data = data,
                              inits = inits_fun())
#Rmodel$initializeInfo

cmodel<-nimble::compileNimble(Rmodel)
mcmc<-buildMCMC(Rmodel, monitors = params)
cmcmc   <- compileNimble(mcmc, project = Rmodel)
samples <- runMCMC(cmcmc, niter = 10000, nburnin = 1000, nchains = 3, samplesAsCodaMCMC = T)

coda.samples<-coda::mcmc(samples)

summ<-summary(samples)
summ
samples.df<-as.data.frame(as.matrix(samples))
summ$statistics[805]

bayesplot::mcmc_trace(coda.samples)
coda::gelman.diag(coda.samples)

