## called in measurements_demo.Rmd file
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
as.matrix(age_ij)[52,6]
saveRDS(age_ij, file = "./data/age_ij.rds")

length_ij<-ij%>%
  dplyr::select(ID, j, length)%>%
  group_by(ID)%>%
  tidyr::pivot_wider(names_from = j, values_from = length)%>%
  ungroup()%>%
  dplyr::select(-ID)

colnames(length_ij)<-NULL
as.matrix(length_ij)[52,6]
saveRDS(age_ij, file = "./data/length_ij.rds")

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
      #Schofield et al. 2013 eq. 1
      y[i,j] ~ dnorm(x[i,j], tau = tau_y) #constrain to positive
      
      #Schofield heirarchical
      Linf[i] ~ dnorm(beta_Linf, tau = tau_Linf) #constrain to positive
      logit(K[i]) ~ dnorm(beta_K, tau = tau_K) # K = e^-k, constrain to 0, 1

      #Schofield et al. 2013 eq. 3
      x[i,j] <- Linf[i]*(1 - (K[i]^(age[i,j] + t0p)))
      
    }
  }
  #prior
  
  beta_Linf ~ dnorm(0, pow(1000, -2))
  beta_K ~ dlogis(0, 1)
  #beta_t0p ~ dlnorm(0, 1)
  
  tau_y <- pow(sigma_y, -2)
  sigma_y ~ T(dt(0, pow(50, -2), df= 3), 0, Inf)
  
  tau_Linf <- pow(sigma_Linf, -2)
  sigma_Linf ~ T(dt(0, pow(50, -2), df= 3), 0, Inf)
  
  tau_K <- pow(sigma_K, -2)
  sigma_K ~ T(dt(0, pow(5, -2), df= 3), 0, Inf)
  
  log(t0p) ~ dnorm(0, 0.1) # throw a negative on here to get t0 proper

  })

inits_fun = function(){
  
  return(list(beta_Linf = 0, beta_K = 0))
  
}


params = c("Linf", "K", "t0p", "sigma_y", "beta_Linf", "sigma_Linf", "beta_K", "sigma_K")

data = list(y = as.matrix(length_ij), 
            age = as.matrix(age_ij))

##### one option for running model ---- 

# Rmodel <- nimble::nimbleModel(code = model,
#                               constants = list(n = n, max_ind = max_ind_j$j),
#                               data = data,
#                               inits = inits_fun())
# #Rmodel$initializeInfo
# 
# cmodel<-nimble::compileNimble(Rmodel)
# mcmc<-buildMCMC(Rmodel, monitors = params)
# cmcmc   <- compileNimble(mcmc, project = Rmodel)
# samples <- runMCMC(cmcmc, niter = 10000, nburnin = 1000, nchains = 3, samplesAsCodaMCMC = T)

#### second option ----

samples<-nimble::nimbleMCMC(code = model,
                            constants = list(n = n, max_ind = max_ind_j$j),
                            data = data,
                            inits = inits_fun(),
                            niter = 10000,
                            nburnin = 1000, 
                            nchains = 3,
                            monitors = params,
                            samplesAsCodaMCMC = T
)

#### -----
coda.samples<-coda::mcmc(samples)

coda::gelman.diag(coda.samples)
bayesplot::mcmc_trace(coda.samples, pars = c("beta_Linf", "beta_K", "sigma_Linf", "sigma_K", "sigma_y", "t0p"))


summ<-summary(samples)
head(summ$statistics)
mean_summ_df<-as.data.frame(summ$statistics)
samples.df<-as.data.frame(as.matrix(samples))
