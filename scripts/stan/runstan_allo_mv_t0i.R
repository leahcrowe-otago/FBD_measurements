library(cmdstanr)
library(posterior)
library(bayesplot)
library(dplyr)
library(ggplot2)

# data ----
i = 34
#IDs
ij_ID = readRDS(file = './data/Measurements/ij_ID.rds')
n_ind = nrow(ij_ID)

### reg model ----
#both bhdf & tl
ij_b = readRDS(file = './data/Measurements/ij_1.rds')
ij_z = readRDS(file = './data/Measurements/ij_2.rds')
ij_y = readRDS(file = './data/Measurements/ij_3.rds')

ij_all<-ij_b%>%
  bind_rows(ij_z)%>%
  bind_rows(ij_y)

median(ij_all$length, na.rm = T)
median(ij_all$BHDF, na.rm = T)

#ij_b = as.matrix(ij_1)
ij_b[i,]
N_b = nrow(ij_b)
N_b
length(unique(ij_b$ind))
#bhdf only
ij_z = ij_z%>%
  dplyr::select(-length)
#ij_z = as.matrix(ij_z)
ij_z[i,]
N_z = nrow(ij_z)
N_z

#tl only
ij_y = ij_y%>%
  dplyr::select(-BHDF)
#ij_y = as.matrix(ij_y)
ij_y[i,]
N_y = nrow(ij_y)
N_y

N_b+N_z+N_y
N_b/(N_b+N_z+N_y)

### sex/pod model ----
obs_sex = readRDS(file = './data/Measurements/ij_obs_sex.rds')
n_k<-obs_sex%>%filter(SEX != "X")%>%distinct(ID)%>%nrow()
n_obs_k<-obs_sex%>%filter(SEX != "X")%>%nrow()
n_obs = nrow(obs_sex)
sex = ij_ID%>%
  dplyr::select(SEX)%>%
  filter(SEX != "X")%>%
  mutate(SEX = case_when(
    SEX == "F" ~ 0,
    SEX == "M" ~ 1))
pod = ij_ID%>%
  dplyr::select(POD)%>%
  mutate(POD = case_when(
    POD == "DOUBTFUL" ~ 0,
    POD == "DUSKY" ~ 1))
#params matrix dims
J=4

# model ----
# z matrix and z_t are the individual deviations from the population parameters

set_cmdstan_path(path = "C:/Users/leahm/cmdstan-2.34.1")

stan_fit = cmdstan_model("./scripts/stan/vb_mod_all_t0i.stan") #reg model
#stan_fit = cmdstan_model("./scripts/stan/vb_mod_all_t0i_sex.stan") #sex/pod effects model

# initial values ----

init_vb = function(){
  
  mu_t = rlnorm(1)
  z_t = rnorm(n_ind)
  t0p = rlnorm(n_ind, log(mu_t), 0.1)
  sigma_t = rlnorm(1, 0, 0.1)
  
  mu = c(rnorm(1, mean(ij_b$length, na.rm = TRUE), 0.2),
         rlnorm(1, 0, 0.1),
         rnorm(1, mean(ij_b$BHDF, na.rm = TRUE), 0.2),
         rlnorm(1, 0, 0.1))
  sigma = rlnorm(J)
  L = diag(J)
  z = matrix(rnorm(n_ind*J), n_ind, J)
  sigma_obs = rlnorm(2)
  rho_obs = runif(1)
  
  #sex/pod
  beta = rlnorm(J)
  gamma = rlnorm(J)
  pi = runif(1)
    
  return(list(#obs/sex
              beta = beta, gamma = gamma, pi = pi,
              #both versions
              z_t = z_t, mu_t = mu_t, t0p = t0p, sigma_t = sigma_t,
              L = L, mu = mu, sigma = sigma, z = z, 
              sigma_obs = sigma_obs, rho_obs = rho_obs
              ))
}

# run stan ----

fit_vb <- stan_fit$sample(
  data = list(
    
  ### both models ----
    n_ind = n_ind,
    J = J,
  ### reg model ----
    ind_b = ij_b$ind,
    data_b = ij_b%>%dplyr::select(length, BHDF),
    age_b = ij_b$age,
    
    ind_y = ij_y$ind,
    y_y = ij_y$length,
    age_y = ij_y$age,
    
    ind_z = ij_z$ind,
    z_z = ij_z$BHDF,
    age_z = ij_z$age,
    
    N_b = N_b,
    N_z = N_z,
    N_y = N_y,
  ### sex/pod model ----
    n_k = n_k,
    n_obs_k = n_obs_k,
    n_obs = n_obs,
    id = obs_sex$ind,
    age = obs_sex$age,
    type = obs_sex$type,
    val = obs_sex%>%dplyr::select(length, BHDF)%>%replace(is.na(.), 0),
    sex = sex$SEX,
    pod = pod$POD
  
  ),
  init = init_vb,
  chains = 4,
  iter_warmup = 1500,
  iter_sampling = 15000,
  thin = 1,
  save_warmup = FALSE,
  max_treedepth = 10,
  parallel_chains = 4,
  refresh = 100
)

# results -----

parout = as_draws_df(fit_vb$draws(c(#sex/pod model
                                    #"beta","gamma","pi",
                                    #both models
                                    "mu","mu_t","sigma","sigma_t","sigma_obs","rho_obs",
                                    "varcov_par",
                                    "corr[1,2]","corr[1,3]","corr[1,4]",
                                    "corr[2,3]","corr[2,4]",
                                    "corr[3,4]",
                                    "Lobs[1,1]","Lobs[2,1]","Lobs[2,2]",
                                    "mu_pred")))

saveRDS(parout, file = paste0("./results/parout_",Sys.Date(),".rds"))
#saveRDS(parout, file = paste0("./results/parout_",Sys.Date(),"sexpod.rds"))

#trace plot
mcmc_trace(parout)+theme_bw()

#summary
as.data.frame(summary(parout))

#draws of par from posterior
parindout = as_draws_df(fit_vb$draws(c("par")))
t0pindout = as_draws_df(fit_vb$draws(c("t0p")))
#save par for individual plotting
saveRDS(parindout, file = paste0("./results/parindout_",Sys.Date(),".rds"))
saveRDS(t0pindout, file = paste0("./results/t0pindout_",Sys.Date(),".rds"))
