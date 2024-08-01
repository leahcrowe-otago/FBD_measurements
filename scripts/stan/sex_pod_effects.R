library(cmdstanr)
library(posterior)
library(bayesplot)
library(dplyr)
library(ggplot2)

# data ----
#IDs
ij_ID = readRDS(file = './data/Measurements/ij_ID.rds')
n_ind = nrow(ij_ID)

#both bhdf & tl
ij_b = readRDS(file = './data/Measurements/ij_1.rds')
#bhdf only
ij_z = readRDS(file = './data/Measurements/ij_2.rds')
#total length only
ij_y = readRDS(file = './data/Measurements/ij_3.rds')

ij_all<-ij_b%>%
  bind_rows(ij_z)%>%
  bind_rows(ij_y)

### sex/pod model ----
obs_sex<-ij_all%>%
  left_join(ij_ID, by = c("ind"))%>%
  mutate(type = case_when(
    !is.na(length) & !is.na(BHDF) ~ 1,
    is.na(length) & !is.na(BHDF) ~ 2,
    !is.na(length) & is.na(BHDF) ~ 3
  ))%>%
  arrange(SEX,type,ID,obs)

n_k<-obs_sex%>%filter(SEX != "X")%>%distinct(ID)%>%nrow()
n_obs_k<-obs_sex%>%filter(SEX != "X")%>%nrow()
n_obs = nrow(obs_sex)

sex<-ij_ID%>%
  filter(SEX != "X")%>%
  dplyr::select(SEX)%>%
  mutate(SEX = case_when(
    SEX == "F" ~ 0,
    SEX == "M" ~ 1))
sex$SEX
  
pod = ij_ID%>%
  dplyr::select(POD)%>%
  mutate(POD = case_when(
    POD == "DOUBTFUL" ~ 0,
    POD == "DUSKY" ~ 1))
pod$POD
#params matrix dims
J=4

# model ----
# z matrix are the individual deviations from the population parameters

set_cmdstan_path(path = "C:/Users/leahm/cmdstan-2.34.1")

stan_fit = cmdstan_model("./scripts/stan/vb_mod_allo_t0_sexpod.stan") #sex/pod effects model

# initial values ----

init_vb = function(){
  
  mu = c(rnorm(1, mean(ij_b$length, na.rm = TRUE), 0.2),
         rlnorm(1, 0, 0.1),
         rnorm(1, mean(ij_b$BHDF, na.rm = TRUE), 0.2),
         rlnorm(1, 0, 0.1))
  sigma =  rlnorm(J)
  L = diag(J)
  z = matrix(rnorm(n_ind*J), n_ind, J)
  sigma_obs = rlnorm(2)
  rho_obs = runif(1)

  beta = rnorm(J,0,0.1)
  gamma = rnorm(J,0,0.1)
  pi = runif(1)
  t0p = rlnorm(1)
  
  return(list(
    beta = beta, gamma = gamma, pi = pi,
    t0p = t0p, L = L, mu = mu, sigma = sigma, z = z, 
    sigma_obs = sigma_obs, rho_obs = rho_obs
  ))

  
}

# run stan ----

fit_vb <- stan_fit$sample(
  data = list(
    
    ### both models ----
    n_ind = n_ind,
    J = J,
    ### sex/pod model ----
    n_k = n_k,
    n_obs_k = n_obs_k,
    n_obs = n_obs,
    id = as.integer(obs_sex$ind),
    age = obs_sex$age,
    type = obs_sex$type,
    val = as.matrix(obs_sex%>%dplyr::select(length, BHDF)%>%replace(is.na(.), 9999)),
    sex = sex$SEX,
    pod = pod$POD
    
  ),
  init = init_vb,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 10000,
  thin = 1,
  save_warmup = FALSE,
  max_treedepth = 10,
  parallel_chains = 4,
  refresh = 100
)

# results -----

parout = as_draws_df(fit_vb$draws(c(#sex/pod model
  "beta","gamma","pi","t0p",
  #both models
  "mu","sigma","sigma_obs","rho_obs",
  "corr[1,2]","corr[1,3]","corr[1,4]",
  "corr[2,3]","corr[2,4]",
  "corr[3,4]",
  "Lobs[1,1]","Lobs[2,1]","Lobs[2,2]",
  "mu_pred", "varcov_par"
)))

saveRDS(parout, file = paste0("./results/parout_",Sys.Date(),"sexpod.rds"))

#trace plot
mcmc_trace(parout)+theme_bw()
cater<-bayesplot::mcmc_intervals(parout, pars = c("beta[1]","beta[2]","beta[3]","beta[4]",
                                                  "gamma[1]","gamma[2]","gamma[3]","gamma[4]"))
bayesplot::mcmc_dens(parout, pars = c("beta[1]","beta[2]","beta[3]","beta[4]","gamma[1]","gamma[2]","gamma[3]","gamma[4]"))
ggplot2::ggsave(paste0("./Figures/sexpod_model.png"), cater, device = "png", dpi = 700, width = 100, height = 100, units = 'mm')
#summary
as.data.frame(summary(parout))

#draws of par from posterior
#sex/pod model
par_kindout = as_draws_df(fit_vb$draws(c("par_k")))
par_u_0indout = as_draws_df(fit_vb$draws(c("par_u_0")))
par_u_1indout = as_draws_df(fit_vb$draws(c("par_u_1")))
#both

#save par for individual plotting
saveRDS(par_kindout, file = paste0("./results/par_kindout_",Sys.Date(),".rds"))
saveRDS(par_u_0indout, file = paste0("./results/par_u_0indout_",Sys.Date(),".rds"))
saveRDS(par_u_1indout, file = paste0("./results/par_u_1indout_",Sys.Date(),".rds"))

## Results ----
date = "2024-06-21" # from runstan_allo_mv_sex, Appendix II
paroutin_sexpod = readRDS(file = paste0('./results/parout_',date,'sexpod.rds'))
summ_paroutin_sexpod<-as.data.frame(summary(paroutin_sexpod))
summ_paroutin_sexpod%>%
  filter(grepl("mu",variable))

#### Table S3 ----
summ_sexpod<-summ_paroutin_sexpod%>%
  mutate(Median = round(median, 3),
         '5th percentile' = format(round(q5,3)),
         '95th percentile' = format(round(q95,3)))%>%
  dplyr::select(variable, Median,'5th percentile','95th percentile')

saveRDS(summ_sexpod, file = "./results/summtable_sexpod_t0.rds")

hist(paroutin_sexpod$`beta[1]`)
hist(paroutin_sexpod$`beta[2]`)
hist(paroutin_sexpod$`beta[3]`)
hist(paroutin_sexpod$`beta[4]`)

hist(paroutin_sexpod$`gamma[1]`)
hist(paroutin_sexpod$`gamma[2]`)
hist(paroutin_sexpod$`gamma[3]`)
hist(paroutin_sexpod$`gamma[4]`)

parindout_in_k = readRDS(file = paste0('./results/par_kindout_',date,'.rds'))
parindout_in_u0 = readRDS(file = paste0('./results/par_u_0indout_',date,'.rds'))
parindout_in_u1 = readRDS(file = paste0('./results/par_u_1indout_',date,'.rds'))
