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
#bhdf only
ij_z = readRDS(file = './data/Measurements/ij_2.rds')
#total length only
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

# model ----
# z matrix and z_t are the individual deviations from the population parameters

set_cmdstan_path(path = "C:/Users/leahm/cmdstan-2.34.1")

stan_fit = cmdstan_model("./scripts/stan/vb_mod_allo_t0.stan") #reg model

# initial values ----

init_vb = function(){
  
  t0p = rlnorm(1)
  mu = c(rnorm(1, mean(ij_b$length, na.rm = TRUE), 0.2),
         rlnorm(1, 0, 0.1),
         rnorm(1, mean(ij_b$BHDF, na.rm = TRUE), 0.2),
         rlnorm(1, 0, 0.1))
  sigma = rlnorm(J)
  L = diag(J)
  z = matrix(rnorm(n_ind*J), n_ind, J)
  sigma_obs = rlnorm(2)
  rho_obs = runif(1)

  return(list(
    t0p = t0p, L = L, mu = mu, sigma = sigma, z = z,
    sigma_obs = sigma_obs, rho_obs = rho_obs
  ))
}

# run stan ----

fit_vb <- stan_fit$sample(
  data = list(
    
    n_ind = n_ind,
    J = J,

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
    N_y = N_y
  
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

parout = as_draws_df(fit_vb$draws(c(
                                    "mu","sigma","sigma_obs","rho_obs",
                                    "corr[1,2]","corr[1,3]","corr[1,4]",
                                    "corr[2,3]","corr[2,4]",
                                    "corr[3,4]","t0p",
                                    "Lobs[1,1]","Lobs[2,1]","Lobs[2,2]",
                                    "varcov_par","mu_pred"
                                    )))

saveRDS(parout, file = paste0("./results/parout_",Sys.Date(),".rds"))

#trace plot
mcmc_trace(parout)+theme_bw()

#summary
as.data.frame(summary(parout))

#draws of par from posterior
parindout = as_draws_df(fit_vb$draws(c("par")))

#save par for individual plotting
saveRDS(parindout, file = paste0("./results/parindout_",Sys.Date(),".rds"))
