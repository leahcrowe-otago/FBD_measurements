library(cmdstanr)
library(posterior)
library(bayesplot)
library(latex2exp)
library(dplyr)
i = 18
### get data
#both bhdf & tl
ij_b = readRDS(file = 'ij_1.rds')
#ij_b = as.matrix(ij_1)
ij_b[i,]
N_b = nrow(ij_b)
N_b
#bhdf only
ij_z = readRDS(file = 'ij_2.rds')%>%
  dplyr::select(-length)
#ij_z = as.matrix(ij_z)
ij_z[i,]
N_z = nrow(ij_z)
N_z

#tl only
ij_y = readRDS(file = 'ij_3.rds')%>%
  dplyr::select(-BHDF)
#ij_y = as.matrix(ij_y)
ij_y[i,]
N_y = nrow(ij_y)
N_y

#IDs
ij_ID = readRDS(file = 'ij_ID.rds')

n_ind<-nrow(ij_ID)

#params matrix dims
J=4

# ----

a = "

data {
  int<lower=0> n_ind; // number individuals
  int<lower=0> N_b; // number of times when tl and bhdf measurements were both taken
  int<lower=0> N_z; // number of times when only bhdf measurements were taken
  int<lower=0> N_y; // number of times when only tl measurements were taken
  int<lower=0> J; // number of params (Ly, ky, Lz, kz), dim of multivariate matrix
  
  array[N_b] int ind_b; //the individuals
  array[N_y] int ind_y; //the individuals
  array[N_z] int ind_z; //the individuals
  
  vector[N_b] age_b; // the age observations when tl and bhdf measurements were both taken
  vector[N_b] y_b; // the total length observations
  vector[N_b] z_b; // the bhdf length observations
  
  vector[N_y] age_y; // the age observations when only tl measurements were taken
  vector[N_y] y_y; // the total length observations when only tl measurements were taken
  
  vector[N_z] age_z; // the age observations when only bhdf measurements were taken
  vector[N_z] z_z; // the bhdf length observationswhen only bhdf measurements were taken
  
}

parameters {
  real<lower=0> t0p;
  matrix[n_ind, 4] par;
  vector<lower=0>[J] sigma;
  cholesky_factor_corr[J] L;
  vector[J] mu;
  vector<lower=0>[2] sigma_obs;

}

transformed parameters{

}


model {


 for(i in 1:n_ind){
      par[i] ~ multi_normal_cholesky(mu,diag_pre_multiply(sigma,L));
      }


 for (i in 1:N_y){
      y_y[i] ~ normal(par[ind_y[i],1]*(1-inv_logit(par[ind_y[i],2])^(age_y[i] + t0p)), sigma_obs[1]);
      }
      
 for (i in 1:N_z) {
      z_z[i] ~ normal(par[ind_z[i],3]*(1-inv_logit(par[ind_z[i],4])^(age_z[i] + t0p)), sigma_obs[2]);
      } 
    
 for (i in 1:N_b){
      y_b[i] ~ normal(par[ind_b[i],1]*(1-inv_logit(par[ind_b[i],2])^(age_b[i] + t0p)), sigma_obs[1]);
      z_b[i] ~ normal(par[ind_b[i],3]*(1-inv_logit(par[ind_b[i],4])^(age_b[i] + t0p)), sigma_obs[2]);
      }


//priors
  t0p ~ lognormal(0, 10);
  mu ~ normal(0, 1000);
  sigma ~ student_t(3, 0, 25);
  L ~ lkj_corr_cholesky(1);
  
  for (i in 1:2){
  sigma_obs[i] ~ student_t(3, 0, 50);
  }

}


generated quantities {

  corr_matrix[J] Sigma;
  Sigma = multiply_lower_tri_self_transpose(L);

}

"

set_cmdstan_path(path = "C:/Users/leahm/cmdstan-2.34.1")

cat(a, file = paste0(cmdstan_path(),"/thesis/vb/vb_mod_allo.stan"))

file = file.path(cmdstan_path(), "/thesis/vb/vb_mod_allo.stan")

stan_fit = cmdstan_model(file)

init_vb = function(){
  
  t0p = rlnorm(1)
  
  mu = c(rnorm(1,mean(ij_NA$length, na.rm = TRUE), 0.2),
         rlnorm(1, 0, 0.1),
         rnorm(1,mean(ij_NA$BHDF, na.rm = TRUE), 0.2),
         rlnorm(1, 0, 0.1))
  
  sigma = rlnorm(J)
  
  L = diag(J)
  
  par = matrix(NA,n_ind, J)
  for(j in 1:J){
    if(j == 1 | j == 3){
      par[,j] = rnorm(n_ind,mu[j], 0.2)  
    } else {
      par[,j] = rlnorm(n_ind, log(mu[j]), 0.1)
    }
  } 
  
  sigma_obs = rlnorm(2)




  return(list(t0p = t0p, 
              L = L, mu = mu, sigma = sigma, par = par, 
              sigma_obs = sigma_obs
              ))
}

fit_vb <- stan_fit$sample(
  data = list(
    
    ind_b = ij_b$ind,
    y_b = ij_b$length,
    z_b = ij_b$BHDF,
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
    
    n_ind = n_ind,
    J = J
    
  ),
  init = init_vb,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 5000,
  thin = 1,
  save_warmup = FALSE,
  max_treedepth = 10,
  parallel_chains = 4,
  refresh = 100
)


library(ggplot2)
parout = as_draws_df(fit_vb$draws(c("mu","sigma","sigma_obs","t0p")))
mcmc_trace(parout)+theme_bw()
ggplot2::ggsave("traceplot_allo.png", device = "png", dpi = 300, height = 200, width = 300, units = 'mm')

summary(parout)

Lyout = as_draws_df(fit_vb$draws(c("Ly")))
kyout_logit = as_draws_df(fit_vb$draws(c("ky_logit")))
Lzout = as_draws_df(fit_vb$draws(c("Lz")))
kzout_logit = as_draws_df(fit_vb$draws(c("kz_logit")))


### plot a couple of individuals
induse = c(1, 12, 14, 50)
ngrid = 101
agegrid = seq(from = 0, to = max(ageuse), length.out = ngrid)
nind = length(induse)
pdf('indplots_bhdf.pdf', height = 8, width = 8)
par(mfrow = c(2,2), mar = c(4, 4, 1, 1))
for(i in 1:nind){

  plot(ageNA[induse[i],], obsNA_tl[induse[i],], pch = 20, xlim = c(0,max(ageuse)), ylim = c(0, max(obsuse_tl)), xlab = "Age", ylab = "Length")
  points(ageNA[induse[i],], obsNA_bhdf[induse[i],], pch = 20, xlim = c(0,max(ageuse)), ylim = c(0, max(obsuse_tl)), xlab = "Age", ylab = "Length")
  
  tmp_y = matrix(NA,nrow(Lyout), ngrid)
  tmp_z = matrix(NA,nrow(Lzout), ngrid)
  for(j in 1:ngrid){
    #inverse logit kyout/kzout #exp(x)/(1+exp(x))
    tmp_y[,j] = Lyout[[induse[i]]]*(1-(exp(kyout_logit[[induse[i]]])/(1+exp(kyout_logit[[induse[i]]])))^(parout[["t0p"]] + agegrid[j]))
    tmp_z[,j] = Lzout[[induse[i]]]*(1-(exp(kzout_logit[[induse[i]]])/(1+exp(kzout_logit[[induse[i]]])))^(parout[["t0p"]] + agegrid[j]))
      }  
  quan_y = apply(tmp_y, 2, quantile, c(0.05, 0.5, 0.95))
  quan_z = apply(tmp_z, 2, quantile, c(0.05, 0.5, 0.95))
  
  lines(agegrid, quan_y[2,], col = "blue", lty = 1)  
  lines(agegrid, quan_y[1,], col = "blue", lty = 2)
  lines(agegrid, quan_y[3,], col = "blue", lty = 2)
  
  lines(agegrid, quan_z[2,], col = "red", lty = 1)  
  lines(agegrid, quan_z[1,], col = "red", lty = 2)
  lines(agegrid, quan_z[3,], col = "red", lty = 2)
}
dev.off()


kyout<-exp(kyout_logit)/(1+exp(kyout_logit))
kyout$
summary(kyout)

mcmc_intervals(Lyout, outer_size = 0.5, inner_size = 1, point_size = 2)
mcmc_intervals(Lzout, outer_size = 0.5, inner_size = 1, point_size = 2)
mcmc_intervals(kyout_logit, outer_size = 0.5, inner_size = 1, point_size = 2)
mcmc_intervals(kzout_logit, outer_size = 0.5, inner_size = 1, point_size = 2)

##########

