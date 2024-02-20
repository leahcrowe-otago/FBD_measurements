library(cmdstanr)
library(posterior)
library(bayesplot)
library(latex2exp)
i = 785
### get data
#both bhdf & tl
ij_b = readRDS(file = 'ij_1.rds')
#ij_b = as.matrix(ij_1)
ij_b[i,]
N_b = nrow(ij_b)
N_b
#bhdf only
ij_z = readRDS(file = 'ij_2.rds')
#ij_z = as.matrix(ij_z)
ij_z[i,]
N_z = nrow(ij_z)
N_z

#tl only
ij_y = readRDS(file = 'ij_3.rds')
#ij_y = as.matrix(ij_y)
ij_y[i,]
N_y = nrow(ij_y)
N_y

ij_NA<-ij_z%>%
  bind_rows(ij_b)%>%
  bind_rows(ij_y)

ij<-ij_NA
ij[is.na(ij)] = 0

ij[i+1,]
N_ij<-nrow(ij)

n_ind<-length(unique(ij$ind))
max(ij$ind)

#this might have to happen separately?
ij%>%
  group_by(ind)%>%
  mutate(j = max(obs))


ij%>%
  filter(BHDF > 0)%>%
  distinct(ind)%>%
  tally()

# ----

a = "

data {
  int<lower=0> n_ind; // number individuals
  int<lower=0> N_ij; // number of measurement rows
  int<lower=0> N_b; // number of times when tl and bhdf measurements were taken
  int<lower=0> N_z; // number of times when only bhdf measurements were taken
  int<lower=0> N_y; // number of times when only tl measurements were taken
  array[N_ij] int ind; //the individuals
  vector[N_ij] age; // the age observations
  vector[N_ij] y; // the total length observations
  vector[N_ij] z; // the bhdf length observations
}

parameters {
  real<lower=0> t0p;
  matrix[n_ind, 4] par;
  vector<lower=0>[n_ind] sigma;
  cholesky_factor_corr[n_ind] L;
  vector[n_ind] mu;


}


model {

 for(i in 1:n_ind){
    
    for (j in 1:4){
      par[i,j] ~ multi_normal_cholesky(mu[i], diag_matrix(sigma[i])*L);
    }
 }
 
  
 for(i in 1:N_ij){
    if (i <= N_z) {
      z[i] ~ normal(par[ind[i],3]*(1-inv_logit(par[ind[i],4])^(age[i] + t0p)), sigma[ind[i]]);
    
    } else if (i > N_b+N_z){
      y[i] ~ normal(par[ind[i],1]*(1-inv_logit(par[ind[i],2])^(age[i] + t0p)), sigma[ind[i]]);
      
    } else {
      y[i] ~ normal(par[ind[i],1]*(1-inv_logit(par[ind[i],2])^(age[i] + t0p)), sigma[ind[i]]);
      z[i] ~ normal(par[ind[i],3]*(1-inv_logit(par[ind[i],4])^(age[i] + t0p)), sigma[ind[i]]);
      
    }
    }

  t0p ~ lognormal(0, 10);
  mu ~ normal(0, 1000);
  sigma ~ student_t(3, 0, 25);
  L ~ lkj_corr_cholesky(1);

}


generated quantities {


}

"

set_cmdstan_path(path = "C:/Users/leahm/cmdstan-2.34.1")

cat(a, file = paste0(cmdstan_path(),"/thesis/vb/vb_mod_allo.stan"))

file = file.path(cmdstan_path(), "/thesis/vb/vb_mod_allo.stan")

stan_fit = cmdstan_model(file)

init_vb = function(){
  
  beta_Ly = rnorm(1,mean(ij_NA$length, na.rm = TRUE), 0.2)
  beta_ky = rnorm(1, 0, 0.2)

  eps_Ly = rnorm(n_ind)
  eps_ky = rnorm(n_ind)
  eps_Lz = rnorm(n_ind)
  eps_kz = rnorm(n_ind)
  
  sigma_Ly = rlnorm(1)
  sigma_ky = rlnorm(1)
  sigma_Lz = rlnorm(1)
  sigma_kz = rlnorm(1)
  
  sigma_y = rlnorm(1)
  sigma_z = rlnorm(1)

  alpha_0 = rnorm(1, 0, 0.2)
  alpha_1 = rnorm(1,mean(ij_NA$BHDF, na.rm = TRUE)/mean(ij_NA$length, na.rm = TRUE), 0.2)
  sigma_a = rlnorm(1)
  
  gamma_0 = rnorm(1, 0, 0.2)
  gamma_1 = rnorm(1, 0, 0.2)
  sigma_g = rlnorm(1)
  
  t0p = rlnorm(1)
  
  mu = rnorm(n_ind, mean(ij_NA$length, na.rm = TRUE), 0.2)
  sigma = rlnorm(n_ind)
  L = diag(4)


  return(list(beta_Ly = beta_Ly, beta_ky = beta_ky, 
              eps_Ly = eps_Ly, eps_ky = eps_ky, eps_Lz = eps_Lz, eps_kz = eps_kz, 
              sigma_Ly = sigma_Ly, sigma_ky = sigma_ky, sigma_Lz = sigma_Lz, sigma_kz = sigma_kz,
              sigma_y = sigma_y, sigma_z = sigma_z, 
              alpha_0 = alpha_0, alpha_1 = alpha_1, sigma_a = sigma_a, 
              gamma_0 = gamma_0, gamma_1 = gamma_1, sigma_g = sigma_g, 
              t0p = t0p, 
              L = L, mu = mu, sigma = sigma
              ))
}

fit_vb <- stan_fit$sample(
  data = list(
    
    y = ij$length,
    z = ij$BHDF,
    age = ij$age,
    N_b = N_b,
    N_z = N_z,
    N_y = N_y,
    N_ij = N_ij,
    n_ind = n_ind,
    ind = ij$ind
    
  ),
  init = init_vb,
  chains = 4,
  iter_warmup = 10,
  iter_sampling = 50,
  thin = 1,
  save_warmup = FALSE,
  max_treedepth = 10,
  parallel_chains = 4,
  refresh = 100
)


library(ggplot2)
parout = as_draws_df(fit_vb$draws(c("beta_Ly", "beta_ky","t0p","sigma_y","sigma_Ly","sigma_ky","alpha_0","alpha_1","sigma_Lz", "gamma_0", "gamma_1", "sigma_kz")))
mcmc_trace(parout)+theme_bw()
ggplot2::ggsave("traceplot_allo.png", device = "png", dpi = 300, height = 200, width = 300, units = 'mm')

summary(parout)

Lyout = as_draws_df(fit_vb$draws(c("Ly")))
kyout_logit = as_draws_df(fit_vb$draws(c("ky_logit")))
Lzout = as_draws_df(fit_vb$draws(c("Lz")))
kzout_logit = as_draws_df(fit_vb$draws(c("kz_logit")))
#bout = as_draws_df(fit_vb$draws(c("b")))

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

