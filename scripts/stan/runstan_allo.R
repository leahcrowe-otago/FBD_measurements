library(cmdstanr)
library(posterior)
library(bayesplot)
library(latex2exp)

### get data

age = readRDS(file = 'age_ij.rds')
ageNA = as.matrix(age)
ageuse = ageNA
ageuse[is.na(ageNA)] = 0

obs_bhdf = readRDS(file = 'bhdf_ij.rds')
obsNA_bhdf = as.matrix(obs_bhdf)
obsuse_bhdf = obsNA_bhdf
obsuse_bhdf[is.na(obsNA_bhdf)] = 0

obs_tl = readRDS(file = 'length_ij.rds')
obsNA_tl = as.matrix(obs_tl)
obsuse_tl = obsNA_tl
obsuse_tl[is.na(obsNA_tl)] = 0

n = nrow(age)
m = ncol(age)
o = ncol(age)
m
o

max_ind = rep(NA,n)
for(i in 1:n){
  max_ind[i] = max(which(!is.na(age[i,])))
}

# ----

a = '

data {
  int<lower=0> n; // number of individuals  
  int<lower=0> m; // dimension of matrix
  array[n] int<lower = 0> max_ind; // the maximum number of observations for each individual
  array[n, m] real<lower=0> y; // the total length observations
  array[n, m] real<lower=0> z; // the bh to df length observations
  array[n, m] real<lower=0> age; // the age observations
}

parameters {
  real beta_Ly; 
  real beta_ky;
  vector[n] eps_Ly;  
  vector[n] eps_ky;
  real<lower=0> sigma_Ly;
  real<lower=0> sigma_ky;
  real<lower=0> sigma_y;
  real<lower=0> sigma_z;
  real alpha_0;
  real gamma_0;
  real<lower=0> alpha_1;
  real<lower=0> gamma_1;
  real<lower=0> sigma_a;
  real<lower=0> sigma_g;
  real<lower=0> t0p;
}

transformed parameters {
  vector[n] ky;
  vector[n] Ly;
  vector[n] kz;
  vector[n] Lz;
  
  for(i in 1:n){
    Ly[i] = beta_Ly + eps_Ly[i]*sigma_Ly;
    ky[i] = inv_logit(beta_ky + eps_ky[i]*sigma_ky);
 }
}

model {

  for (i in 1:n) {
    Lz[i] ~ normal(alpha_0 + alpha_1*Ly[i], sigma_a);
    kz[i] ~ normal(gamma_0 + gamma_1*ky[i], sigma_g);
    
    for(j in 1:max_ind[i]){
      y[i,j] ~ normal(Ly[i]*(1-ky[i]^(age[i,j] + t0p)), sigma_y);
      z[i,j] ~ normal(Lz[i]*(1-kz[i]^(age[i,j] + t0p)), sigma_z);
      
    }
  }
  
  beta_Ly ~ normal(0, 1000);
  beta_ky ~ logistic(0, 1);
  
  eps_Ly ~ std_normal();
  eps_ky ~ std_normal();
  
  sigma_Ly ~ student_t(3, 0, 50);
  sigma_ky ~ student_t(3, 0, 5);

  sigma_y ~ student_t(3, 0, 50);
  sigma_z ~ student_t(3, 0, 50);

  alpha_0 ~ std_normal();
  alpha_1 ~ std_normal();
  sigma_a ~ student_t(3, 0, 50);
  
  gamma_0 ~ std_normal();
  gamma_1 ~ std_normal();
  sigma_g ~ student_t(3, 0, 50);

  t0p ~ lognormal(0, 10);

}
'

set_cmdstan_path(path = "C:/Users/leahm/cmdstan-2.34.1")

cat(a, file = paste0(cmdstan_path(),"/thesis/vb/vb_mod_allo.stan"))

file = file.path(cmdstan_path(), "/thesis/vb/vb_mod_allo.stan")

stan_fit = cmdstan_model(file)

init_vb = function(){
  
  beta_Ly = rnorm(1,mean(obsNA_tl, na.rm = TRUE), 0.2)
  beta_ky = rnorm(1, 0, 0.2)
  eps_Ly = rnorm(n)
  eps_ky = rnorm(n)
  sigma_Ly = rlnorm(1)
  sigma_ky = rlnorm(1)
  
  sigma_y = rlnorm(1)
  sigma_z = rlnorm(1)

  alpha_0 = rnorm(1, 0, 0.2)
  alpha_1 = rlnorm(1)
  sigma_a = rlnorm(1)
  
  gamma_0 = rnorm(1, 0, 0.2)
  gamma_1 = rlnorm(1)
  sigma_g = rlnorm(1)
  
  t0p = rlnorm(1)
  
  return(list(beta_Ly = beta_Ly, beta_ky = beta_ky, eps_Ly = eps_Ly, eps_ky = eps_ky, sigma_ky = sigma_ky, sigma_Ly = sigma_Ly, 
              sigma_y = sigma_y, sigma_z = sigma_z, 
              alpha_0 = alpha_0, alpha_1 = alpha_1, sigma_a = sigma_a, 
              gamma_0 = gamma_0, gamma_1 = gamma_1, sigma_g = sigma_g, 
              t0p = t0p))
}

fit_vb <- stan_fit$sample(
  data = list(
    n = n,
    m = m,
    max_ind = max_ind,
    y = obsuse_tl,
    z = obsuse_bhdf,
    age = ageuse
  ),
  init = init_vb,
  chains = 1,
  iter_warmup = 1000,
  iter_sampling = 5000,
  thin = 1,
  save_warmup = FALSE,
  max_treedepth = 10,
  parallel_chains = 4,
  refresh = 100
)

parout = as_draws_df(fit_vb$draws(c("beta_L", "beta_k", "t0p","sigma_y","sigma_L","sigma_k","alpha_0","alpha_1")))
mcmc_trace(parout)+theme_bw()
ggplot2::ggsave("traceplot_allo.png", device = "png")

summary(parout)

Lout = as_draws_df(fit_vb$draws(c("L")))
kout = as_draws_df(fit_vb$draws(c("k")))
#bout = as_draws_df(fit_vb$draws(c("b")))

### plot a couple of individuals
induse = c(1, 12, 14, 50)
ngrid = 101
agegrid = seq(from = 0, to = max(ageuse), length.out = ngrid)
nind = length(induse)
pdf('indplots_bhdf.pdf', height = 8, width = 8)
par(mfrow = c(2,2), mar = c(4, 4, 1, 1))
for(i in 1:nind){
  plot(ageNA[induse[i],], obsNA[induse[i],], pch = 20, xlim = c(0,max(ageuse)), ylim = c(0, max(obsuse)), xlab = "Age", ylab = "Length")
  
  tmp = matrix(NA,nrow(Lout), ngrid)
  for(j in 1:ngrid){
    tmp[,j] = Lout[[induse[i]]]*(1-kout[[induse[i]]]^(parout[["t0p"]] + agegrid[j]))
  }  
  quan = apply(tmp, 2, quantile, c(0.05, 0.5, 0.95))
  
  lines(agegrid, quan[2,], col = "blue", lty = 1)  
  lines(agegrid, quan[1,], col = "blue", lty = 2)
  lines(agegrid, quan[3,], col = "blue", lty = 2)
}
dev.off()

# mcmc_intervals(kout, outer_size = 0.5, inner_size = 1, point_size = 2)

ageNA[induse[4],]

## BHDF ----
kout_bhdf<-summary(kout)
Lout_bhdf<-summary(Lout)
ID_i<-readRDS("BHDF_ij_ID.rds")

nrow(kout_bhdf)
nrow(Lout_bhdf)
nrow(ID_i)

k_bhdf<-kout_bhdf%>%
  bind_cols(ID = ID_i)

L_bhdf<-Lout_bhdf%>%
  bind_cols(ID = ID_i)

## TL ----
kout_tl<-summary(kout)
Lout_tl<-summary(Lout)
ID_i<-readRDS("length_ij_ID.rds")

nrow(kout_tl)
nrow(Lout_tl)
nrow(ID_i)

k_tl<-kout_tl%>%
  bind_cols(ID = ID_i)

L_tl<-Lout_tl%>%
  bind_cols(ID = ID_i)


## compare ----

young_birthyear<-2019

kcompare<-k_tl%>%
  left_join(lifehist, by = c(ID = 'NAME'))%>%
  left_join(k_bhdf, by = "ID")%>%
  mutate(k_diff = mean.x/mean.y)%>%
  mutate(AGE = case_when(
    BIRTH_YEAR >= young_birthyear ~ "Young",
    TRUE ~ "Old"
  ))%>%
  dplyr::select(ID, mean.x, mean.y, k_diff, SEX, AGE)%>%
  mutate(mean_k_diff = mean(k_diff),
         sd_k_diff = sd(k_diff))

kcompare%>%
 
  dplyr::select(ID, mean.x, mean.y, k_diff, SEX)%>%
  group_by(SEX)%>%
  mutate(mean_k_diff_sex = mean(k_diff),
         sd_k_diff_sex = sd(k_diff))

kcompare%>%
  arrange(k_diff)

ggplot(kcompare)+
  geom_boxplot(aes(x = AGE,y = k_diff, color = AGE))+
  geom_point(aes(x = AGE, y = k_diff, color = AGE))

## bhdf

lcompare<-L_tl%>%
  left_join(lifehist, by = c(ID = 'NAME'))%>%
  left_join(L_bhdf, by = "ID")%>%
  mutate(L_diff = mean.x/mean.y)%>%
  mutate(AGE = case_when(
    BIRTH_YEAR >= young_birthyear ~ "Young",
    TRUE ~ "Old"
  ))%>%
  dplyr::select(ID, mean.x, mean.y, L_diff, SEX, AGE)%>%
  mutate(mean_L_diff = mean(L_diff),
         sd_L_diff = sd(L_diff))

lcompare%>%
  dplyr::select(ID, mean.x, mean.y, L_diff, SEX)%>%
  group_by(SEX)%>%
  mutate(mean_L_diff_sex = mean(L_diff),
         sd_L_diff_sex = sd(L_diff))



ggplot(lcompare)+
  geom_boxplot(aes(x = "L_tl/L_bhdf", y = L_diff, color = SEX))

ggplot(lcompare)+
  geom_boxplot(aes(x = AGE, y = L_diff, color = AGE))+
  geom_point(aes(x = AGE, y = L_diff, color = AGE))

         