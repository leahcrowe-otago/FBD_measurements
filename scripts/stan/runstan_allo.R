library(cmdstanr)
library(posterior)
library(bayesplot)
library(latex2exp)
i = 52
### get data
age = readRDS(file = 'age_ijboth.rds')
ageNA = as.matrix(age)
ageuse = ageNA
ageuse[is.na(ageNA)] = 0
ageuse[i,]

obs_bhdf = readRDS(file = 'bhdf_ijboth.rds')
obsNA_bhdf = as.matrix(obs_bhdf)
obsuse_bhdf = obsNA_bhdf
obsuse_bhdf[is.na(obsNA_bhdf)] = 0
obsuse_bhdf[i, ]

obs_tl = readRDS(file = 'length_ijboth.rds')
obsNA_tl = as.matrix(obs_tl)
obsuse_tl = obsNA_tl
obsuse_tl[is.na(obsNA_tl)] = 0
obsuse_tl[i, ]

n = nrow(age)
m = ncol(age)
o = ncol(age)
m
o



max_ind = rep(NA,n)
for(i in 1:n){
  max_ind[i] = max(which(!is.na(age[i,])))
}

obsuse_tl[1:max_ind[i]]
max_ind[i]
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
  vector[n] eps_Lz;
  vector[n] eps_kz;
  real<lower=0> sigma_Ly;
  real<lower=0> sigma_ky;
  real<lower=0> sigma_y;
  real<lower=0> sigma_z;
  real alpha_0;
  real gamma_0;
  real alpha_1;
  real gamma_1;
  real<lower=0> sigma_a;
  real<lower=0> sigma_g;
  real<lower=0> t0p;
  //real<lower=0> t0p_z;


}

transformed parameters {
  vector[n] ky_logit;
  vector[n] Ly;
  vector[n] kz_logit;
  vector[n] Lz;
  
  for(i in 1:n){
    Ly[i] = beta_Ly + eps_Ly[i]*sigma_Ly;
    ky_logit[i] = beta_ky + eps_ky[i]*sigma_ky;
    Lz[i] = alpha_0 + alpha_1*Ly[i] + eps_Lz[i]*sigma_a;
    kz_logit[i] = gamma_0 + gamma_1*ky_logit[i] + eps_kz[i]*sigma_g;
  }

}

model {
  vector[n] ky;
  vector[n] kz;

  for (i in 1:n) {
    
    ky[i] = inv_logit(ky_logit[i]);
    kz[i] = inv_logit(kz_logit[i]);
    
    for(j in 1:max_ind[i]){
      y[i,j] ~ normal(Ly[i]*(1-ky[i]^(age[i,j] + t0p)), sigma_y);
      z[i,j] ~ normal(Lz[i]*(1-kz[i]^(age[i,j] + t0p)), sigma_z);
      
    }
  }
  
  beta_Ly ~ normal(0, 1000);
  beta_ky ~ logistic(0, 1);
  
  eps_Ly ~ std_normal();
  eps_Lz ~ std_normal();
  eps_ky ~ std_normal();
  eps_kz ~ std_normal();
  
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
  //t0p_z ~ lognormal(0, 10);

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
  eps_Lz = rnorm(n)
  eps_kz = rnorm(n)
  sigma_Ly = rlnorm(1)
  sigma_ky = rlnorm(1)
  sigma_Lz = rlnorm(1)
  sigma_kz = rlnorm(1)
  
  sigma_y = rlnorm(1)
  sigma_z = rlnorm(1)

  alpha_0 = rnorm(1, 0, 0.2)
  alpha_1 = rnorm(1,mean(obsNA_bhdf, na.rm = TRUE)/mean(obsNA_tl, na.rm = TRUE), 0.2)
  sigma_a = rlnorm(1)
  
  gamma_0 = rnorm(1, 0, 0.2)
  gamma_1 = rnorm(1, 0, 0.2)
  sigma_g = rlnorm(1)
  
  t0p = rlnorm(1)
  #t0p_z = rlnorm(1)

  return(list(beta_Ly = beta_Ly, beta_ky = beta_ky, 
              eps_Ly = eps_Ly, eps_ky = eps_ky, eps_Lz = eps_Lz, eps_kz = eps_kz, 
              sigma_Ly = sigma_Ly, sigma_ky = sigma_ky, sigma_Lz = sigma_Lz, sigma_kz = sigma_kz,
              sigma_y = sigma_y, sigma_z = sigma_z, 
              alpha_0 = alpha_0, alpha_1 = alpha_1, sigma_a = sigma_a, 
              gamma_0 = gamma_0, gamma_1 = gamma_1, sigma_g = sigma_g, 
              t0p = t0p))#, t0p_z = t0p_z))
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
parout = as_draws_df(fit_vb$draws(c("beta_Ly", "beta_ky","t0p","sigma_y","sigma_Ly","sigma_ky","alpha_0","alpha_1","sigma_a", "gamma_0", "gamma_1", "sigma_g")))
mcmc_trace(parout)+theme_bw()
ggplot2::ggsave("traceplot_allo.png", device = "png", dpi = 300, height = 200, width = 300, units = 'mm')

summary(parout)

Lyout = as_draws_df(fit_vb$draws(c("Ly")))
kyout_logit = as_draws_df(fit_vb$draws(c("ky_logit")))
Lzout = as_draws_df(fit_vb$draws(c("Lz")))
kzout_logit = as_draws_df(fit_vb$draws(c("kz_logit")))
#bout = as_draws_df(fit_vb$draws(c("b")))

mcmc_intervals(Lyout) + theme_minimal()
mcmc_intervals(Lyout[110]) + theme_minimal()
mcmc_intervals(Lzout) + theme_minimal()
mcmc_intervals(Lzout[110]) + theme_minimal()

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

ggplot()

# mcmc_intervals(kout, outer_size = 0.5, inner_size = 1, point_size = 2)

ageNA[induse[4],]

## BHDF ----
kyout_bhdf<-summary(kyout)
Lyout_bhdf<-summary(Lyout)
ID_i<-readRDS("BHDF_ij_ID.rds")

nrow(kyout_bhdf)
nrow(Lyout_bhdf)
nrow(ID_i)
library(dplyr)
ky_bhdf<-kyout_bhdf%>%
  bind_cols(ID = ID_i)

Ly_bhdf<-Lyout_bhdf%>%
  bind_cols(ID = ID_i)

## TL ----
kyout_tl<-summary(kyout)
Lyout_tl<-summary(Lyout)
ID_i<-readRDS("length_ij_IDboth.rds")
ID_i[110,]
nrow(kyout_tl)
nrow(Lyout_tl)
nrow(ID_i)

ky_tl<-kyout_tl%>%
  bind_cols(ID = ID_i)

Ly_tl<-Lyout_tl%>%
  bind_cols(ID = ID_i)


## compare ----

young_birthyear<-2019

kcompare<-ky_tl%>%
  left_join(lifehist, by = c(ID = 'NAME'))%>%
  left_join(ky_bhdf, by = "ID")%>%
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


dnorm(-0.0599714+(2.86325*2.70193), 1.28604)
log(-0.06)
