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

a = '

data {
  int<lower=0> n_ind; // number individuals
  int<lower=0> N_ij; // number of measurement rows
  int<lower=0> N_b; // number of timse when tl and bhdf measurements were taken
  int<lower=0> N_z; // number of times when only bhdf measurements were taken
  int<lower=0> N_y; // number of times when only tl measurements were taken
  array[N_ij] int ind;
  vector [N_ij] age; // the age observations
  vector [N_ij] y; // the total length observations
  vector [N_ij] z; // the bh to df length observations
}

parameters {
  real beta_Ly;
  real beta_ky;
  vector[n_ind] eps_Ly;
  vector[n_ind] eps_ky;
  vector[n_ind] eps_Lz;
  vector[n_ind] eps_kz;
  real<lower=0> sigma_Ly;
  real<lower=0> sigma_ky;
  real<lower=0> sigma_y;
  real<lower=0> sigma_z;
  real alpha_0;
  real gamma_0;
  real alpha_1;
  real gamma_1;
  real<lower=0> sigma_Lz;
  real<lower=0> sigma_kz;
  real<lower=0> t0p;

}

transformed parameters {
  vector[n_ind] ky_logit;
  vector[n_ind] Ly;
  vector[n_ind] kz_logit;
  vector[n_ind] Lz;

for(i in 1:N_ij){
    Ly[ind[i]] = beta_Ly + eps_Ly[ind[i]]*sigma_Ly;
    ky_logit[ind[i]] = beta_ky + eps_ky[ind[i]]*sigma_ky;
    Lz[ind[i]] = alpha_0 + alpha_1*Ly[ind[i]] + eps_Lz[ind[i]]*sigma_Lz;
    kz_logit[ind[i]] = gamma_0 + gamma_1*ky_logit[ind[i]] + eps_kz[ind[i]]*sigma_kz;
}
}

model {
  vector[n_ind] ky;
  vector[n_ind] kz;

  for(i in 1:N_ij){    
    
    ky[ind[i]] = inv_logit(ky_logit[ind[i]]);
    kz[ind[i]] = inv_logit(kz_logit[ind[i]]);

    if (i <= N_z) {

      z[i] ~ normal(Lz[ind[i]]*(1-kz[ind[i]]^(age[i] + t0p)), sigma_z);
    
    } else if (i > N_b+N_z){
      y[i] ~ normal(Ly[ind[i]]*(1-ky[ind[i]]^(age[i] + t0p)), sigma_y);
      
    } else {
      y[i] ~ normal(Ly[ind[i]]*(1-ky[ind[i]]^(age[i] + t0p)), sigma_y);
      z[i] ~ normal(Lz[ind[i]]*(1-kz[ind[i]]^(age[i] + t0p)), sigma_z);
      
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
  sigma_Ly ~ student_t(3, 0, 50);
  
  gamma_0 ~ std_normal();
  gamma_1 ~ std_normal();
  sigma_kz ~ student_t(3, 0, 50);

  t0p ~ lognormal(0, 10);
  //t0p_z ~ lognormal(0, 10);

}

generated quantities {


}

'

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

  return(list(beta_Ly = beta_Ly, beta_ky = beta_ky, 
              eps_Ly = eps_Ly, eps_ky = eps_ky, eps_Lz = eps_Lz, eps_kz = eps_kz, 
              sigma_Ly = sigma_Ly, sigma_ky = sigma_ky, sigma_Lz = sigma_Lz, sigma_kz = sigma_kz,
              sigma_y = sigma_y, sigma_z = sigma_z, 
              alpha_0 = alpha_0, alpha_1 = alpha_1, sigma_a = sigma_a, 
              gamma_0 = gamma_0, gamma_1 = gamma_1, sigma_g = sigma_g, 
              t0p = t0p))
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
  iter_warmup = 1000,
  iter_sampling = 5000,
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

# mcmc_intervals(kout, outer_size = 0.5, inner_size = 1, point_size = 2)

##########

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
ID_i<-readRDS("length_ij_ID.rds")

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

as.matrix[n, m]

B[j,1:3]
