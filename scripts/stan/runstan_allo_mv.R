library(cmdstanr)
library(posterior)
library(bayesplot)
library(latex2exp)
library(dplyr)

# data ----

i = 18
#both bhdf & tl
ij_b = readRDS(file = 'ij_1.rds')
#ij_b = as.matrix(ij_1)
ij_b[i,]
N_b = nrow(ij_b)
N_b
length(unique(ij_b$ind))
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

N_b+N_z+N_y
N_b/(N_b+N_z+N_y)

#params matrix dims
J=4

# model ----

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
  matrix[N_b, 2] data_b; // [1] = the total length observations, [2] = the bhdf length observations
  
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
  real<lower=0,upper=1> rho_obs;

  
 }

transformed parameters{

  matrix[2,2] varcov;
  cholesky_factor_cov[2] Lobs;

  varcov[1,1] = sigma_obs[1]^2;
  varcov[2,2] = sigma_obs[2]^2;
  varcov[1,2] = sigma_obs[1]*sigma_obs[2]*rho_obs;
  varcov[2,1] = varcov[1,2];
  Lobs = cholesky_decompose(varcov);

 }

model {

 matrix[N_b,2] obs_mean;

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
      obs_mean[i,1] = par[ind_b[i],1]*(1-inv_logit(par[ind_b[i],2])^(age_b[i] + t0p));
      obs_mean[i,2] = par[ind_b[i],3]*(1-inv_logit(par[ind_b[i],4])^(age_b[i] + t0p));
      data_b[i] ~ multi_normal_cholesky (obs_mean[i], Lobs);

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

  corr_matrix[J] corr;
  //below is L*L'
  corr = multiply_lower_tri_self_transpose(L);

 }

"

set_cmdstan_path(path = "C:/Users/leahm/cmdstan-2.34.1")

cat(a, file = paste0(cmdstan_path(),"/thesis/vb/vb_mod_allo.stan"))

file = file.path(cmdstan_path(), "/thesis/vb/vb_mod_allo.stan")

stan_fit = cmdstan_model(file)

# initial values ----

init_vb = function(){
  
  t0p = rlnorm(1)
  
  mu = c(rnorm(1,mean(ij_b$length, na.rm = TRUE), 0.2),
         rlnorm(1, 0, 0.1),
         rnorm(1,mean(ij_b$BHDF, na.rm = TRUE), 0.2),
         rlnorm(1, 0, 0.1))
  
  sigma = rlnorm(J)
  
  L = diag(J)
  
  par = matrix(NA,n_ind, J)
  for(j in 1:J){
    if(j == 1 | j == 3){
      par[,j] = rnorm(n_ind, mu[j], 0.2)  
    } else {
      par[,j] = rlnorm(n_ind, log(mu[j]), 0.1)
    }
  } 
  
  sigma_obs = rlnorm(2)
  
  rho_obs = runif(1)

  return(list(t0p = t0p, 
              L = L, mu = mu, sigma = sigma, par = par, 
              sigma_obs = sigma_obs, rho_obs = rho_obs
              ))
}

# run stan ----

fit_vb <- stan_fit$sample(
  data = list(
    
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
    
    n_ind = n_ind,
    J = J
    
  ),
  init = init_vb,
  chains = 4,
  iter_warmup = 100,
  iter_sampling = 1000,
  thin = 1,
  save_warmup = FALSE,
  max_treedepth = 10,
  parallel_chains = 4,
  refresh = 100
)

 #saveRDS(fit_vb, file = paste0("fit_vb_",Sys.Date(),".rds"))
getwd()

my_dir <- 'C:/Users/leahm/OneDrive - University of Otago/Documents/git-otago/FBD_measurements/stan_output/'
fit_vb$save_output_files(basename = "fit_vb", random = FALSE)
fit_vb$save_data_file(dir = ".", basename = "fit_vb", timestamp = TRUE, random = FALSE)

# results -----

parout = as_draws_df(fit_vb$draws(c("mu","sigma","sigma_obs","t0p","rho_obs","corr[1,2]","corr[1,3]","corr[1,4]","corr[2,3]","corr[2,4]","corr[3,4]","Lobs[1,1]","Lobs[2,1]","Lobs[2,2]")))
saveRDS(parout, file = paste0("parout_",Sys.Date(),"_2.rds"))

#trace plot
library(ggplot2)
mcmc_trace(parout)+theme_bw()
ggplot2::ggsave("traceplot_allo.png", device = "png", dpi = 300, height = 200, width = 300, units = 'mm')
#density plot
bayesplot::mcmc_dens(parout)
#summary
as.data.frame(summary(parout))

#draws of par from posterior
parindout = as_draws_df(fit_vb$draws(c("par")))

#save par for individual plotting
saveRDS(parindout, file = paste0("parindout_",Sys.Date(),"_2.rds"))

# read in results ----
date = "2024-03-16"
parout_in = readRDS(file = paste0('parout_',date,'.rds'))
bayesplot::mcmc_dens(parout_in)
bayesplot::mcmc_trace(parout_in)+theme_bw()
summ_paroutin<-as.data.frame(summary(parout_in))

max(summ_paroutin$rhat)
min(summ_paroutin$ess_bulk)
min(summ_paroutin$ess_tail)
max(summ_paroutin$ess_bulk)
max(summ_paroutin$ess_tail)

as.data.frame(parout_in)
parindout_in = readRDS(file = paste0('./parindout_',date,'.rds'))
parindout_in_summ<-summary(parindout_in)

### if 
csv_files<-fit_vb$output_files()

fit_vb_in<-read_cmdstan_csv(
  csv_files,
  variables = c("mu")
  #sampler_diagnostics = NULL,
  #format = getOption("cmdstanr_draws_format", NULL)
)


# report results ----

## proportional relationship between z and y ----

x<-sample(1:n_ind,1, replace=F) 
x

### Ly given Lz ----

Lz = ij_b$BHDF[x]

rho_sigma_for_Ly = mean(parout$rho_obs)*(mean(parout$`sigma[1]`)/mean(parout$`sigma[3]`))

alpha_0_for_Ly = mean(parout$`mu[1]`) - rho_sigma_for_Ly*mean(parout$`mu[3]`)
alpha_1_for_Ly = rho_sigma_for_Ly*Lz

alpha_0_for_Ly+alpha_1_for_Ly

sd_for_Ly = (1-mean(parout$rho_obs)^2)*mean(parout$`sigma[1]`)^2 
sd_for_Ly

ij_b$length[x]

### plot a couple of individuals
induse = c(1, 12, 14, 50)
ngrid = 101
agegrid = seq(from = 0, to = max(ij_b$age), length.out = ngrid)
nind = length(induse)
pdf('indplots_bhdf.pdf', height = 8, width = 8)
par(mfrow = c(2,2), mar = c(4, 4, 1, 1))

for(i in 1:n_ind){
i = 28
  #x = ij_b%>%filter(ind == i)%>%select(age)
  #y = ij_b%>%filter(ind == i)%>%select(length)
  #z = ij_b%>%filter(ind == i)%>%select(BHDF)
  
  #plot(x$age, y$length, pch = 20, xlim = c(0,max(ij_b$age)), ylim = c(0, max(ij_b$length)), xlab = "Age", ylab = "Length", col = "blue")
  #points(x$age, z$BHDF, pch = 20, xlim = c(0,max(ij_b$age)), ylim = c(0, max(ij_b$length)), xlab = "Age", ylab = "Length", col = "red")

  Ly = parindout_in[[paste0('par[',i,',1]')]]
  ky = parindout_in[[paste0('par[',i,',2]')]]
  Lz = parindout_in[[paste0('par[',i,',3]')]]
  kz = parindout_in[[paste0('par[',i,',4]')]]
  
  tmp_y = matrix(NA,length(Ly),ngrid)
  tmp_z = matrix(NA,length(Lz),ngrid)
  
  for(j in 1:ngrid){
    #inverse logit kyout/kzout #exp(x)/(1+exp(x))
    tmp_y[,j] = Ly*(1-(exp(ky)/(1+exp(ky))^(mean(parout[["t0p"]]) + agegrid[j])))
    tmp_z[,j] = Lz*(1-(exp(kz)/(1+exp(kz))^(mean(parout[["t0p"]]) + agegrid[j])))
  }  
  
  quan_y = apply(tmp_y, 2, quantile, c(0.05, 0.5, 0.95), na.rm = T)
  #quan_z = apply(tmp_z, 2, quantile, c(0.05, 0.5, 0.95), na.rm = T)
  
  #matrix(NA, length(quan_y*n_ind), ngrid)
  
  #lines(agegrid, quan_y[2,], col = "blue", lty = 1)  
  #lines(agegrid, quan_y[1,], col = "blue", lty = 2)
  #lines(agegrid, quan_y[3,], col = "blue", lty = 2)
  
  #lines(agegrid, quan_z[2,], col = "red", lty = 1)  
  #lines(agegrid, quan_z[1,], col = "red", lty = 2)
  #lines(agegrid, quan_z[3,], col = "red", lty = 2)


  
  }


  

dev.off()

parindout
#### need to adjust for new model outputs above
mcmc_intervals(parindout_in[1:143], outer_size = 0.5, inner_size = 1, point_size = 2)
mcmc_intervals(parindout_in[144:286], outer_size = 0.5, inner_size = 1, point_size = 2)
mcmc_intervals(parindout_in[287:429], outer_size = 0.5, inner_size = 1, point_size = 2)
mcmc_intervals(parindout_in[430:ncol(parindout_in)], outer_size = 0.5, inner_size = 1, point_size = 2)

## L and k by birthyear----

parindout_in_summ$variable

id_parsumm<-parindout_in_summ%>%
  mutate(param = as.factor(
    case_when(
    grepl(",1]", variable) ~ 'Ly',
    grepl(",2]", variable) ~ 'ky',
    grepl(",3]", variable) ~ 'Lz',
    grepl(",4]", variable) ~ 'kz'
  )))%>%
  mutate(ind = as.numeric(stringr::str_extract(substr(variable, 5, nchar(variable)), '[^,]+')))%>%
  left_join(ij_ID, by = 'ind')
  

  
###
t0p = median(parout_in$t0p)
quantile(parout_in$t0p)
ngrid = 101

agegrid = seq(from = -2, to = max(ij_b$age), length.out = ngrid)


ind_median<-id_parsumm%>%
  group_by(ind)%>%
  tidyr::pivot_wider(names_from = "param", values_from = "median")%>%
  group_by(ind)%>%
  tidyr::fill(Ly,ky,Lz,kz, .direction = "downup")%>%
  mutate(ky = -1*log(exp(ky)/(1+exp(ky))),
         kz = -1*log(exp(kz)/(1+exp(kz))))#%>%
  distinct(ind, ID, year_zero, age_value, SEX, POD, Ly, ky, Lz, kz)



age_vb_y = matrix(NA,nrow(ind_median),ngrid)
age_vb_byr = matrix(NA,nrow(ind_median),ngrid)

for (i in 1:nrow(ind_median)){

  for(j in 1:ngrid){
    #inverse logit kyout/kzout #exp(x)/(1+exp(x))
    age_vb_y[i,j] = ind_median$Ly[i]*(1-exp(-ind_median$ky[i]*(t0p + agegrid[j])))
    age_vb_byr[i,j] = agegrid[j] + ind_median$year_zero[i]
    }  
}

by_df<-transform(expand.grid(i = seq(nrow(age_vb_byr)), j = seq(ncol(age_vb_byr))), year_by = c(age_vb_byr))
y_est<-transform(expand.grid(i = seq(nrow(age_vb_y)), j = seq(ncol(age_vb_y))), y_est = c(age_vb_y))

growest_plot<-by_df%>%
  left_join(y_est)%>%
  arrange(i,j)%>%
  mutate(agegrid = rep(agegrid,143))%>%
  left_join(ij_ID, by = c("i" = "ind"))%>%
  group_by(ID)%>%
  mutate(est_diff = lead(y_est)-y_est)%>%
  dplyr::rename("Pod" = "POD")


ggplot(growest_plot)+
  geom_histogram(mapping = aes(x = est_diff), binwidth = 0.001)

growest_plot%>%
  filter(est_diff <= 0.006)%>%
  group_by(ID)%>%
  mutate(rank = rank(j))%>%
  filter(rank == 1)%>%
  ungroup()%>%
  mutate(mean_limit = mean(agegrid))%>%
  distinct(mean_limit)

vbgc_zero<-ggplot(growest_plot)+
  geom_line(aes(x = agegrid, y = y_est, group = as.factor(i), color = Pod), alpha = 0.6)+
  coord_cartesian(xlim=c(0, 45))+
  coord_cartesian(ylim=c(0, 3.5))+
  #facet_wrap(~POD)+
  theme_bw()+
  xlab("Age (years)")+
  ylab("Total length estimate (m)")+
  theme(legend.position = "bottom")

#ggplot2::ggsave("./Figures/pod_vbplot.png", device = "png", dpi = 300, height = 150, width = 200, units = 'mm')

# ggplot(growest_plot)+
#   geom_line(aes(x = agegrid, y = y_est, group = as.factor(i)), alpha = 0.3)+
#   coord_cartesian(xlim=c(0, 45))+
#   coord_cartesian(ylim=c(0, 3.5))+
#   facet_wrap(~age_value)+
#   theme_bw()+
#   xlab("Age")+
#   ylab("Total length estimate (m)")

vbgc_by<-ggplot(growest_plot)+
  geom_line(aes(x = year_by, y = y_est, group = as.factor(i), color = Pod, linetype = age_value), alpha = 0.8)+
  xlim(c(1980,2050))+
  coord_cartesian(ylim=c(0, 3.5))+
  theme_bw()+
  xlab("Birth year")+
  ylab("")+
  theme(legend.position = "bottom",
        legend.title = element_blank())

ggpubr::ggarrange(vbgc_zero, vbgc_by, common.legend = T, legend = "bottom", widths = c(1,2), labels = "auto")

ggplot2::ggsave("./Figures/vbgcplot.png", device = "png", dpi = 300, height = 75, units = 'mm')

ky_box<-ggplot(ind_median)+
  geom_boxplot(aes(x = as.factor(SEX), y = ky, fill = SEX), alpha = 0.6)+
  theme_bw()+
  xlab("")+
  ylab("Estimated growth rate (K)")

maxLy_box<-ggplot(ind_median)+
  geom_boxplot(aes(x = as.factor(SEX), y = Ly, fill = SEX), alpha = 0.6)+
  theme_bw()+
  xlab("")+
  ylab("Estimated max length (m)")

box<-ggpubr::ggarrange(maxLy_box, ky_box, common.legend = T, legend = "none", ncol = 1)
box
ggplot2::ggsave("./Figures/box.png", device = "png", dpi = 300, height = 150, width = 75, units = 'mm')

ind_mean%>%
  group_by(POD)%>%
  tally()

ind_median%>%
  group_by(POD)%>%
  summarise(median = median(Ly), max = max(Ly), min = min(Ly))

ind_median%>%
  group_by(POD)%>%
  summarise(median = median(ky), max = max(ky), min = min(ky))

##
####need to transform ky and q5 and q95
id_parsumm%>%
  group_by(ind)%>%
  tidyr::pivot_wider(names_from = "param", values_from = "median")%>%
  group_by(ind)%>%
  tidyr::fill(Ly,ky,Lz,kz, .direction = "downup")%>%
  mutate(ky = -1*log(exp(ky)/(1+exp(ky))),
         kz = -1*log(exp(kz)/(1+exp(kz))))%>%
  distinct(ind, ID, year_zero, age_value, SEX, POD, Ly, ky, Lz, kz)

uncertainty<-id_parsumm%>%filter(param == "Ly" | param == "ky")%>%arrange(param, year_zero)%>%group_by(year_zero)%>%
  mutate(rank = rank(ind))%>%
  mutate(year_cat = as.numeric(year_zero+(rank*0.1)))%>%
  dplyr::rename(`Birth year` = "age_value")%>%
  mutate(median = case_when(
      param == "ky" ~ -1*log(exp(median)/(1+exp(median))),
      TRUE ~ median),
    q5 = case_when(
      param == "ky" ~ -1*log(exp(q5)/(1+exp(q5))),
      TRUE ~ q5),
    q95 = case_when(
      param == "ky" ~ -1*log(exp(q95)/(1+exp(q95))),
      TRUE ~ q95))
head(uncertainty)

ggplot(uncertainty%>%filter(param == "Ly"))+
  geom_point(aes(x = year_cat, y = median, group = as.factor(i), color = POD, shape = `Birth year`), size = 2, alpha = 0.6)+
  geom_linerange(aes(x = year_cat, ymin = q95, ymax = q5, group = as.factor(i), color = POD))+
  facet_wrap(~param, ncol = 1)+
  theme_bw()+
  theme(legend.position = "bottom")+
  xlab("")+
  ylab("Length (m)")%>%
  scale_x_continuous(n.breaks = 38)

uncertainty%>%filter(param == "ky")%>%filter(q5 == max(q5))

