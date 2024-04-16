library(cmdstanr)
library(posterior)
library(bayesplot)
library(latex2exp)
library(dplyr)
library(rstan)

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
# z matrix and z_t are the individual deviations from the population parameters

stan_fit = cmdstan_model("./scripts/stan/vb_mod_all_t0i.stan")

# initial values ----

init_vb = function(){
  
  mu_t = rlnorm(1)
  z_t = rnorm(n_ind)
  # t0p = rlnorm(n_ind, log(mu_t), 0.1)
  sigma_t = rlnorm(1, 0, 0.1)
  
  mu = c(rnorm(1, mean(ij_b$length, na.rm = TRUE), 0.2),
         rlnorm(1, 0, 0.1),
         rnorm(1, mean(ij_b$BHDF, na.rm = TRUE), 0.2),
         rlnorm(1, 0, 0.1))
  
  sigma = rlnorm(J)
  
  L = diag(J)
  
  # par = matrix(NA, n_ind, J)
  # 
  # for(j in 1:J){
  #   if(j == 1 | j == 3){
  #     par[,j] = rnorm(n_ind, mu[j], 0.2)  
  #   } else {
  #     par[,j] = rlnorm(n_ind, log(mu[j]), 0.1)
  #   }
  # } 
  
  par = matrix(rnorm(n_ind*J), n_ind, J)
  
  sigma_obs = rlnorm(2)
  
  rho_obs = runif(1)

  return(list(z_t = z_t, 
              mu_t = mu_t, sigma_t = sigma_t,
              L = L, mu = mu, sigma = sigma, z = par, 
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
  iter_warmup = 1000,
  iter_sampling = 10000,
  thin = 1,
  save_warmup = FALSE,
  max_treedepth = 10,
  parallel_chains = 4,
  refresh = 100
)

# results -----

parout = as_draws_df(fit_vb$draws(c("mu","sigma","sigma_obs","mu_t","sigma_t","rho_obs","corr[1,2]","corr[1,3]","corr[1,4]","corr[2,3]","corr[2,4]","corr[3,4]","Lobs[1,1]","Lobs[2,1]","Lobs[2,2]")))
saveRDS(parout, file = paste0("parout_",Sys.Date(),".rds"))

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
t0pindout = as_draws_df(fit_vb$draws(c("t0p")))
#save par for individual plotting
saveRDS(parindout, file = paste0("parindout_",Sys.Date(),".rds"))
saveRDS(t0pindout, file = paste0("t0pindout_",Sys.Date(),".rds"))
# read in results ----
date = "2024-04-16"
# remember k is logit(k)
#all non individual based params
parout_in = readRDS(file = paste0('parout_',date,'.rds'))
as.data.frame(parout_in)
bayesplot::mcmc_dens(parout_in)
bayesplot::mcmc_trace(parout_in)+theme_bw()
summ_paroutin<-as.data.frame(summary(parout_in))

mu12<-mcmc_scatter(
  as.matrix(parout_in),
  pars = c("mu[1]", "mu[2]"))

mu34<-mcmc_scatter(
  as.matrix(parout_in),
  pars = c("mu[3]", "mu[4]"))

mu13<-mcmc_scatter(
  as.matrix(parout_in),
  pars = c("mu[1]", "mu[3]"))

mu24<-mcmc_scatter(
  as.matrix(parout_in),
  pars = c("mu[2]", "mu[4]"))

max(summ_paroutin$rhat)
min(summ_paroutin$ess_bulk)
min(summ_paroutin$ess_tail)
max(summ_paroutin$ess_bulk)
max(summ_paroutin$ess_tail)

#individual based params
parindout_in = readRDS(file = paste0('./parindout_',date,'.rds'))
parindout_in_summ<-summary(parindout_in)

t0pindout_in = readRDS(file = paste0('./t0pindout_',date,'.rds'))
t0pindout_in_summ<-summary(t0pindout_in)
# report results ----

### plot a couple of individuals ----
induse = c(1, 12, 99, 90)
ngrid = 101
agegrid = seq(from = -2, to = max(ij_b$age), length.out = ngrid)
nind = length(induse)
pdf('indplots_bhdf.pdf', height = 8, width = 8)
par(mfrow = c(2,2), mar = c(4, 4, 1, 1))

for(i in 1:nind){

  i = induse[i]
  print(i)
  x = ij_b%>%filter(ind == i)%>%select(age)
  y = ij_b%>%filter(ind == i)%>%select(length)
  z = ij_b%>%filter(ind == i)%>%select(BHDF)
  
  plot(x$age, y$length, pch = 20, xlim = c(0,max(ij_b$age)), ylim = c(0, max(ij_b$length)), xlab = "Age", ylab = "Length", col = "blue")
  points(x$age, z$BHDF, pch = 20, xlim = c(0,max(ij_b$age)), ylim = c(0, max(ij_b$length)), xlab = "Age", ylab = "Length", col = "red")

  Ly = parindout_in[[paste0('par[',i,',1]')]]
  ky = parindout_in[[paste0('par[',i,',2]')]]
  Lz = parindout_in[[paste0('par[',i,',3]')]]
  kz = parindout_in[[paste0('par[',i,',4]')]]
  t0p = t0pindout_in[[paste0('t0p[',i,']')]]
  
  tmp_y = matrix(NA,length(Ly),ngrid)
  tmp_z = matrix(NA,length(Lz),ngrid)
  
  for(j in 1:ngrid){
    #inverse logit kyout/kzout 1/(1+exp(-k))
    tmp_y[,j] = Ly*(1-(1/(1+exp(-ky))^(t0p + agegrid[j])))
    tmp_z[,j] = Lz*(1-(1/(1+exp(-kz))^(t0p + agegrid[j])))
  }  
  
  quan_y = apply(tmp_y, 2, quantile, c(0.05, 0.5, 0.95), na.rm = T)
  quan_z = apply(tmp_z, 2, quantile, c(0.05, 0.5, 0.95), na.rm = T)
  
  matrix(NA, length(quan_y*n_ind), ngrid)
  
  lines(agegrid, quan_y[2,], col = "blue", lty = 1)  
  lines(agegrid, quan_y[1,], col = "blue", lty = 2)
  lines(agegrid, quan_y[3,], col = "blue", lty = 2)
  
  lines(agegrid, quan_z[2,], col = "red", lty = 1)  
  lines(agegrid, quan_z[1,], col = "red", lty = 2)
  lines(agegrid, quan_z[3,], col = "red", lty = 2)


  
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

id_t0psumm<-t0pindout_in_summ%>%
  mutate(param = "t0p")%>%
  mutate(ind = as.numeric(stringr::str_extract(substr(variable, 5, nchar(variable)), '[^]]+')))%>%
  left_join(ij_ID, by = 'ind')

#rename logit_ks
id_parsumm<-parindout_in_summ%>%
  mutate(param = as.factor(
    case_when(
    grepl(",1]", variable) ~ 'Ly',
    grepl(",2]", variable) ~ 'logit_ky',
    grepl(",3]", variable) ~ 'Lz',
    grepl(",4]", variable) ~ 'logit_kz'
  )))%>%
  mutate(ind = as.numeric(stringr::str_extract(substr(variable, 5, nchar(variable)), '[^,]+')))%>%
  left_join(ij_ID, by = 'ind')%>%
  bind_rows(id_t0psumm)
  
###
#t0p = median(parout_in$t0p)
#quantile(parout_in$t0p)
ngrid = 101

agegrid = seq(from = -2, to = max(ij_b$age), length.out = ngrid)

#inverse logit ks to just get k on 0-1 scale, k = e^-K
ind_median<-id_parsumm%>%
  group_by(ind)%>%
  tidyr::pivot_wider(names_from = "param", values_from = "median")%>%
  group_by(ind)%>%
  tidyr::fill(Ly,logit_ky,Lz,logit_kz,t0p, .direction = "downup")%>%
  #inverse logit kyout/kzout #exp(x)/(1+exp(x)), k = exp(-K)
  mutate(ky = 1/(1+exp(-logit_ky)),
         kz = 1/(1+exp(-logit_kz)))%>%
  distinct(ind, ID, year_zero, age_value, SEX, POD, Ly, logit_ky, ky, Lz, logit_kz, kz, t0p)

summary(ind_median)

age_vb_y = matrix(NA,nrow(ind_median),ngrid)
age_vb_byr = matrix(NA,nrow(ind_median),ngrid)

for (i in 1:nrow(ind_median)){
  
  for(j in 1:ngrid){
    age_vb_y[i,j] = ind_median$Ly[i]*(1-ind_median$ky[i]^(ind_median$t0p[i] + agegrid[j]))
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
  filter(est_diff <= 0.01)%>%
  group_by(ID)%>%
  mutate(rank = rank(j))%>%
  filter(rank == 1)%>%
  ungroup()%>%
  mutate(mean_limit = mean(agegrid))%>%
  distinct(mean_limit)

vbgc_zero<-ggplot(growest_plot)+
  geom_line(aes(x = agegrid, y = y_est, group = as.factor(i), color = Pod), alpha = 0.6)+
  coord_cartesian(xlim=c(0, 41), ylim=c(0, 3.5))+
  #facet_wrap(~SEX)+
  theme_bw()+
  xlab("Age (years)")+
  ylab("Total length estimate (m)")+
  theme(legend.position = "bottom")

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
  xlim(c(1980,2065))+
  coord_cartesian(ylim=c(0, 3.5))+
  theme_bw()+
  xlab("Birth year")+
  ylab("")+
  theme(legend.position = "bottom",
        legend.title = element_blank())

ggpubr::ggarrange(vbgc_zero, vbgc_by, common.legend = T, legend = "bottom", widths = c(1,2), labels = "auto")

ggplot2::ggsave(paste0("./Figures/vbgcplot_",Sys.Date(),".png"), device = "png", dpi = 700, width = 250, height = 125, units = 'mm')

ind_median%>%
  group_by(POD)%>%
  summarise(median = median(Ly), max = max(Ly), min = min(Ly))

#k = e^-K
ind_median%>%
  group_by(POD)%>%
  summarise(median = median(ky), max = max(ky), min = min(ky))

##
####need to transform ky and q5 and q95
wide<-id_parsumm%>%
  group_by(ind)%>%
  tidyr::pivot_wider(names_from = "param", values_from = "median")%>%
  group_by(ind)%>%
  tidyr::fill(Ly,logit_ky,Lz,logit_kz, .direction = "downup")%>%
  #inverse logit + log to get K out
  mutate(Ky = -1*log(1/(1+exp(-logit_ky))), #K
         ky =  1/(1+exp(-logit_ky)), # k = e^-K
         Kz = -1*log(1/(1+exp(-logit_kz))), #K
         kz =  1/(1+exp(-logit_kz)))%>% # k = e^-K
  distinct(ind, ID, year_zero, age_value, SEX, POD, Ly, logit_ky, Ky, ky, Lz, logit_kz, Kz, kz)

summary(wide)

uncertainty<-id_parsumm%>%filter(param == "Ly" | param == "logit_ky")%>%ungroup()%>%arrange(param, year_zero)%>%group_by(year_zero)%>%
  mutate(rank = percent_rank(ind))%>%
  mutate(year_cat = as.numeric(year_zero+(rank)))%>%
  dplyr::rename(`Birth year` = "age_value")#%>%
  # mutate(median = case_when(
  #     param == "ky" ~ -1*log(exp(median)/(1+exp(median))),
  #     TRUE ~ median),
  #   q5 = case_when(
  #     param == "ky" ~ -1*log(exp(q5)/(1+exp(q5))),
  #     TRUE ~ q5),
  #   q95 = case_when(
  #     param == "ky" ~ -1*log(exp(q95)/(1+exp(q95))),
  #     TRUE ~ q95))

head(uncertainty)

ggplot(uncertainty%>%filter(param == "Ly"))+
  geom_point(aes(x = year_cat, y = median, group = as.factor(i), color = POD, shape = `Birth year`), size = 2, alpha = 0.6)+
  geom_linerange(aes(x = year_cat, ymin = q95, ymax = q5, group = as.factor(i), color = POD))+
  facet_wrap(~param, ncol = 1)+
  theme_bw()+
  theme(legend.position = "bottom")+
  xlab("")+
  ylab("Length (m)")+
  scale_x_continuous(breaks = seq(1984, 2024, 1))+
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -90, vjust = -0.5))

ggplot2::ggsave(paste0("./Figures/Ly_CI_", Sys.Date() ,".png"), device = "png", dpi = 700, height = 100, width = 300, units = 'mm')

ggplot(uncertainty%>%filter(param == "logit_ky"))+
  geom_point(aes(x = year_cat, y = median, group = as.factor(i), color = POD, shape = `Birth year`), size = 2, alpha = 0.6)+
  geom_linerange(aes(x = year_cat, ymin = q95, ymax = q5, group = as.factor(i), color = POD))+
  facet_wrap(~param, ncol = 1)+
  theme_bw()+
  theme(legend.position = "bottom")+
  xlab("")+
  ylab("Growth rate")+
  scale_x_continuous(breaks = seq(1984, 2024, 1))+
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -90, vjust = -0.5))

ggplot2::ggsave(paste0("./Figures/ky_CI_", Sys.Date() ,".png"), device = "png", dpi = 700, height = 100, width = 300, units = 'mm')

ggplot(wide)+
  geom_point(aes(x = logit_ky, y = Ly, color = "logit(ky) ~ Ly", shape = age_value))+
  #geom_point(aes(x = logit_kz, y = Lz, color = "logit(kz) ~ Lz", shape = age_value))+
  geom_point(aes(x = logit_kz, y = Ly, color = "logit(kz) ~ Ly", shape = age_value))
  #geom_point(aes(x = logit_ky, y = Lz, color = "logit(ky) ~ Lz"))

ggplot(wide, aes(x = -Ky, y = Ly, color = SEX))+
  geom_point()

ggplot(id_parsumm%>%
         filter(param == "logit_kz" | param == "Ly")%>%
         tidyr::pivot_wider(names_from = "param", values_from = "median"))+
  geom_point(mapping = aes(x = logit_kz, y = Ly, color = age_value))

ggplot(wide, aes(x = -Kz, y = Ly, color = age_value))+
  geom_point()

ggplot(ind_median, aes(x = kz, y = Lz, color = age_value))+
  geom_point()

ggplot(ind_median, aes(x = Lz, y = Ly, color = age_value))+
  geom_point()

ggplot(ind_median, aes(x = ky, y = kz, color = age_value))+
  geom_point()


## proportional relationship between z and y ----

x<-sample(1:n_ind,1, replace=F) 
x

### Ly given Lz ----

(median(parout$`sigma[1]`)/median(parout$`sigma[3]`))


z_grid = seq(from = min(ij_b$BHDF), max(ij_b$BHDF), length.out = 101)
y_est = NULL

(median(parout$`sigma[3]`)/median(parout$`sigma[1]`))

(1 - median(parout$rho_obs)^2) * median(parout$`sigma[1]`)^2
rho_sigma_for_Ly = median(parout$rho_obs)*(median(parout$`sigma[1]`)/median(parout$`sigma[3]`))
alpha_0_for_Ly = median(parout$`mu[1]`) - rho_sigma_for_Ly*median(parout$`mu[3]`)

for (i in 1:length(z_grid)){
  i = 1
  #alpha_1_for_Ly = rho_sigma_for_Ly*z_grid[i]
  
  #y_est[i]<-alpha_0_for_Ly+alpha_1_for_Ly
  y_est[i] = median(parout$`mu[1]`)+((median(parout$`sigma[1]`)/median(parout$`sigma[3]`))*median(parout$rho_obs)*(z_grid[i] - median(parout$`mu[3]`)))
  
}

y_est

sd_for_Ly = (1-median(parout$rho_obs)^2)*median(parout$`sigma[1]`)^2 
sd_for_Ly

lm(ij_b$BHDF ~ ij_b$length)

ggplot(ij_b, aes(x=BHDF, y = length))+
  geom_point()+ 
  geom_smooth(method='lm', formula= y~x)+
  #ylim(c(0,3.5))+
  #xlim(c(0,1.5))+
  geom_abline(intercept = 1.86, slope = 1.13, color = "red")
#geom_point(mapping = aes(x = z_grid, y = y_est), color = "blue")

ggplot()+
  geom_point(aes(x = parout$`mu[3]`, y = parout$`mu[1]`))

ggplot()+
  geom_point(aes(x = parout$`mu[4]`, y = parout$`mu[2]`))
