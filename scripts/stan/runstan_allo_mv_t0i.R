library(cmdstanr)
library(posterior)
library(bayesplot)
library(latex2exp)
library(dplyr)

# data ----
i = 34

#both bhdf & tl
ij_b = readRDS(file = 'ij_1.rds')
ij_z = readRDS(file = 'ij_2.rds')
ij_y = readRDS(file = 'ij_3.rds')

ij_all<-ij_b%>%
  bind_rows(ij_z)%>%
  bind_rows(ij_y)

ij_all_old<-ij_b_old%>%
  bind_rows(ij_z_old)%>%
  bind_rows(ij_y_old)

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

#IDs
ij_ID = readRDS(file = 'ij_ID.rds')

n_ind<-nrow(ij_ID)

N_b+N_z+N_y
N_b/(N_b+N_z+N_y)

#params matrix dims
J=4

# model ----
# z matrix and z_t are the individual deviations from the population parameters

set_cmdstan_path(path = "C:/Users/leahm/cmdstan-2.34.1")

stan_fit = cmdstan_model("./scripts/stan/vb_mod_all_t0i.stan")

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
  
  # par = matrix(NA, n_ind, J)
  # 
  # for(j in 1:J){
  #   if(j == 1 | j == 3){
  #     par[,j] = rnorm(n_ind, mu[j], 0.2)  
  #   } else {
  #     par[,j] = rlnorm(n_ind, log(mu[j]), 0.1)
  #   }
  # } 
  
  z = matrix(rnorm(n_ind*J), n_ind, J)
  
  sigma_obs = rlnorm(2)
  
  rho_obs = runif(1)

  return(list(
              z_t = z_t, mu_t = mu_t, t0p = t0p, sigma_t = sigma_t,
              L = L, mu = mu, sigma = sigma, z = z, 
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
  iter_warmup = 1500,
  iter_sampling = 15000,
  thin = 1,
  save_warmup = FALSE,
  max_treedepth = 10,
  parallel_chains = 4,
  refresh = 100
)

# results -----

parout = as_draws_df(fit_vb$draws(c("mu","mu_t","sigma","sigma_t","sigma_obs","rho_obs",
                                    "varcov_par",
                                    "corr[1,2]","corr[1,3]","corr[1,4]",
                                    "corr[2,3]","corr[2,4]",
                                    "corr[3,4]",
                                    "Lobs[1,1]","Lobs[2,1]","Lobs[2,2]")))
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
date = "2024-04-30"
# remember k is logit(k)
#all non individual based params
parout_in = readRDS(file = paste0('parout_',date,'.rds'))
as.data.frame(parout_in)
bayesplot::mcmc_dens(parout_in)
bayesplot::mcmc_trace(parout_in)+theme_bw()
summ_paroutin<-as.data.frame(summary(parout_in))

max(summ_paroutin$rhat)
min(summ_paroutin$ess_bulk)
min(summ_paroutin$ess_tail)
max(summ_paroutin$ess_bulk)
max(summ_paroutin$ess_tail)

max(summary(parout)$rhat)
min(summary(parout)$ess_bulk)
min(summary(parout)$ess_tail)
#individual based params
parindout_in = readRDS(file = paste0('./parindout_',date,'.rds'))
parindout_in_summ<-summary(parindout_in)

t0pindout_in = readRDS(file = paste0('./t0pindout_',date,'.rds'))
t0pindout_in_summ<-summary(t0pindout_in)
# report results ----

### plot a couple of individuals ----
induse = c(42,95,23)
#induse = ceiling(runif(4, 0, 143))

ngrid = 101
agegrid = seq(from = -2, to = max(ij_b$age), length.out = ngrid)
nind = length(induse)
#pdf('indplots_bhdf.pdf', height = 8, width = 8)
#par(mfrow = c(2,2), mar = c(4, 4, 1, 1))
label = c("a","b","c")

for(i in 1:nind){
  k = label[i]
  i = induse[i]
  print(i)
  x = ij_all%>%filter(ind == i)%>%select(age)
  y = ij_all%>%filter(ind == i)%>%select(length)
  z = ij_all%>%filter(ind == i)%>%select(BHDF)
  
  #plot(x$age, y$length, pch = 20, xlim = c(0,max(ij_b$age)), ylim = c(0, max(ij_b$length)), xlab = "Age", ylab = "Length", col = "blue")
  #points(x$age, z$BHDF, pch = 20, xlim = c(0,max(ij_b$age)), ylim = c(0, max(ij_b$length)), xlab = "Age", ylab = "Length", col = "red")

  Ly = parindout_in[[paste0('par[',i,',1]')]]
  ky = parindout_in[[paste0('par[',i,',2]')]]
  Lz = parindout_in[[paste0('par[',i,',3]')]]
  kz = parindout_in[[paste0('par[',i,',4]')]]
  t0p = t0pindout_in[[paste0('t0p[',i,']')]]
  
  # mu_1 = parout[['mu[1]']]
  # mu_2 = parout[['mu[2]']]
  # mu_3 = parout[['mu[3]']]
  # mu_4 = parout[['mu[4]']]
  # mu_t = parout[['mu_t']]
  
  tmp_y = matrix(NA,length(Ly),ngrid)
  tmp_z = matrix(NA,length(Lz),ngrid)
  
  for(j in 1:ngrid){
    #inverse logit kyout/kzout 1/(1+exp(-k))
    tmp_y[,j] = Ly*(1-(1/(1+exp(-ky))^(t0p + agegrid[j])))
    tmp_z[,j] = Lz*(1-(1/(1+exp(-kz))^(t0p + agegrid[j])))
  }  
  
  # for(j in 1:ngrid){
  #   #inverse logit kyout/kzout 1/(1+exp(-k))
  #   tmp_y[,j] = mu_1*(1-(1/(1+exp(-mu_2))^(mu_t + agegrid[j])))
  #   tmp_z[,j] = mu_3*(1-(1/(1+exp(-mu_4))^(mu_t + agegrid[j])))
  # } 
  
  quan_y = apply(tmp_y, 2, quantile, c(0.05, 0.5, 0.95), na.rm = T)
  quan_z = apply(tmp_z, 2, quantile, c(0.05, 0.5, 0.95), na.rm = T)
  
  matrix(NA, length(quan_y*n_ind), ngrid)
  
  #lines(agegrid, quan_y[2,], col = "blue", lty = 1)  
  #lines(agegrid, quan_y[1,], col = "blue", lty = 2)
  #lines(agegrid, quan_y[3,], col = "blue", lty = 2)
  
  #lines(agegrid, quan_z[2,], col = "red", lty = 1)  
  #lines(agegrid, quan_z[1,], col = "red", lty = 2)
  #lines(agegrid, quan_z[3,], col = "red", lty = 2)

  ggplot()+
      geom_point(aes(x = x$age, y = y$length, color = "TL"), size = 2.5, alpha = 0.8)+
      geom_point(aes(x = x$age, y = z$BHDF, color = "BHDF"), size = 2.5, alpha = 0.8)+
      geom_line(aes(x = agegrid, y = quan_y[1,], color = "TL"), linetype = "dashed")+
      geom_line(aes(x = agegrid, y = quan_y[2,], color = "TL"))+
      geom_line(aes(x = agegrid, y = quan_y[3,], color = "TL"), linetype = "dashed")+
      geom_line(aes(x = agegrid, y = quan_z[1,], color = "BHDF"), linetype = "dashed")+
      geom_line(aes(x = agegrid, y = quan_z[2,], color = "BHDF"))+
      geom_line(aes(x = agegrid, y = quan_z[3,], color = "BHDF"), linetype = "dashed")+
      coord_cartesian(ylim=c(0, 3.5))+
      theme_bw()+
      ylab("Length (m)")+
      xlab("Age (years)")+
      theme(legend.position = "none")+
      scale_color_viridis_d(begin = 0.1, end = 0.8)+
      labs(title = k)
  
  ggplot2::ggsave(paste0("./Figures/indplot_",i,".png"), device = "png", dpi = 700, width = 100, height = 100, units = 'mm')
  
  }

#dev.off()

parindout

#### need to adjust for new model outputs above
mcmc_intervals(parindout_in[1:143], outer_size = 0.5, inner_size = 1, point_size = 2)
mcmc_intervals(parindout_in[144:286], outer_size = 0.5, inner_size = 1, point_size = 2)
mcmc_intervals(parindout_in[287:429], outer_size = 0.5, inner_size = 1, point_size = 2)
mcmc_intervals(parindout_in[430:ncol(parindout_in)], outer_size = 0.5, inner_size = 1, point_size = 2)

## L and k by birthyear----

#parindout_in_summ$variable

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
  geom_line(aes(x = agegrid, y = y_est, group = as.factor(i), color = Pod), alpha = 0.6, linewidth = 0.2)+
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
  geom_line(aes(x = year_by, y = y_est, group = as.factor(i), color = Pod, linetype = age_value), alpha = 0.8, linewidth = 0.2)+
  xlim(c(1980,2065))+
  coord_cartesian(ylim=c(0, 3.5))+
  theme_bw()+
  xlab("Birth year")+
  ylab("")+
  theme(legend.position = "bottom",
        legend.title = element_blank())

vbgc_by+
  facet_wrap(~age_value)

growest_plot%>%
  group_by(age_value)%>%
  filter(y_est == max(y_est))

ggpubr::ggarrange(vbgc_zero, vbgc_by, common.legend = T, legend = "bottom", widths = c(1,2), labels = "auto")

ggplot2::ggsave(paste0("./Figures/vbgcplot.png"), device = "png", dpi = 700, width = 250, height = 100, units = 'mm')

ind_median%>%
  group_by(POD)%>%
  summarise(median = median(Ly), max = max(Ly), min = min(Ly))

#k = e^-K
ind_median%>%
  group_by(POD)%>%
  summarise(median = median(ky), max = max(ky), min = min(ky))


ind_median%>%
  group_by(age_value)%>%
  summarise(median = median(Ly), max = max(Ly), min = min(Ly))

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
  geom_point(aes(x = year_cat, y = median, group = as.factor(i), color = POD, shape = `Birth year`), size = 1.5, alpha = 0.6)+
  geom_linerange(aes(x = year_cat, ymin = q95, ymax = q5, group = as.factor(i), color = POD), linewidth = 0.2)+
  facet_wrap(~param, ncol = 1)+
  theme_bw()+
  theme(legend.position = "bottom")+
  xlab("")+
  ylab("Length (m)")+
  scale_x_continuous(breaks = seq(1984, 2024, 1))+
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -90, vjust = -0.5))

ggplot2::ggsave(paste0("./Figures/Ly_CI.png"), device = "png", dpi = 700, height = 100, width = 300, units = 'mm')

ggplot(uncertainty%>%filter(param == "logit_ky")%>%mutate(param = "ky"))+
  geom_point(aes(x = year_cat, y = 1/(1+exp(-median)), group = as.factor(i), color = POD, shape = `Birth year`), size = 1.5, alpha = 0.6)+
  geom_linerange(aes(x = year_cat, ymin = 1/(1+exp(-q95)), ymax = 1/(1+exp(-q5)), group = as.factor(i), color = POD), linewidth = 0.2)+
  facet_wrap(~param, ncol = 1)+
  theme_bw()+
  theme(legend.position = "bottom")+
  xlab("")+
  ylab("Growth rate")+
  scale_x_continuous(breaks = seq(1984, 2024, 1))+
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -90, vjust = -0.5))

ggplot2::ggsave(paste0("./Figures/ky_CI.png"), device = "png", dpi = 700, height = 100, width = 300, units = 'mm')


# box plots ----
## sex ----

sex_box<-ggplot(id_parsumm%>%filter(param != "t0p")%>%mutate(median = case_when(
  param == "logit_ky" ~ 1/(1+exp(-median)),
  param == "logit_kz" ~ 1/(1+exp(-median)),
  TRUE ~ median
), 
param = case_when(
  param == "logit_ky" ~ "ky",
  param == "logit_kz" ~ "kz",
  TRUE ~ param
),
group = case_when(
  param == "ky" | param == "kz" ~ "Growth parameters",
  TRUE ~ param
)))+
  geom_boxplot(aes(x = param, y = median, fill = SEX), alpha = 0.6)+
  facet_wrap(~group, scales = "free")+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(legend.position = "bottom")+
  scale_fill_viridis_d()

## pod ----
pod_box<-ggplot(id_parsumm%>%filter(param != "t0p")%>%mutate(median = case_when(
  param == "logit_ky" ~ 1/(1+exp(-median)),
  param == "logit_kz" ~ 1/(1+exp(-median)),
  TRUE ~ median
), 
 param = case_when(
   param == "logit_ky" ~ "ky",
   param == "logit_kz" ~ "kz",
   TRUE ~ param
 ),
group = case_when(
  param == "ky" | param == "kz" ~ "Growth parameters",
  TRUE ~ param
)))+
  geom_boxplot(aes(x = param, y = median, fill = POD), alpha = 0.6)+
  facet_wrap(~group, scales = "free")+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(legend.position = "bottom")

ggpubr::ggarrange(sex_box, pod_box, ncol = 1, labels = "auto")

ggplot2::ggsave(paste0("./Figures/boxplots.png"), device = "png", dpi = 700, height = 200, width = 200, units = 'mm')

### Ly given Lz ----

muLy<-summ_paroutin%>%
  filter(variable == "mu[1]")

muky<-summ_paroutin%>%
  filter(variable == "mu[2]")

muLz<-summ_paroutin%>%
  filter(variable == "mu[3]")

mukz<-summ_paroutin%>%
  filter(variable == "mu[4]")

sigma1<-summ_paroutin%>%
  filter(variable == "sigma[1]")

sigma3<-summ_paroutin%>%
  filter(variable == "sigma[3]")

sigma_obs1<-summ_paroutin%>%
  filter(variable == "sigma_obs[1]")

sigma_obs2<-summ_paroutin%>%
  filter(variable == "sigma_obs[2]")

rho_obs<-summ_paroutin%>%
  filter(variable == "rho_obs")

corr13<-summ_paroutin%>%
  filter(variable == "corr[1,3]")

vc11<-summ_paroutin%>%
  filter(variable == "Lobs[1,1]")

vc12<-summ_paroutin%>%
  filter(variable == "Lobs[2,1]")

vc22<-summ_paroutin%>%
  filter(variable == "Lobs[2,2]")

vc11<-summ_paroutin%>%
  filter(variable == "varcov_par[1,1]")

vc12<-summ_paroutin%>%
  filter(variable == "varcov_par[1,3]")

vc22<-summ_paroutin%>%
  filter(variable == "varcov_par[3,3]")

varcov = matrix(NA, 2, 2)

varcov[1,1] = vc11$median
sigma1$median^2

varcov[2,2] = vc22$median
sigma3$median^2

varcov[1,2] = vc12$median
corr13$median*sigma1$median*sigma3$median

varcov[2,1] = varcov[1,2]

#regression coefficient / slope
reg_co<-(varcov[1,2]/varcov[2,2])
slope = reg_co
slope

y_intercept<-muLy$median - reg_co*muLz$median
y_intercept
###
#conditional variance
con_var<-varcov[1,1]-slope*varcov[2,1]

ij_est<-ij_all%>%
  mutate(Length_est = slope*BHDF + y_intercept)%>%
  mutate(`Age bin` = case_when(
    age < 10 ~ "<10",
    age >= 10 & age < 20 ~ "10-20",
    age >= 20 & age < 30 ~ "20-30",
    age >= 30 & age < 40 ~ "30-40",
    age >= 40 ~ "40+"
  ))

summary(ij_est)

ij_est

#look at Ly and Lz to see if this equation works
Lyest<-id_parsumm%>%
  filter(param == "Lz" | param == "Ly")%>%
  distinct(param, ind, ID, mean, median, sd, mad)%>%
  mutate(Length_est = case_when(
    param == "Lz" ~ slope*median + y_intercept,
    TRUE ~ NA))

x = Lyest%>%filter(param == "Lz")%>%dplyr::select(ind, median)
y = Lyest%>%filter(param == "Ly")%>%dplyr::select(ind, median)

join<-x%>%
  left_join(y, by = "ind")

linear_plot<-ggplot()+
  geom_point(ij_est, mapping = aes(x = BHDF, y = length, color = `Age bin`), alpha = 0.6, size = 2)+
  #geom_point(ij_est, mapping = aes(x = BHDF, y = Length_est))+
  geom_point(join, mapping = aes(x = median.x, y = median.y), alpha = 0.6, size = 2)+
  geom_abline(mapping = aes(slope = slope, intercept = y_intercept), size = 1.5)+
  #geom_abline(mapping = aes(slope = 2.408, intercept = 0.72521, color = "Vivier et al. 2023"), size = 1.5)+
  #geom_abline(mapping = aes(slope = 3.1314, intercept = 0.070626, color = "Cheney et al. 2018"), size = 1.5)+
  #geom_abline(mapping = aes(slope = 3.17, intercept = 0.050583, color = "van Aswegen et al. 2019"), size = 1.5)+
  scale_color_viridis_d(option="plasma", end = 0, begin = 0.88)+
  theme_bw()+
  theme(legend.position = "bottom")+
  ylab("TL (m)")+
  xlab("BHDF (m)")+
  xlim(c(0.4, 1.2))+
  ylim(c(1.3, 3.4))
  
fbd_allo<-ggMarginal(linear_plot)

lm<-lm(ij_est$length ~ ij_est$BHDF)
lm_intercept = lm$coefficients[[1]]
lm_slope = lm$coefficients[[2]]

compare_plot<-ggplot()+
  geom_point(ij_est, mapping = aes(x = BHDF, y = length),alpha = 0.1, size = 2)+
  #geom_point(ij_est, mapping = aes(x = BHDF, y = Length_est))+
  #geom_point(join, mapping = aes(x = median.x, y = median.y), alpha = 0.6, size = 2)+
  geom_abline(mapping = aes(slope = slope, intercept = y_intercept), size = 2)+
  geom_abline(mapping = aes(slope = 2.408, intercept = 0.72521, color = "Vivier et al. 2023"), size = 1.5, alpha = 0.6)+
  geom_abline(mapping = aes(slope = 3.1314, intercept = 0.070626, color = "Cheney et al. 2018"), size = 1.5, alpha = 0.8)+
  geom_abline(mapping = aes(slope = 3.17, intercept = 0.050583, color = "van Aswegen et al. 2019"), size = 1.5, alpha = 0.8)+
  scale_color_viridis_d(end = 0.4, begin = 0.99)+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_blank())+
  ylab("TL (m)")+
  xlab("BHDF (m)")+
  xlim(c(0.4, 1.2))+
  ylim(c(1.3, 3.4))

compare_plot

alloplot<-ggpubr::ggarrange(linear_plot, compare_plot, labels = "auto")

ggplot2::ggsave(paste0("./Figures/alloplot.png"), device = "png", dpi = 700, height = 125, width = 250, units = 'mm')
