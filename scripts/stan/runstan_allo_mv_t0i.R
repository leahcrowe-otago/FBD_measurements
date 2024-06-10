library(cmdstanr)
library(posterior)
library(bayesplot)
library(latex2exp)
library(dplyr)

# data ----
i = 34
#IDs
ij_ID = readRDS(file = './data/Measurements/Data for review/ij_ID.rds')
n_ind = nrow(ij_ID)

### reg model ----
#both bhdf & tl
ij_b = readRDS(file = './data/Measurements/Data for review/ij_1.rds')
ij_z = readRDS(file = './data/Measurements/Data for review/ij_2.rds')
ij_y = readRDS(file = './data/Measurements/Data for review/ij_3.rds')

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

### sex/pod model ----
obs_sex = readRDS(file = './data/Measurements/Data for review/obs_sex.rds')
n_k<-obs_sex%>%filter(SEX != "X")%>%distinct(ID)%>%nrow()
n_obs_k<-obs_sex%>%filter(SEX != "X")%>%nrow()
n_obs = nrow(obs_sex)
sex = ij_ID%>%
  dplyr::select(SEX)%>%
  filter(SEX != "X")%>%
  mutate(SEX = case_when(
    SEX == "F" ~ 0,
    SEX == "M" ~ 1))
pod = ij_ID%>%
  dplyr::select(POD)%>%
  mutate(POD = case_when(
    POD == "DOUBTFUL" ~ 0,
    POD == "DUSKY" ~ 1))
#params matrix dims
J=4

# model ----
# z matrix and z_t are the individual deviations from the population parameters

set_cmdstan_path(path = "C:/Users/leahm/cmdstan-2.34.1")

stan_fit = cmdstan_model("./scripts/stan/vb_mod_all_t0i.stan") #ms model
#stan_fit = cmdstan_model("./scripts/stan/vb_mod_all_t0i_sex.stan") #sex/pod effects model

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
  z = matrix(rnorm(n_ind*J), n_ind, J)
  sigma_obs = rlnorm(2)
  rho_obs = runif(1)
  
  #sex/pod
  beta = rlnorm(J)
  gamma = rlnorm(J)
  pi = runif(1)
    
  return(list(#obs/sex
              beta = beta, gamma = gamma, pi = pi,
              #both versions
              z_t = z_t, mu_t = mu_t, t0p = t0p, sigma_t = sigma_t,
              L = L, mu = mu, sigma = sigma, z = z, 
              sigma_obs = sigma_obs, rho_obs = rho_obs
              ))
}

# run stan ----

fit_vb <- stan_fit$sample(
  data = list(
    
  ### both models ----
    n_ind = n_ind,
    J = J,
  ### reg model ----
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
  ### sex/pod model ----
    n_k = n_k,
    n_obs_k = n_obs_k,
    n_obs = n_obs,
    id = obs_sex$ind,
    age = obs_sex$age,
    type = obs_sex$type,
    val = obs_sex%>%dplyr::select(length, BHDF)%>%replace(is.na(.), 0),
    sex = sex$SEX,
    pod = pod$POD
  
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

parout = as_draws_df(fit_vb$draws(c(#sex/pod model
                                    #"beta","gamma","pi",
                                    #both models
                                    "mu","mu_t","sigma","sigma_t","sigma_obs","rho_obs",
                                    "varcov_par",
                                    "corr[1,2]","corr[1,3]","corr[1,4]",
                                    "corr[2,3]","corr[2,4]",
                                    "corr[3,4]",
                                    "Lobs[1,1]","Lobs[2,1]","Lobs[2,2]",
                                    "mu_pred")))
#saveRDS(parout, file = paste0("./results/parout_",Sys.Date(),"sexpod.rds"))
saveRDS(parout, file = paste0("./results/parout_",Sys.Date(),".rds"))

#trace plot
library(ggplot2)
mcmc_trace(parout)+theme_bw()
ggplot2::ggsave(paste0("./results/traceplot_allo",Sys.Date(),".png"), device = "png", dpi = 500, height = 200, width = 300, units = 'mm')
#density plot
#bayesplot::mcmc_dens(parout)
#summary
as.data.frame(summary(parout))
parout$`mu_pred[1]`
#draws of par from posterior
parindout = as_draws_df(fit_vb$draws(c("par")))
t0pindout = as_draws_df(fit_vb$draws(c("t0p")))
#save par for individual plotting
saveRDS(parindout, file = paste0("./results/parindout_",Sys.Date(),".rds"))
saveRDS(t0pindout, file = paste0("./results/t0pindout_",Sys.Date(),".rds"))

# read in results ----
date = "2024-05-17"
# remember k is logit(k)
#all non individual based params
parout_in = readRDS(file = paste0('./results/parout_',date,'.rds'))
par_df<-as.data.frame(parout_in)
bayesplot::mcmc_dens(parout_in)
bayesplot::mcmc_trace(parout_in)+theme_bw()
summ_paroutin<-as.data.frame(summary(parout_in))

summ<-summ_paroutin%>%
  # mutate(median = case_when(
  #   variable == "mu_t" ~ exp(median),
  #   TRUE ~ median),
  #   q5 = case_when(
  #     variable == "mu_t" ~ exp(q5),
  #     TRUE ~ q5),
  #   q95 = case_when(
  #     variable == "mu_t" ~ exp(q95),
  #     TRUE ~ q95
  #   ))%>%
  mutate(`90%CI` = paste0(as.character(format(round(q5,3)),3),"–",as.character(format(round(q95,3)),3)),
         Median = round(median, 3))%>%
  dplyr::select(variable, Median,`90%CI`)

saveRDS(summ, file = "./results/summtable.rds")

max(summ_paroutin$rhat)
min(summ_paroutin$ess_bulk)
min(summ_paroutin$ess_tail)
max(summ_paroutin$ess_bulk)
max(summ_paroutin$ess_tail)

#max(summary(parout)$rhat)
#min(summary(parout)$ess_bulk)
#min(summary(parout)$ess_tail)
#individual based params
parindout_in = readRDS(file = paste0('./results/parindout_',date,'.rds'))
parindout_in_summ<-summary(parindout_in)

t0pindout_in = readRDS(file = paste0('./results/t0pindout_',date,'.rds'))
t0pindout_in_summ<-summary(t0pindout_in)
# report results ----

### plot a couple of individuals ----
induse = c(42,95,23)
#SN9 lineage
#induse = c(50,21,20,123,106,3,101)
#induse = ceiling(runif(4, 0, 143))
#induse = 4 #sunshine
ngrid = 101
agegrid = seq(from = -2, to = max(ij_b$age), length.out = ngrid)
nind = length(induse)
#pdf('indplots_bhdf.pdf', height = 8, width = 8)
#par(mfrow = c(2,2), mar = c(4, 4, 1, 1))
label = c("a","b","c")
#label = "b" #for sunshine
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
      scale_color_grey(start = 0.6, end = 0.1)+
      labs(title = k)
  
  ggplot2::ggsave(paste0("./Figures/indplot_",i,".png"), device = "png", dpi = 700, width = 100, height = 100, units = 'mm')
  
  }

# sunshine_plot<-sunshine+
#   geom_point(aes(x = c(1,1) , y = c(2.07,0.65)), color = "darkorange", size = 2.5, alpha = 0.8)
# 
# ggplot2::ggsave(paste0("./Figures/indplot_sunshine.png"), sunshine_plot, device = "png", dpi = 700, width = 100, height = 100, units = 'mm')

#dev.off()

parindout

#### need to adjust for new model outputs above
library(bayesplot)
#mcmc_intervals(parindout_in[1:143], outer_size = 0.5, inner_size = 1, point_size = 2)
#mcmc_intervals(parindout_in[144:286], outer_size = 0.5, inner_size = 1, point_size = 2)
#mcmc_intervals(parindout_in[287:429], outer_size = 0.5, inner_size = 1, point_size = 2)
#mcmc_intervals(parindout_in[430:ncol(parindout_in)], outer_size = 0.5, inner_size = 1, point_size = 2)

### L and k by birthyear----

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

agegrid = seq(from = -2, to = max(ij_b$age+20), length.out = ngrid)

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

ind_median%>%
  filter(year_zero <= 2013)

age_vb_y = matrix(NA,nrow(ind_median),ngrid)
age_vb_byr = matrix(NA,nrow(ind_median),ngrid)
age_vb_z = matrix(NA,nrow(ind_median),ngrid)

for (i in 1:nrow(ind_median)){
  
  for(j in 1:ngrid){
    age_vb_y[i,j] = ind_median$Ly[i]*(1-ind_median$ky[i]^(ind_median$t0p[i] + agegrid[j]))
    age_vb_byr[i,j] = agegrid[j] + ind_median$year_zero[i]
    age_vb_z[i,j] = ind_median$Lz[i]*(1-ind_median$kz[i]^(ind_median$t0p[i] + agegrid[j]))
  }  
}

by_df<-transform(expand.grid(i = seq(nrow(age_vb_byr)), j = seq(ncol(age_vb_byr))), year_by = c(age_vb_byr))
y_est<-transform(expand.grid(i = seq(nrow(age_vb_y)), j = seq(ncol(age_vb_y))), est = c(age_vb_y))
z_est<-transform(expand.grid(i = seq(nrow(age_vb_z)), j = seq(ncol(age_vb_z))), est = c(age_vb_z))
#y_est_BHDF<-transform(expand.grid(i = seq(nrow(age_vb_y_BHDF)), j = seq(ncol(age_vb_y_BHDF))), y_est = c(age_vb_y_BHDF))


growest_plot_TL<-by_df%>%
  left_join(y_est)%>%
  arrange(i,j)%>%
  mutate(agegrid = rep(agegrid,143))%>%
  left_join(ij_ID, by = c("i" = "ind"))%>%
  group_by(ID)%>%
  mutate(est_diff = lead(est)-est)%>%
  dplyr::rename("Pod" = "POD")%>%
  mutate(Length = "TL")

growest_plot_BHDF<-by_df%>%
  left_join(z_est)%>%
  arrange(i,j)%>%
  mutate(agegrid = rep(agegrid,143))%>%
  left_join(ij_ID, by = c("i" = "ind"))%>%
  group_by(ID)%>%
  mutate(est_diff = lead(est)-est)%>%
  dplyr::rename("Pod" = "POD")%>%
  mutate(Length = "BHDF")

growest_plot<-growest_plot_TL%>%
  bind_rows(growest_plot_BHDF)%>%
  dplyr::rename("Birth year" = age_value)

ggplot(growest_plot)+
  geom_histogram(mapping = aes(x = est_diff), binwidth = 0.001)

growest_plot%>%
  filter(est_diff <= 0.01)%>%
  group_by(ID, length)%>%
  mutate(rank = rank(j))%>%
  filter(rank == 1)%>%
  ungroup()%>%
  mutate(mean_limit = mean(agegrid))%>%
  distinct(mean_limit)

mu_pred_L1<-summ_paroutin%>%filter(variable == "mu_pred[1]")
mu_pred_L2<-summ_paroutin%>%filter(variable == "mu_pred[3]")

vbgc_zero<-ggplot(growest_plot)+
  annotate("rect", xmin = -5, xmax = 55, ymin = mu_pred_L1$q5, ymax =  mu_pred_L1$q95, fill = "red", alpha = 0.1)+
  annotate("rect", xmin = -5, xmax = 55, ymin = mu_pred_L2$q5, ymax =  mu_pred_L2$q95, fill = "red", alpha = 0.1)+
  geom_path(aes(x = agegrid, y = est, group = interaction(i, Length), color = Length), alpha = 0.6, linewidth = 0.2)+
  coord_cartesian(xlim=c(0, 41), ylim=c(0, 3.5))+
  geom_hline(yintercept = mu_pred_L1$median, color = "red")+
  geom_hline(yintercept = mu_pred_L2$median, color = "red")+
  geom_hline(yintercept = mu_pred_L1$q5, color = "red", linetype = "dashed")+
  geom_hline(yintercept = mu_pred_L2$q5, color = "red", linetype = "dashed")+
  geom_hline(yintercept = mu_pred_L1$q95, color = "red", linetype = "dashed")+
  geom_hline(yintercept = mu_pred_L2$q95, color = "red", linetype = "dashed")+
  #facet_wrap(~SEX)+
  theme_bw()+
  xlab("Age (years)")+
  ylab("Length (m)")+
  theme(legend.position = "none")+
  scale_color_grey(start = 0.6, end = 0.1)
  #geom_vline(xintercept = 7, color = "red")+
  #annotate("text", x=6.4, y=0, label="7", angle=0, color = "red")

ggplot2::ggsave(paste0("./Figures/vbgc_zero.png"), vbgc_zero, device = "png", dpi = 700, width = 250, height = 100, units = 'mm')

# ggplot(growest_plot)+
#   geom_line(aes(x = agegrid, y = y_est, group = as.factor(i)), alpha = 0.3)+
#   coord_cartesian(xlim=c(0, 45))+
#   coord_cartesian(ylim=c(0, 3.5))+
#   facet_wrap(~age_value)+
#   theme_bw()+
#   xlab("Age")+
#   ylab("Total length estimate (m)")

vbgc_by<-ggplot(growest_plot%>%filter(Length == "TL"))+
  geom_line(aes(x = year_by, y = est, group = interaction(i, Length), color = `Birth year`), alpha = 0.6, linewidth = 0.3)+
  xlim(c(1980,2050))+
  coord_cartesian(ylim=c(0, 3.5))+
  theme_bw()+
  xlab("Year")+
  ylab("")+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c('black','#ffb000'))
  #scale_color_viridis_d(begin = 0.1, end = 0.8)

ab<-ggpubr::ggarrange(vbgc_zero, vbgc_by, common.legend = F, legend = "bottom", widths = c(1,1.5), labels = "auto")

ggplot2::ggsave(paste0("./Figures/vbgcplot.png"), ab, device = "png", dpi = 700, width = 300, height = 150, units = 'mm')

ind_median%>%
  group_by(POD)%>%
  summarise(median = round(median(Ly),3), max = max(Ly), min = min(Ly))

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

uncertainty<-id_parsumm%>%ungroup()%>%arrange(param, year_zero)%>%group_by(year_zero)%>%
  mutate(rank = percent_rank(ind))%>%
  mutate(year_cat = as.numeric(year_zero+(rank)))%>%
  dplyr::rename(`Birth year` = "age_value",'Sex' = "SEX","Pod" = "POD")
  
uncertainty_Ly<-ggplot(uncertainty%>%filter(param == "Ly"))+
  geom_linerange(aes(x = year_cat, ymin = q95, ymax = q5, group = as.factor(i), color = Sex), linewidth = 0.5)+
  geom_point(aes(x = year_cat, y = median, group = as.factor(i), color = Sex, shape = Pod), size = 2, alpha = 0.6)+
  #facet_wrap(~param, ncol = 1, scales = "free")+
  theme_bw()+
  theme(legend.position = "bottom")+
  xlab("First year")+
  ylab("Length (m)")+
  scale_x_continuous(breaks = seq(1984, 2024, 1))+
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -90, vjust = -0.5))+
  #facet_wrap(~POD, ncol = 1)+
  scale_color_viridis_d(begin = 0, end = 0.9)
  #scale_color_manual(values = c('#648fff','#dc267f','#ffb000'))

ggplot2::ggsave(paste0("./Figures/Ly_CI.png"), uncertainty_Ly, device = "png", dpi = 700, height = 200, width = 300, units = 'mm')

curve_plot<-ggpubr::ggarrange(ab,uncertainty_Ly, ncol = 1, labels = c('','c'), heights = c(2,1.25))

ggplot2::ggsave(paste0("./Figures/curve_plots.png"), curve_plot, device = "png", dpi = 700, height = 225, width = 300, units = 'mm')

ggplot(uncertainty%>%filter(param == "logit_ky")%>%mutate(param = "ky"))+
  geom_point(aes(x = year_cat, y = median, group = as.factor(i), color = POD, shape = `Birth year`), size = 1.5, alpha = 0.6)+
  geom_linerange(aes(x = year_cat, ymin = median, ymax = median, group = as.factor(i), color = POD), linewidth = 0.2)+
  #facet_wrap(~param, ncol = 1)+
  theme_bw()+
  theme(legend.position = "bottom")+
  xlab("Year")+
  ylab("Growth rate")+
  scale_x_continuous(breaks = seq(1984, 2024, 1))+
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -90, vjust = -0.5))

ggplot2::ggsave(paste0("./Figures/ky_CI.png"), device = "png", dpi = 700, height = 100, width = 300, units = 'mm')

uncertainty%>%
  filter(year_zero < 2013 & param == "Ly")%>%
  mutate(year_group = case_when(
    year_zero <= 2007 ~ "Before",
    TRUE ~ "After"
  ))%>%
  group_by(year_group)%>%dplyr::summarise(mean = mean(median))

## box plots ----
#### sex ----

box<-ind_median%>%
  dplyr::select(-ky, -kz, -t0p)%>%
  tidyr::pivot_longer(cols = c(Ly, logit_ky, Lz, logit_kz), names_to = "param", values_to =  "median")%>%
  mutate(param2 = case_when(
    param == "Ly" ~ "L<sub>1i</sub>",
    param == "Lz" ~ "L<sub>2i</sub>",
    param == "logit_ky" ~ "logit(k<sub>1i</sub>)",
    param == "logit_kz" ~ "logit(k<sub>2i</sub>)"
  ))%>%
  mutate(group = case_when(
    param == "logit_ky" | param == "logit_kz" ~ "Growth rate parameters",
    TRUE ~ param2
  ))%>%
  dplyr::rename('Sex' = 'SEX', "Pod" = "POD")

sex_box<-ggplot(box)+
  geom_boxplot(aes(x = param2, y = median, fill = Sex), alpha = 0.5)+
  facet_wrap(.~as.factor(group), scales = "free", ncol = 4)+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(legend.position = "bottom")+
  scale_fill_viridis_d(begin = 0, end = 0.9)+
  theme(strip.text = ggtext::element_markdown(),
        axis.text.x = ggtext::element_markdown())
  

#### pod ----
pod_box<-ggplot(box)+
  geom_boxplot(aes(x = param2, y = median, fill = Pod), alpha = 0.6)+
  facet_wrap(~group, scales = "free")+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(legend.position = "bottom")+
  theme(strip.text = ggtext::element_markdown(),
        axis.text.x = ggtext::element_markdown())

ggpubr::ggarrange(pod_box, sex_box, ncol = 1, labels = "auto")

ggplot2::ggsave(paste0("./Figures/boxplots.png"), device = "png", dpi = 700, height = 175, width = 150, units = 'mm')

### Ly given Lz ----

muLy<-summ_paroutin%>%
  filter(variable == "mu[1]")
muLy
muky<-summ_paroutin%>%
  filter(variable == "mu[2]")
muky
muLz<-summ_paroutin%>%
  filter(variable == "mu[3]")
muLz
mukz<-summ_paroutin%>%
  filter(variable == "mu[4]")

vc11<-summ_paroutin%>%
  filter(variable == "varcov_par[1,1]")
vc13<-summ_paroutin%>%
  filter(variable == "varcov_par[1,3]")
vc33<-summ_paroutin%>%
  filter(variable == "varcov_par[3,3]")

vc24<-summ_paroutin%>%
  filter(variable == "varcov_par[2,4]")
vc22<-summ_paroutin%>%
  filter(variable == "varcov_par[2,2]")
vc44<-summ_paroutin%>%
  filter(variable == "varcov_par[4,4]")

vc12<-summ_paroutin%>%
  filter(variable == "varcov_par[1,2]")
vc22<-summ_paroutin%>%
  filter(variable == "varcov_par[2,2]")
vc11<-summ_paroutin%>%
  filter(variable == "varcov_par[1,1]")



#regression coefficient / slope
slope_L<-vc13$median/vc33$median
y_intercept_L<-muLy$median - slope_L*muLz$median
#conditional variance
var_L<-vc11$median-slope*vc13$median
sd_L<-sqrt(var_L)

ij_est<-ij_all%>%
  mutate(Length_est = slope_L*BHDF + y_intercept_L)%>%
  mutate(`Age bin` = case_when(
    age < 10 ~ "<10",
    age >= 10 & age < 20 ~ "10-20",
    age >= 20 & age < 30 ~ "20-30",
    age >= 30 & age < 40 ~ "30-40",
    age >= 40 ~ "40+"
  ))

#look at Ly and Lz to see if this equation works
Lyest<-id_parsumm%>%
  filter(param == "Lz" | param == "Ly")%>%
  distinct(param, ind, ID, mean, median, sd, mad)%>%
  mutate(Length_est = case_when(
    param == "Lz" ~ slope_L*median + y_intercept_L,
    TRUE ~ NA))

x = Lyest%>%filter(param == "Lz")%>%dplyr::select(ind, median)
y = Lyest%>%filter(param == "Ly")%>%dplyr::select(ind, median)

join<-x%>%
  left_join(y, by = "ind")

linear_plot<-ggplot()+
  geom_point(ij_est, mapping = aes(x = BHDF, y = length, color = `Age bin`), alpha = 0.6, size = 2)+
  geom_point(join, mapping = aes(x = median.x, y = median.y), alpha = 0.6, size = 2)+
  geom_abline(mapping = aes(slope = slope_L, intercept = y_intercept_L), color = "red")+
  geom_abline(slope = slope_L, intercept = y_intercept_L-sd_L, color = "red", linetype = "dashed") +
  geom_abline(slope = slope_L, intercept = y_intercept_L+sd_L, color = "red", linetype = "dashed") +
  scale_color_viridis_d(option="plasma", end = 0, begin = 0.88)+
  theme_bw()+
  theme(legend.position = c(0.85,0.25))+
  ylab("TL (m)")+
  xlab("BHDF (m)")+
  coord_fixed(ratio = 0.4)

ggplot2::ggsave(paste0("./Figures/cond_BHDF_TL.png"), linear_plot, device = "png", dpi = 700, height = 150, width = 150, units = 'mm')

## three plots ----
library(ggExtra)

slope_L
y_intercept_L
sd_L
eq_L<-expression(paste(hat(y) == 2.13*x + 0.93,", ",sigma == 0.07))#,round(slope_L,3),"* L['2i'] "))#,round(y_intercept_L,3),", \U03c3 = ",round(var_L,3))

LzLy<-ggplot()+
  geom_point(ind_median, mapping = aes(x = Lz, y = Ly), color = "black", alpha = 0.9)+
  geom_abline(slope = slope_L, intercept = y_intercept_L-sd_L, color = "red", linetype = "dashed") +
  geom_abline(slope = slope_L, intercept = y_intercept_L+sd_L, color = "red", linetype = "dashed") +
  geom_abline(mapping = aes(slope = slope_L, intercept = y_intercept_L), color = "red", alpha = 0.9)+
  theme_bw()+
  xlab(bquote(L['2i']))+
  ylab(bquote(L['1i']))+
  annotate(geom = "text", x=0.89, y=3.3, label=eq_L, parse = F, size = 3)+
  coord_fixed(ratio = 0.4)

###
#regression coefficient / slope
slope_k<-vc24$median/vc22$median
slope_k
y_intercept_k<-mukz$median - slope_k*muky$median
y_intercept_k
#conditional variance
var_k<-vc44$median-slope_k*vc24$median
sd_k<-sqrt(var_k)
sd_k

eq_k<-expression(paste(hat(y) == 1.2,'0',x - 0.17,", ",sigma == 0.14))#,round(slope_L,3),"* L['2i'] "))#,round(y_intercept_L,3),", \U03c3 = ",round(var_L,3))

kykz<-ggplot()+
  geom_point(ind_median, mapping = aes(x = logit_ky, y = logit_kz), color = "black", alpha = 0.9)+
  geom_abline(slope = slope_k, intercept = y_intercept_k-sd_k, color = "red", linetype = "dashed") +
  geom_abline(slope = slope_k, intercept = y_intercept_k+sd_k, color = "red", linetype = "dashed") +
  geom_abline(mapping = aes(slope = slope_k, intercept = y_intercept_k), color = "red", alpha = 0.9)+
  theme_bw()+
  xlab(bquote(logit(k['1i'])))+
  ylab(bquote(logit(k['2i'])))+
  annotate(geom = "text", x=0.32, y=1.6, label=eq_k, parse = F, size = 3)+
  xlim(c(-0.1, 1.6))+
  ylim(c(-0.1, 1.6))+
  coord_fixed(ratio = 1)

###
#regression coefficient / slope
slope_Lk<-vc12$median/vc22$median
slope_Lk
y_intercept_Lk<-muLy$median - slope_Lk*muky$median
y_intercept_Lk
#conditional variance
var_Lk<-vc11$median-slope_Lk*vc12$median
sd_Lk<-sqrt(var_Lk)
sd_Lk

eq_Lk<-expression(paste(hat(y) == 0.26*x + 2.74,", ",sigma == 0.1,'0'))#,round(slope_L,3),"* L['2i'] "))#,round(y_intercept_L,3),", \U03c3 = ",round(var_L,3))

kyLy<-ggplot()+
  geom_point(ind_median, mapping = aes(x = logit_ky, y = Ly), color = "black", alpha = 0.9)+
  geom_abline(slope = slope_Lk, intercept = y_intercept_Lk-sd_Lk, color = "red", linetype = "dashed") +
  geom_abline(slope = slope_Lk, intercept = y_intercept_Lk+sd_Lk, color = "red", linetype = "dashed") +
  geom_abline(mapping = aes(slope = slope_Lk, intercept = y_intercept_Lk), color = "red", alpha = 0.9)+
  theme_bw()+
  xlab(bquote(logit(k['1i'])))+
  ylab(bquote(L['1i']))+
  annotate(geom = "text", x=0.50, y=3.3, label=eq_Lk, parse = F, size = 3)

params_plot<-ggpubr::ggarrange(kyLy, kykz, LzLy, labels = "auto")
ggplot2::ggsave(paste0("./Figures/params.png"), params_plot, device = "png", dpi = 700, height = 200, width = 200, units = 'mm', bg="white")

summ_paroutin%>%
  filter(grepl("mu",variable))

summ_paroutin%>%
  filter(variable == "mu_pred[2]")%>%
  mutate(across(2:q95, ~(1/(1+exp(-.x)))))

summ_paroutin%>%
  filter(variable == "mu_pred[4]")%>%
  mutate(across(2:q95, ~(1/(1+exp(-.x)))))

         