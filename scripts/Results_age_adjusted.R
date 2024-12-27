library(bayesplot)
library(latex2exp)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(cmdstanr)

# data ----
#IDs
ij_ID = readRDS(file = './data/Measurements/ij_ID.rds')
n_ind = nrow(ij_ID)

#both bhdf & tl
ij_b = readRDS(file = './data/Measurements/ij_1.rds')
#bhdf only
ij_z = readRDS(file = './data/Measurements/ij_2.rds')
#total length only
ij_y = readRDS(file = './data/Measurements/ij_3.rds')

ij_all<-ij_b%>%
  bind_rows(ij_z)%>%
  bind_rows(ij_y)%>%
  mutate(age_add = 0)

#run below to adjust minimum age
age_est_add<-5 # add time to age (years), use same value from the model run
source('./scripts/stan/adjust_age_minyr.R', local = TRUE, verbose = F)$value

#these ij_* have been changed in the above script
ij_all_5<-ij_b%>%
  bind_rows(ij_z)%>%
  bind_rows(ij_y)%>%
  mutate(age_add = 5)

ij_all<-ij_all%>%
  bind_rows(ij_all_5)

ij_ID%>%
  group_by(POD,SEX)%>%
  tally()

# read in results ----
# from runstan_allo_mv_t0
date = "2024-06-21"
date_5 = "2024-12-13_5" # add 5 years to unknown birth year animals
date_10 = "2024-12-13_10" # add 5 years to unknown birth year animals

# k is logit(k)
# all population-level based params ----
parout_in = readRDS(file = paste0('./results/parout_',date,'.rds'))
parout_in_5 = readRDS(file = paste0('./results/parout_',date_5,'.rds'))
parout_in_10 = readRDS(file = paste0('./results/parout_',date_10,'.rds'))

par_df<-as.data.frame(parout_in)
par_df_5<-as.data.frame(parout_in_5)
par_df_10<-as.data.frame(parout_in_10)
#bayesplot::mcmc_dens(parout_in)
#bayesplot::mcmc_trace(parout_in)+theme_bw()
summ_paroutin<-as.data.frame(summary(parout_in))
summ_paroutin_5<-as.data.frame(summary(parout_in_5))

#### Table S2 ----
summ_5<-summ_paroutin_5%>%
  mutate(Median = round(median, 3),
         '5th percentile' = format(round(q5,3)),
         '95th percentile' = format(round(q95,3)))%>%
  dplyr::select(variable, Median,'5th percentile','95th percentile')

saveRDS(summ_5, file = "./results/summtable_p5.rds")

#### Rhat, ESS ----
max(summ_paroutin_5$rhat)
min(summ_paroutin_5$ess_bulk)
min(summ_paroutin_5$ess_tail)
max(summ_paroutin_5$ess_bulk)
max(summ_paroutin_5$ess_tail)

t0p_est<-summ_paroutin%>%
  filter(variable == "t0p")
t0p = t0p_est$median

t0p_est_5<-summ_paroutin_5%>%
  filter(variable == "t0p")
t0p_5 = t0p_est_5$median

# individual-level estimates ----
parindout_in = readRDS(file = paste0('./results/parindout_',date,'.rds'))%>%mutate(age_add = 0)
parindout_in_summ<-summary(parindout_in)

parindout_in_5 = readRDS(file = paste0('./results/parindout_',date_5,'.rds'))%>%mutate(age_add = 5)
parindout_in_summ_5<-summary(parindout_in_5)

parindout_in<-parindout_in%>%
  bind_rows(parindout_in_5)

# report results ----
 

### Fig. 3, VBGCs, L and k by birthyear----

#rename logit_ks
id_parsumm<-parindout_in_summ_5%>%
  mutate(param = as.factor(
    case_when(
      grepl(",1]", variable) ~ 'Ly',
      grepl(",2]", variable) ~ 'logit_ky',
      grepl(",3]", variable) ~ 'Lz',
      grepl(",4]", variable) ~ 'logit_kz'
    )))%>%
  mutate(ind = as.numeric(stringr::str_extract(substr(variable, 5, nchar(variable)), '[^,]+')))%>%
  left_join(ij_ID, by = 'ind')

ngrid = 101

agegrid = seq(from = -2, to = max(ij_b$age+20), length.out = ngrid)

#inverse logit ks to just get k on 0-1 scale, k = e^-K
ind_median<-id_parsumm%>%
  group_by(ind)%>%
  tidyr::pivot_wider(names_from = "param", values_from = "median")%>%
  group_by(ind)%>%
  tidyr::fill(Ly,logit_ky,Lz,logit_kz, .direction = "downup")%>% #t0p
  #logit(k) = log(k/(1-k)), k = exp(-K)
  mutate(ky = 1/(1+exp(-logit_ky)),
         kz = 1/(1+exp(-logit_kz)))%>%
  mutate(Ky = -(log(ky)),
         Kz = -(log(kz)),
         exp_K = exp(-(-(log(kz)))))%>%
  distinct(ind, ID, year_zero, age_value, SEX, POD, Ly, Lz, logit_ky, ky,  Ky, logit_kz, kz, Kz, exp_K)

## Fig. 5, box plots ----

box<-ind_median%>%
  dplyr::select(-ky, -kz)%>% #t0p
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

box_Ly_F<-box%>%filter(Sex == "F" & param == "Ly")%>%ungroup()
plot_Ly_F<-quantile(box_Ly_F$median, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
box_Ly_M<-box%>%filter(Sex == "M" & param == "Ly")%>%ungroup()
plot_Ly_M<-quantile(box_Ly_M$median, probs = c(0.05, 0.30, 0.5, 0.70, 0.95))

#IQR1.5
IQR1.5_F<-(plot_Ly_F[4]-plot_Ly_F[2])*1.5
plot_Ly_F[2]-IQR1.5_F
plot_Ly_F[4]+IQR1.5_F

IQR1.5_M<-(plot_Ly_M[4]-plot_Ly_M[2])*1.5
plot_Ly_M[2]-IQR1.5_M
plot_Ly_M[4]+IQR1.5_M

#### Fig. 5a, pod ----
pod_box<-ggplot(box)+
  geom_boxplot(aes(x = param2, y = median, fill = Pod), alpha = 0.6)+
  facet_wrap(~group, scales = "free")+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(legend.position = "bottom")+
  theme(strip.text = ggtext::element_markdown(),
        axis.text.x = ggtext::element_markdown())

ggplot_build(pod_box)$data

#### Fig. 5b, sex ----
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

ggplot_build(sex_box)$data
#### Fig. 5, together ----
ggpubr::ggarrange(pod_box, sex_box, ncol = 1, labels = "auto")

ggplot2::ggsave(paste0("./Figures/plus_age/boxplots_p5.png"), device = "png", dpi = 300, height = 175, width = 150, units = 'mm')

## Caterpillar plots

library(bayesplot)
library(cmdstanr)
library(dplyr)

options(width = 200)

# all population-level based params ----
parout_in_0 = readRDS(file = paste0('./results/parout_2024-06-21.rds'))#%>%mutate(age_add = 0)
parout_in_5 = readRDS(file = paste0('./results/parout_2024-12-13_5.rds'))#%>%mutate(age_add = 5)
parout_in_10 = readRDS(file = paste0('./results/parout_2024-12-13_10.rds'))#%>%mutate(age_add = 10)
## caterpillar plots ----

p0<-apply(parout_in_0, 2, quantile, c(0.05,0.25, 0.5, 0.75, 0.95), na.rm = T)[,1:18]
p5<-apply(parout_in_5, 2, quantile, c(0.05,0.25, 0.5, 0.75, 0.95), na.rm = T)[,1:18]
p10<-apply(parout_in_10, 2, quantile, c(0.05,0.25, 0.5, 0.75, 0.95), na.rm = T)[,1:18]

ggplot()+
  geom_linerange(mapping = aes(xmin = matrix(p0[1,]), xmax = matrix(p0[5,]), y = c(seq(18,1,-1))), size = 1)+
  geom_linerange(mapping = aes(xmin = matrix(p5[1,]), xmax = matrix(p5[5,]), y = c(seq(18.2,1.2,-1))), size = 1, color = "red")+
  geom_linerange(mapping = aes(xmin = matrix(p10[1,]), xmax = matrix(p10[5,]), y = c(seq(18.4,1.4,-1))), size = 1, color = "blue")+
  geom_linerange(mapping = aes(xmin = matrix(p0[2,]), xmax = matrix(p0[4,]), y = c(seq(18,1,-1))), size = 3, alpha = 0.8)+
  geom_linerange(mapping = aes(xmin = matrix(p5[2,]), xmax = matrix(p5[4,]), y = c(seq(18.2,1.2,-1))), size = 3, color = "red", alpha = 0.8)+
  geom_linerange(mapping = aes(xmin = matrix(p10[2,]), xmax = matrix(p10[4,]), y = c(seq(18.4,1.4,-1))), size = 3, color = "blue", alpha = 0.8)+
  geom_point(mapping = aes(x = matrix(p0[3,]), y = c(seq(18,1,-1))), size = 2, fill = "white", color = "black", shape = 21, alpha = 0.7)+
  geom_point(mapping = aes(x = matrix(p5[3,]), y = c(seq(18.2,1.2,-1))), size = 2, fill = "white", color = "red", shape = 21, alpha = 0.7)+
  geom_point(mapping = aes(x = matrix(p10[3,]), y = c(seq(18.4,1.4,-1))), size = 2, fill = "white", color = "blue", shape = 21, alpha = 0.7)+
  theme_bw()+
  xlab('')+
  ylab('')+
  scale_y_continuous(breaks=c(seq(18.2,1.2,-1)),
                     labels=c(expression(mu[L[1]]),expression(mu[k[1]]),expression(mu[L[2]]),expression(mu[k[2]]),
                              expression(sigma[L[1]]),expression(sigma[k[1]]),expression(sigma[L[2]]),expression(sigma[k[2]]),
                              expression(sigma[y[1]]),expression(sigma[y[2]]),expression(rho[y1y2]), expression(rho[L1k1]),
                              expression(rho[L1L2]),expression(rho[L1k2]),expression(rho[k1L2]),expression(rho[k1k2]),
                              expression(rho[L2k2]),expression(t[0])))+
  xlim(c(0,3.5))+
  theme(panel.grid.minor.y = element_blank() )

ggplot2::ggsave(paste0("./Figures/plus_age/catplots510.png"), device = "png", dpi = 400, height = 225, width = 200, units = 'mm')

names(parout_in_0)[1:37]
seq(12.1, 1.1, by = -1)
