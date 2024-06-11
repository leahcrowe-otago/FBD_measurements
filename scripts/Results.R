library(bayesplot)
library(latex2exp)
library(dplyr)
library(ggplot2)
library(ggExtra)

# read in results ----
date = "2024-05-17"

# k is logit(k)
# all population-level based params ----
parout_in = readRDS(file = paste0('./results/parout_',date,'.rds'))
par_df<-as.data.frame(parout_in)
#bayesplot::mcmc_dens(parout_in)
#bayesplot::mcmc_trace(parout_in)+theme_bw()
summ_paroutin<-as.data.frame(summary(parout_in))

summ_paroutin%>%
  filter(grepl("mu",variable))

#### Table S2 ----
summ<-summ_paroutin%>%
  mutate(`90%CI` = paste0(as.character(format(round(q5,3)),3),"–",as.character(format(round(q95,3)),3)),
         Median = round(median, 3))%>%
  dplyr::select(variable, Median,`90%CI`)

saveRDS(summ, file = "./results/summtable.rds")

#### Rhat, ESS ----
max(summ_paroutin$rhat)
min(summ_paroutin$ess_bulk)
min(summ_paroutin$ess_tail)
max(summ_paroutin$ess_bulk)
max(summ_paroutin$ess_tail)

# individual-level estimates ----
parindout_in = readRDS(file = paste0('./results/parindout_',date,'.rds'))
parindout_in_summ<-summary(parindout_in)

t0pindout_in = readRDS(file = paste0('./results/t0pindout_',date,'.rds'))
t0pindout_in_summ<-summary(t0pindout_in)

# report results ----

### Fig. S7, plot a couple of individuals ----
induse = c(42,95,23)
#induse = 4 #sunshine

ngrid = 101
agegrid = seq(from = -2, to = max(ij_b$age), length.out = ngrid)
nind = length(induse)

label = c("a","b","c")
#label = "b" #for sunshine
for(i in 1:nind){
  k = label[i]
  i = induse[i]
  print(i)
  x = ij_all%>%filter(ind == i)%>%select(age)
  y = ij_all%>%filter(ind == i)%>%select(length)
  z = ij_all%>%filter(ind == i)%>%select(BHDF)
  
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

### Fig. 2, VBGCs, L and k by birthyear----

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

#### Fig. 2a ----
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
  theme_bw()+
  xlab("Age (years)")+
  ylab("Length (m)")+
  theme(legend.position = "none")+
  scale_color_grey(start = 0.6, end = 0.1)

ggplot2::ggsave(paste0("./Figures/vbgc_zero.png"), vbgc_zero, device = "png", dpi = 700, width = 250, height = 100, units = 'mm')

#### Fig. 2b ----
vbgc_by<-ggplot(growest_plot%>%filter(Length == "TL"))+
  geom_line(aes(x = year_by, y = est, group = interaction(i, Length), color = `Birth year`), alpha = 0.6, linewidth = 0.3)+
  xlim(c(1980,2050))+
  coord_cartesian(ylim=c(0, 3.5))+
  theme_bw()+
  xlab("Year")+
  ylab("")+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c('black','#ffb000'))

ab<-ggpubr::ggarrange(vbgc_zero, vbgc_by, common.legend = F, legend = "bottom", widths = c(1,1.5), labels = "auto")

ggplot2::ggsave(paste0("./Figures/vbgcplot.png"), ab, device = "png", dpi = 700, width = 300, height = 150, units = 'mm')

## individual median summaries
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

#### Fig. 2c ----
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
  scale_color_viridis_d(begin = 0, end = 0.9)

ggplot2::ggsave(paste0("./Figures/Ly_CI.png"), uncertainty_Ly, device = "png", dpi = 700, height = 200, width = 300, units = 'mm')

#### Fig 2, together ----
curve_plot<-ggpubr::ggarrange(ab,uncertainty_Ly, ncol = 1, labels = c('','c'), heights = c(2,1.25))

ggplot2::ggsave(paste0("./Figures/curve_plots.png"), curve_plot, device = "png", dpi = 700, height = 225, width = 300, units = 'mm')

# ggplot(uncertainty%>%filter(param == "logit_ky")%>%mutate(param = "ky"))+
#   geom_point(aes(x = year_cat, y = median, group = as.factor(i), color = Sex, shape = Pod), size = 1.5, alpha = 0.6)+
#   geom_linerange(aes(x = year_cat, ymin = median, ymax = median, group = as.factor(i), color = Sex), linewidth = 0.2)+
#   #facet_wrap(~param, ncol = 1)+
#   theme_bw()+
#   theme(legend.position = "bottom")+
#   xlab("Year")+
#   ylab("Growth rate")+
#   scale_x_continuous(breaks = seq(1984, 2024, 1))+
#   theme(panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -90, vjust = -0.5))
# 
# ggplot2::ggsave(paste0("./Figures/ky_CI.png"), device = "png", dpi = 700, height = 100, width = 300, units = 'mm')

#summary on birth year cutoff
uncertainty%>%
  filter(year_zero < 2013 & param == "Ly")%>%
  mutate(year_group = case_when(
    year_zero <= 2007 ~ "Before",
    TRUE ~ "After"
  ))%>%
  group_by(year_group)%>%dplyr::summarise(mean = mean(median))

## Fig. 4, box plots ----

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

#### Fig. 4a, pod ----
pod_box<-ggplot(box)+
  geom_boxplot(aes(x = param2, y = median, fill = Pod), alpha = 0.6)+
  facet_wrap(~group, scales = "free")+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(legend.position = "bottom")+
  theme(strip.text = ggtext::element_markdown(),
        axis.text.x = ggtext::element_markdown())

#### Fig. 4b, sex ----
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

#### Fig. 4, together ----
ggpubr::ggarrange(pod_box, sex_box, ncol = 1, labels = "auto")

ggplot2::ggsave(paste0("./Figures/boxplots.png"), device = "png", dpi = 700, height = 175, width = 150, units = 'mm')

### Fig. 3, y given x

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

## Fig. 3, three conditional distribution plots ----

#### Fig. 3a, L_1 | logit(k_1) ----
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

#### Fig. 3b, logit(k_2) | logit(k_1) ----
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

#### Fig. 3c, Ly given Lz ----
#regression coefficient / slope
slope_L<-vc13$median/vc33$median
y_intercept_L<-muLy$median - slope_L*muLz$median
#conditional variance
var_L<-vc11$median-slope*vc13$median
sd_L<-sqrt(var_L)
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

#### Fig. 3, together ----
params_plot<-ggpubr::ggarrange(kyLy, kykz, LzLy, labels = "auto")
ggplot2::ggsave(paste0("./Figures/params.png"), params_plot, device = "png", dpi = 700, height = 200, width = 200, units = 'mm', bg="white")
