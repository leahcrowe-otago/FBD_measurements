library(condMVNorm)

# works for p(Ly | Lz) -----

vc11<-summ_paroutin%>%
  filter(variable == "varcov_par[1,1]")

vc12<-summ_paroutin%>%
  filter(variable == "varcov_par[1,3]")

vc22<-summ_paroutin%>%
  filter(variable == "varcov_par[3,3]")

varcov = matrix(NA, 2, 2)

varcov[1,1] = vc11$median
varcov[2,2] = vc22$median
varcov[1,2] = vc12$median
varcov[2,1] = varcov[1,2]

varcov

mu1 = muLy$median
mu2 = muLz$median

# doesn't work, trying to do p(y | z) -----
vc11<-summ_paroutin%>%
  filter(variable == "Lobs[1,1]")

vc12<-summ_paroutin%>%
  filter(variable == "Lobs[2,1]")

vc22<-summ_paroutin%>%
  filter(variable == "Lobs[2,2]")

varcov[1,1] = vc11$median^2
varcov[2,2] = vc12$median^2+vc22$median^2
varcov[1,2] = vc11$median*vc12$median
varcov[2,1] = varcov[1,2]

s1<-vc11$median
s2<-sqrt(varcov[2,2])

vc11$median*vc12$median/(s1*s2)

mu1 = ij_all$length[1]
mu2 = ij_all$BHDF[1]

mu1 + rho_obs$median*s1/s2*(ij_all$BHDF[8] - mu2)

# data and formula ----

has_z<-ij_all%>%filter(!is.na(BHDF))

Lyest = NULL
convar = NULL

for (i in 1:nrow(has_z)){
  cond_est=condMVN(mean = c(mu1, mu2), sigma = varcov, dependent.ind = 1, given.ind = 2, X.given = has_z$BHDF[i])
  #print(cond_est)
  
  Lyest[i] = cond_est$condMean
  convar[i] = cond_est$condVar
  }


has_z<-has_z%>%
  mutate(y_est = Lyest,
         var_est = convar)

has_z

ggplot()+
  geom_point(has_z, mapping = aes(x = BHDF, y = y_est))+
  geom_point(has_z, mapping = aes(x = BHDF, y = length), color = "red")+
  geom_point(mapping = aes(x = x$median, y = y$median))

lm<-lm(has_z$y_est ~ has_z$BHDF)
lm_intercept = lm$coefficients[[1]]
lm_slope = lm$coefficients[[2]]

exp<-lm(log(has_z$y_est) ~ has_z$BHDF)
exp_intercept = exp$coefficients[[1]]
exp_slope = exp$coefficients[[2]]


model_est<-has_z%>%
  mutate(lm_est = BHDF*slope + intercept,
         exp_est = exp(exp_intercept)*exp(exp_slope)^BHDF
         )%>%
  arrange(age)

library("ggExtra")
lm_plot<-ggplot()+
  geom_point(model_est, mapping = aes(x = BHDF, y = length))+
  geom_point(model_est, mapping = aes(x = BHDF, y = lm_est))+
  geom_abline(mapping = aes(slope = lm_slope,intercept = lm_intercept, color =  "Linf"), size = 1.5)+
  theme(legend.position = "bottom")

ggMarginal(lm_plot)+
  geom_point(mapping = aes(x = x$median, y = y$median, color = "Linf"))
  
  

