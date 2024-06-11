

data {
  int<lower=0> n_ind; // number individuals
  int<lower=0> n_k; // number of individuals of known sex
  int<lower=0> n_obs_k; // number of observations of known sex
  int<lower=0> n_obs; // number of observations
  
  int<lower=0> J; // number of params (Ly, ky, Lz, kz), dim of multivariate matrix
  
  array[n_obs] int id; // IDs k (known sex), u (unknown sex)
  
  vector[n_obs] age; // age k (known sex), u (unknown sex)
  
  array[n_obs] int type; // type k (known sex), u (unknown sex)
  
  matrix[n_obs,2] val; // length values k (known sex), u (unknown sex)
  
  array[n_k] int sex; // sex
  
  array[n_ind] int pod; // pod
 
 }

parameters {
  vector[n_ind] z_t;
  real<lower=0> mu_t;
  real<lower=0> sigma_t;
  
  matrix[n_ind, 4] z;
  vector<lower=0>[J] sigma;
  cholesky_factor_corr[J] L;
  vector[J] mu; 
  vector[J] beta; // effect of sex = 1
  vector[J] gamma; // effect of pod = 1
  vector<lower=0>[2] sigma_obs;
  real<lower=0,upper=1> rho_obs;
  real<lower=0,upper=1> pi;
  
 }

transformed parameters{
  matrix[n_k, 4] par_k;
  matrix[n_ind, 4] par_u_0;
  matrix[n_ind, 4] par_u_1;
  vector[n_ind] t0p;
  matrix[2,2] varcov;
  cholesky_factor_cov[2] Lobs;
  
  for(i in 1:n_k){
    par_k[i] = (mu + beta*sex[i] + gamma*pod[i] + diag_pre_multiply(sigma,L)*z[i]')';
    par_u_0[i] = mu';
    par_u_1[i] = mu';
  }
  
  for(i in (n_k+1):n_ind){
    par_u_0[i] = (mu + gamma*pod[i] + diag_pre_multiply(sigma,L)*z[i]')';
    par_u_1[i] = (mu + beta + gamma*pod[i] + diag_pre_multiply(sigma,L)*z[i]')';
  }
  
  for(i in 1:n_ind){
    t0p[i] = exp(mu_t + z_t[i]*sigma_t);
  }
  
  varcov[1,1] = sigma_obs[1]^2;
  varcov[2,2] = sigma_obs[2]^2;
  varcov[1,2] = sigma_obs[1]*sigma_obs[2]*rho_obs;
  varcov[2,1] = varcov[1,2];
  Lobs = cholesky_decompose(varcov);

 }

model {

 matrix[n_obs,2] obs_mean;
 matrix[n_obs,2] obs_mean_0;
 matrix[n_obs,2] obs_mean_1;
 matrix[n_ind,2] tmp = rep_matrix(0,n_ind,2); 
 
 for(i in 1:n_ind){
    z[i] ~ std_normal();
    z_t[i] ~ std_normal();
 }
 
 for(i in 1:n_k){
    sex[i] ~ bernoulli(pi);
 }
 
 // observations with known sex

 for(i in 1:n_obs_k){
    if(type[i] == 1){ // both
      obs_mean[i,1] = par_k[id[i],1]*(1-inv_logit(par_k[id[i],2])^(age[i] + t0p[id[i]]));
      obs_mean[i,2] = par_k[id[i],3]*(1-inv_logit(par_k[id[i],4])^(age[i] + t0p[id[i]]));
      target += multi_normal_cholesky_lpdf(val[i]| obs_mean[i], Lobs);
    } else if(type[i] == 2) { // z
      target += normal_lpdf(val[i,2]| par_k[id[i],3]*(1-inv_logit(par_k[id[i],4])^(age[i] + t0p[id[i]])), sigma_obs[2]);
    } else if(type[i] == 3) { // y
      target += normal_lpdf(val[i,1]| par_k[id[i],1]*(1-inv_logit(par_k[id[i],2])^(age[i] + t0p[id[i]])), sigma_obs[1]);
    }
    
 }
 
 // observations with unknown sex
 
 for(i in (n_obs_k+1):n_obs){
    if(type[i] == 1){ // both
      obs_mean_0[i,1] = par_u_0[id[i],1]*(1-inv_logit(par_u_0[id[i],2])^(age[i] + t0p[id[i]]));
      obs_mean_0[i,2] = par_u_0[id[i],3]*(1-inv_logit(par_u_0[id[i],4])^(age[i] + t0p[id[i]]));
      obs_mean_1[i,1] = par_u_1[id[i],1]*(1-inv_logit(par_u_1[id[i],2])^(age[i] + t0p[id[i]]));
      obs_mean_1[i,2] = par_u_1[id[i],3]*(1-inv_logit(par_u_1[id[i],4])^(age[i] + t0p[id[i]]));
      tmp[id[i],1] += multi_normal_cholesky_lpdf(val[i]| obs_mean_0[i], Lobs);
      tmp[id[i],2] += multi_normal_cholesky_lpdf(val[i]| obs_mean_1[i], Lobs);
    } else if(type[i] == 2) { // y
      tmp[id[i],1] += normal_lpdf(val[i,1]| par_u_0[id[i],1]*(1-inv_logit(par_u_0[id[i],2])^(age[i] + t0p[id[i]])), sigma_obs[1]);
      tmp[id[i],2] += normal_lpdf(val[i,1]| par_u_1[id[i],1]*(1-inv_logit(par_u_1[id[i],2])^(age[i] + t0p[id[i]])), sigma_obs[1]);
    } else if(type[i] == 3) { // z
      tmp[id[i],1] += normal_lpdf(val[i,2]| par_u_0[id[i],3]*(1-inv_logit(par_u_0[id[i],4])^(age[i] + t0p[id[i]])), sigma_obs[2]);
      tmp[id[i],2] += normal_lpdf(val[i,2]| par_u_1[id[i],3]*(1-inv_logit(par_u_1[id[i],4])^(age[i] + t0p[id[i]])), sigma_obs[2]);
    }
 }
 
 for(i in (n_k+1):n_ind){
    tmp[i,1] += log(1-pi);
    tmp[i,2] += log(pi);
    
    target += log_sum_exp(tmp[i]);
 }
 

//priors
  //t0p ~ lognormal(0, 10);
  mu_t ~ normal(0, 1000);
  sigma_t ~ student_t(3, 0, 50);
  
  mu ~ normal(0, 1000);
  beta ~ normal(0, 1000);
  gamma ~ normal(0, 1000);
  sigma ~ student_t(3, 0, 50);
  L ~ lkj_corr_cholesky(1);
  
  for (i in 1:2){
    sigma_obs[i] ~ student_t(3, 0, 50);
  }

 }


generated quantities {

  corr_matrix[J] corr;
  cov_matrix[J] varcov_par;
  //below is L*L'
  corr = multiply_lower_tri_self_transpose(L);
  //below is diag(sigma)*LL'*diag(sigma)
  varcov_par = quad_form_diag(corr, sigma);
  cholesky_factor_cov[4] L_covar = cholesky_decompose(varcov_par);
  //ppd
  vector[J] mu_pred = multi_normal_cholesky_rng(mu, L_covar);

 }

