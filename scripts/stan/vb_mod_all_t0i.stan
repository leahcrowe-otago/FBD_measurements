
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

  // real<lower=0> t0p;
  // vector<lower=0>[n_ind] t0p;
  vector[n_ind] z_t;
  real mu_t;
  real<lower=0> sigma_t;
  matrix[n_ind, J] z;
  //matrix[n_ind, J] par;
  vector<lower=0>[J] sigma;
  cholesky_factor_corr[J] L;
  vector[J] mu;
  vector<lower=0>[2] sigma_obs;
  
  real<lower=0,upper=1> rho_obs;

  
 }

transformed parameters{
  matrix[n_ind, J] par;
  vector[n_ind] t0p;
  matrix[2,2] varcov;
  cholesky_factor_cov[2] Lobs;

  for(i in 1:n_ind){
    par[i] = (mu + diag_pre_multiply(sigma,L)*z[i]')';
    t0p[i] = exp(mu_t + z_t[i]*sigma_t);
  }

  varcov[1,1] = sigma_obs[1]^2;
  varcov[2,2] = sigma_obs[2]^2;
  varcov[1,2] = sigma_obs[1]*sigma_obs[2]*rho_obs;
  varcov[2,1] = varcov[1,2];
  Lobs = cholesky_decompose(varcov);

 }

model {

 matrix[N_b,2] obs_mean;

 for(i in 1:n_ind){
      // par[i] ~ multi_normal_cholesky(mu,diag_pre_multiply(sigma,L));
      z[i] ~ std_normal();     
      // t0p[i] ~ lognormal(mu_t, sigma_t);
      z_t[i] ~ std_normal();
      }

 for (i in 1:N_y){
      y_y[i] ~ normal(par[ind_y[i],1]*(1-inv_logit(par[ind_y[i],2])^(age_y[i] + t0p[ind_y[i]])), sigma_obs[1]);
      }
      
 for (i in 1:N_z) {
      z_z[i] ~ normal(par[ind_z[i],3]*(1-inv_logit(par[ind_z[i],4])^(age_z[i] + t0p[ind_z[i]])), sigma_obs[2]);
      } 
    
 for (i in 1:N_b){
      obs_mean[i,1] = par[ind_b[i],1]*(1-inv_logit(par[ind_b[i],2])^(age_b[i] + t0p[ind_b[i]]));
      obs_mean[i,2] = par[ind_b[i],3]*(1-inv_logit(par[ind_b[i],4])^(age_b[i] + t0p[ind_b[i]]));
      data_b[i] ~ multi_normal_cholesky(obs_mean[i], Lobs);

      }


//priors
  //t0p ~ lognormal(0, 10);
  mu_t ~ normal(0, 1000);
  sigma_t ~ student_t(3, 0, 25);
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

