data {
  int<lower=0> Number_Envt_Variables;
  int<lower=0> Number_Obs_Ann;
  int<lower=0> Number_Obs_Cul;
  int<lower=0> Number_Obs_Dir;
  int<lower=0> Number_Obs_Fluv;
  int<lower=0> Number_Obs_Min;
  int<lower=0> Number_Obs_Ste;
  int<lower=0> Number_Obs_Sub;

  matrix[Number_Obs_Ann, Number_Envt_Variables] Envt_Variables_Ann;
  matrix[Number_Obs_Cul, Number_Envt_Variables] Envt_Variables_Cul;
  matrix[Number_Obs_Dir, Number_Envt_Variables] Envt_Variables_Dir;
  matrix[Number_Obs_Fluv, Number_Envt_Variables] Envt_Variables_Fluv;
  matrix[Number_Obs_Min, Number_Envt_Variables] Envt_Variables_Min;
  matrix[Number_Obs_Ste, Number_Envt_Variables] Envt_Variables_Ste;
  matrix[Number_Obs_Sub, Number_Envt_Variables] Envt_Variables_Sub;

  int Cluster_Identity_Ann[Number_Obs_Ann];
  int Cluster_Identity_Cul[Number_Obs_Cul];
  int Cluster_Identity_Dir[Number_Obs_Dir];
  int Cluster_Identity_Fluv[Number_Obs_Fluv];
  int Cluster_Identity_Min[Number_Obs_Min];
  int Cluster_Identity_Ste[Number_Obs_Ste];
  int Cluster_Identity_Sub[Number_Obs_Sub];
}

parameters {
  real<lower=0> tau;
  vector[Number_Envt_Variables] envt_betas;
  real alpha_ann;
  real alpha_cul;
  real alpha_dir;
  real alpha_fluv;
  real alpha_min;
  real alpha_ste;
  real alpha_sub;
}

model {
  vector[Number_Obs_Ann] beta_x_envt_ann = Envt_Variables_Ann * envt_betas; 
  vector[Number_Obs_Cul] beta_x_envt_cul = Envt_Variables_Cul * envt_betas; 
  vector[Number_Obs_Dir] beta_x_envt_dir = Envt_Variables_Dir * envt_betas; 
  vector[Number_Obs_Fluv] beta_x_envt_fluv = Envt_Variables_Fluv * envt_betas; 
  vector[Number_Obs_Min] beta_x_envt_min = Envt_Variables_Min * envt_betas; 
  vector[Number_Obs_Ste] beta_x_envt_ste = Envt_Variables_Ste * envt_betas; 
  vector[Number_Obs_Sub] beta_x_envt_sub = Envt_Variables_Sub * envt_betas; 

  tau ~ normal(0, 150);
  to_vector(envt_betas) ~ normal(0, tau);
  alpha_ann ~ normal(0, tau);
  alpha_cul ~ normal(0, tau);
  alpha_dir ~ normal(0, tau);
  alpha_fluv ~ normal(0, tau);
  alpha_min ~ normal(0, tau);
  alpha_ste ~ normal(0, tau);
  alpha_sub ~ normal(0, tau);

  for (a in 1:Number_Obs_Ann) {
    Cluster_Identity_Ann[a] ~ bernoulli_logit(beta_x_envt_ann[a]' + alpha_ann);
  }
  for (b in 1:Number_Obs_Cul) {
    Cluster_Identity_Cul[b] ~ bernoulli_logit(beta_x_envt_cul[b]' + alpha_cul);
  }
  for (c in 1:Number_Obs_Dir) {
    Cluster_Identity_Dir[c] ~ bernoulli_logit(beta_x_envt_dir[c]' + alpha_dir);
  }
  for (d in 1:Number_Obs_Fluv) {
    Cluster_Identity_Fluv[d] ~ bernoulli_logit(beta_x_envt_fluv[d]' + alpha_fluv);
  }
  for (e in 1:Number_Obs_Min) {
    Cluster_Identity_Min[e] ~ bernoulli_logit(beta_x_envt_min[e]' + alpha_min);
  }
  for (f in 1:Number_Obs_Ste) {
    Cluster_Identity_Ste[f] ~ bernoulli_logit(beta_x_envt_ste[f]' + alpha_ste);
  }
  for (g in 1:Number_Obs_Sub) {
    Cluster_Identity_Sub[g] ~ bernoulli_logit(beta_x_envt_sub[g]' + alpha_sub);
  }
}

generated quantities {

  vector[Number_Obs_Ann] Cluster_Probabilities_Ann;
  vector[Number_Obs_Cul] Cluster_Probabilities_Cul;
  vector[Number_Obs_Dir] Cluster_Probabilities_Dir;
  vector[Number_Obs_Fluv] Cluster_Probabilities_Fluv;
  vector[Number_Obs_Min] Cluster_Probabilities_Min;
  vector[Number_Obs_Ste] Cluster_Probabilities_Ste;
  vector[Number_Obs_Sub] Cluster_Probabilities_Sub;

  vector[Number_Obs_Ann] beta_x_envt_ann = Envt_Variables_Ann * envt_betas; 
  vector[Number_Obs_Cul] beta_x_envt_cul = Envt_Variables_Cul * envt_betas; 
  vector[Number_Obs_Dir] beta_x_envt_dir = Envt_Variables_Dir * envt_betas; 
  vector[Number_Obs_Fluv] beta_x_envt_fluv = Envt_Variables_Fluv * envt_betas; 
  vector[Number_Obs_Min] beta_x_envt_min = Envt_Variables_Min * envt_betas; 
  vector[Number_Obs_Ste] beta_x_envt_ste = Envt_Variables_Ste * envt_betas; 
  vector[Number_Obs_Sub] beta_x_envt_sub = Envt_Variables_Sub * envt_betas; 
  
  for (a in 1:Number_Obs_Ann) {
    Cluster_Probabilities_Ann[a] = inv_logit(beta_x_envt_ann[a] + alpha_ann);
  }
  for (b in 1:Number_Obs_Cul) {
    Cluster_Probabilities_Cul[b] = inv_logit(beta_x_envt_cul[b] + alpha_cul);
  }
  for (c in 1:Number_Obs_Dir) {
    Cluster_Probabilities_Dir[c] = inv_logit(beta_x_envt_dir[c] + alpha_dir);
  }
  for (d in 1:Number_Obs_Fluv) {
    Cluster_Probabilities_Fluv[d] = inv_logit(beta_x_envt_fluv[d] + alpha_fluv);
  }
  for (e in 1:Number_Obs_Min) {
    Cluster_Probabilities_Min[e] = inv_logit(beta_x_envt_min[e] + alpha_min);
  }
  for (f in 1:Number_Obs_Ste) {
    Cluster_Probabilities_Ste[f] = inv_logit(beta_x_envt_ste[f] + alpha_ste);
  }
  for (g in 1:Number_Obs_Sub) {
    Cluster_Probabilities_Sub[g] = inv_logit(beta_x_envt_sub[g] + alpha_sub);
  }
}
