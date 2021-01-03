data {
  int<lower=0> Number_Clusters;
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
  matrix[Number_Envt_Variables, Number_Clusters] envt_betas;
  vector[Number_Clusters] alpha_ann;
  vector[Number_Clusters] alpha_cul;
  vector[Number_Clusters] alpha_dir;
  vector[Number_Clusters] alpha_fluv;
  vector[Number_Clusters] alpha_min;
  vector[Number_Clusters] alpha_ste;
  vector[Number_Clusters] alpha_sub;
}

model {
  matrix[Number_Obs_Ann, Number_Clusters] beta_x_envt_ann = Envt_Variables_Ann * envt_betas; 
  matrix[Number_Obs_Cul, Number_Clusters] beta_x_envt_cul = Envt_Variables_Cul * envt_betas; 
  matrix[Number_Obs_Dir, Number_Clusters] beta_x_envt_dir = Envt_Variables_Dir * envt_betas; 
  matrix[Number_Obs_Fluv, Number_Clusters] beta_x_envt_fluv = Envt_Variables_Fluv * envt_betas; 
  matrix[Number_Obs_Min, Number_Clusters] beta_x_envt_min = Envt_Variables_Min * envt_betas; 
  matrix[Number_Obs_Ste, Number_Clusters] beta_x_envt_ste = Envt_Variables_Ste * envt_betas; 
  matrix[Number_Obs_Sub, Number_Clusters] beta_x_envt_sub = Envt_Variables_Sub * envt_betas; 

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
    Cluster_Identity_Ann[a] ~ categorical_logit(beta_x_envt_ann[a]' + alpha_ann);
  }
  for (b in 1:Number_Obs_Cul) {
    Cluster_Identity_Cul[b] ~ categorical_logit(beta_x_envt_cul[b]' + alpha_cul);
  }
  for (c in 1:Number_Obs_Dir) {
    Cluster_Identity_Dir[c] ~ categorical_logit(beta_x_envt_dir[c]' + alpha_dir);
  }
  for (d in 1:Number_Obs_Fluv) {
    Cluster_Identity_Fluv[d] ~ categorical_logit(beta_x_envt_fluv[d]' + alpha_fluv);
  }
  for (e in 1:Number_Obs_Min) {
    Cluster_Identity_Min[e] ~ categorical_logit(beta_x_envt_min[e]' + alpha_min);
  }
  for (f in 1:Number_Obs_Ste) {
    Cluster_Identity_Ste[f] ~ categorical_logit(beta_x_envt_ste[f]' + alpha_ste);
  }
  for (g in 1:Number_Obs_Sub) {
    Cluster_Identity_Sub[g] ~ categorical_logit(beta_x_envt_sub[g]' + alpha_sub);
  }
}

generated quantities {
  vector[Number_Clusters] temp_storage_ann;
  vector[Number_Clusters] temp_storage_cul;
  vector[Number_Clusters] temp_storage_dir;
  vector[Number_Clusters] temp_storage_fluv;
  vector[Number_Clusters] temp_storage_min;
  vector[Number_Clusters] temp_storage_ste;
  vector[Number_Clusters] temp_storage_sub;
  
  matrix[Number_Obs_Ann, Number_Clusters] Cluster_Probabilities_Ann;
  matrix[Number_Obs_Cul, Number_Clusters] Cluster_Probabilities_Cul;
  matrix[Number_Obs_Dir, Number_Clusters] Cluster_Probabilities_Dir;
  matrix[Number_Obs_Fluv, Number_Clusters] Cluster_Probabilities_Fluv;
  matrix[Number_Obs_Min, Number_Clusters] Cluster_Probabilities_Min;
  matrix[Number_Obs_Ste, Number_Clusters] Cluster_Probabilities_Ste;
  matrix[Number_Obs_Sub, Number_Clusters] Cluster_Probabilities_Sub;

  matrix[Number_Obs_Ann, Number_Clusters] beta_x_envt_ann = Envt_Variables_Ann * envt_betas; 
  matrix[Number_Obs_Cul, Number_Clusters] beta_x_envt_cul = Envt_Variables_Cul * envt_betas; 
  matrix[Number_Obs_Dir, Number_Clusters] beta_x_envt_dir = Envt_Variables_Dir * envt_betas; 
  matrix[Number_Obs_Fluv, Number_Clusters] beta_x_envt_fluv = Envt_Variables_Fluv * envt_betas; 
  matrix[Number_Obs_Min, Number_Clusters] beta_x_envt_min = Envt_Variables_Min * envt_betas; 
  matrix[Number_Obs_Ste, Number_Clusters] beta_x_envt_ste = Envt_Variables_Ste * envt_betas; 
  matrix[Number_Obs_Sub, Number_Clusters] beta_x_envt_sub = Envt_Variables_Sub * envt_betas; 
  
  for (a in 1:Number_Obs_Ann) {
    temp_storage_ann = softmax(beta_x_envt_ann[a]' + alpha_ann);
    Cluster_Probabilities_Ann[a] = temp_storage_ann';
  }
  for (b in 1:Number_Obs_Cul) {
    temp_storage_cul = softmax(beta_x_envt_cul[b]' + alpha_cul);
    Cluster_Probabilities_Cul[b] = temp_storage_cul';
  }
  for (c in 1:Number_Obs_Dir) {
    temp_storage_dir = softmax(beta_x_envt_dir[c]' + alpha_dir);
    Cluster_Probabilities_Dir[c] = temp_storage_dir';
  }
  for (d in 1:Number_Obs_Fluv) {
    temp_storage_fluv = softmax(beta_x_envt_fluv[d]' + alpha_fluv);
    Cluster_Probabilities_Fluv[d] = temp_storage_fluv';
  }
  for (e in 1:Number_Obs_Min) {
    temp_storage_min = softmax(beta_x_envt_min[e]' + alpha_min);
    Cluster_Probabilities_Min[e] = temp_storage_min';
  }
  for (f in 1:Number_Obs_Ste) {
    temp_storage_ste = softmax(beta_x_envt_ste[f]' + alpha_ste);
    Cluster_Probabilities_Ste[f] = temp_storage_ste';
  }
  for (g in 1:Number_Obs_Sub) {
    temp_storage_sub = softmax(beta_x_envt_sub[g]' + alpha_sub);
    Cluster_Probabilities_Sub[g] = temp_storage_sub';
  }
}
