data {
  // Datapoints - Recorded and Those to Be Predicted
  int<lower=1> N_train;                   // Number of actual datapoints
  int<lower=1> N_pred;              // Number of datapoints to be predicted
  int y[N_train];                         // Actual datapoint, output value
  real x_train[N_train];                        // Actual datapoint, index position
  real x_pred[N_pred];              // Datapoint to be predicted, index position
  
  // Model Parameters
	real<lower=0> length_scale_mean;  // Length scale prior mean and sd
	real<lower=0> length_scale_sd;    
	real<lower=0> period_mean;        // Period prior mean and sd
	real<lower=0> period_sd;          
	real<lower=0> alpha_mean;         // Alpha prior mean and sd (controlling variability of GP) 
  real<lower=0> alpha_sd;  
  real<lower=0> overdispersion_mean;  // Sd for negative binomial overdispersion parameter prior
  real<lower=0> overdispersion_sd;  // Sd for negative binomial overdispersion parameter prior
}

transformed data {
  
  int<lower=1> N_tot;     // Length of vector containing both actual datapoints and those to be predicted
  int<lower=1> k;         // A Counter
  real x_tot[N_pred + N_train]; // Vector storing the index positions for both actual datapoints and those to be predicted
  
  N_tot = N_pred + N_train;     
  k = 1;
  
  for (n in 1:N_train) {        // Filling the vector x_tot with both the actual datapoints' index position, 
    x_tot[k] = x_train[n];      // as well as the index position for the datapoints which we're hoping to predict
    k = k + 1;
  }
  for (n in 1:N_pred) {   
    x_tot[k] = x_pred[n];
    k = k + 1;
  }
}

parameters {
  real<lower=1, upper=12> length_scale;
  real<lower=4, upper=18> period;
  real<lower=0> alpha;
  vector[N_tot] eta;
  real<lower=0> overdispersion;
}

transformed parameters {
  vector[N_tot] f;        
  {                                 // think this stops this stuff from being stored by STAN
    matrix[N_tot, N_tot] L;
    matrix[N_tot, N_tot] Sigma;
	  for (i in 1:N_tot) {   
  		for (j in 1:i) {
  			if (i == j) 
  			  Sigma[i,j] = pow(alpha, 2) * exp(-(2/pow(length_scale, 2)) * pow((sin((pi() * (x_tot[i] - x_tot[j]))/period)), 2));
  			else 
  			  Sigma[i,j] = pow(alpha, 2) * exp(-(2/pow(length_scale, 2)) * pow((sin((pi() * (x_tot[i] - x_tot[j]))/period)), 2));
    		Sigma[j,i] = Sigma[i,j];
  		}
		  Sigma[i,i] = Sigma[i,i] + 1e-12;
	  }
    L = cholesky_decompose(Sigma);  // cholesky decomposition of the calculated covariance matrix
    f = L * eta;                    // choleksy decompose * normal random variate is equivalent to sampling from the required MVN
  }
}

// Model Block 
model{
	  // Priors 
	  length_scale ~ normal(length_scale_mean, length_scale_sd);
	  period ~ normal(period_mean, period_sd);
	  alpha ~ normal(alpha_mean, alpha_sd);
    overdispersion ~ normal(overdispersion_mean, overdispersion_sd);
    
    // Sample from Normal for MVN Sampling
	  eta ~ normal(0, 1);
	  
	  // Calculating the Likelihood
	  y ~ neg_binomial_2(exp(f[1:N_train]), overdispersion); 
}
