
data {
  int<lower=0> N_unknown;      // Number of unknown status samples
  real Ab_unknown[N_unknown];  // Ab values for unknown samples
  real Ab_min;                 // maximum value of Ab (for priors constraint)
  real Ab_max;                 // minimum value of Ab (for priors constraint)
}


// The parameters accepted by the model. 
parameters {
  real<lower=0, upper = 1> theta;                  // seroprevalence
  real<lower=Ab_min, upper = Ab_max> mu_Neg;       // Mean AB value for Negatives
  real<lower=0, upper = 3> sigma_Neg;              // Std deviation of mu_N     
  real<lower=log(2), upper = Ab_max> mu_delta;     // delta between mu_N and mu_P, has to be positive to ensure mu_P greater than mu_N
  real<lower=0, upper = 3> sigma_Pos;              // Std deviation of mu_P
}


transformed parameters {
    // Serocatalytic model
    real mu_Pos;
    
    mu_Pos = mu_Neg + mu_delta;  // Mean of the positives
}



// The model to be estimated. 
model {
  // Priors definition
  theta     ~ uniform(0,1);
  mu_Neg    ~ uniform(Ab_min,Ab_max);
  sigma_Neg ~ uniform(0,100);
  mu_delta  ~ uniform(log(2), Ab_max-log(2));  
  sigma_Pos ~ uniform(0,100);


// likelihood for the unknown
  for (l in 1:N_unknown){
    target += log_mix(theta,
                      normal_lpdf(Ab_unknown[l] | mu_Pos, sigma_Pos),
                      normal_lpdf(Ab_unknown[l] | mu_Neg, sigma_Neg));
  }
}
