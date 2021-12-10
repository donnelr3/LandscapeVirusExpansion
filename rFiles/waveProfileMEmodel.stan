// Description: Bayesian regression model of landscape wave-profiles, A/(I+1), by transect distance
// distance FE
// distance^2 FE
// year specific distance effects 
// option of shared intercept across years
// option of random effect on field id#s if there are multiple years (eg. Landscape Experiment Data vs. Landscape Survey Data)
// coded for 1 or two years data only

// note convention throughout is that variables with D_ in front all involve data coming as input
data {                                    // Data block
  int<lower=1> D_Nfields;      // number of fields in transect
  int<lower=1> D_levels;        // Number of years data
  int<lower=1> D_K;              // Number of coeficents per level eg. 3 for quadratic regression
  matrix[D_K,D_Nfields] D_X;                            // Model Matrix 
  matrix[D_levels, D_Nfields] D_y;                        // Response variable (wave-profile)
  int<lower=1, upper=D_Nfields> D_subj[D_Nfields*D_levels]; //subject id (for R.E.)
  int<lower=1> D_intcpt;   // shared intercept over years has value 1; value 2 has separate intercepts
}

parameters {                             // Parameters block
  vector[D_K] beta1; // FE vector level 1
  real beta2[(D_levels>1) ? ((D_intcpt>1) ? (D_K) : (D_K-1)) : 0];  // cond. declar. FE vector level 2 depends on
  real<lower=0> sigma_eps;   // Res sd
  real<lower=0> sigma_id[(D_levels>1) ? 1 : 0];    // RE sd
  real eta[(D_levels>1) ? D_Nfields : 0];    // RE sd
}
 
transformed parameters {                // Trans-Parameters block
  
}

model {                     // Model block
  vector[D_Nfields] ran1;
  vector[D_Nfields*D_levels] yhat;
  real shared_beta2[((D_levels>1)&&(D_intcpt==1)) ? (D_K) : 0];  // conditional dummy vector taking shared intercept yr2
  matrix[D_levels,D_K] beta; // coefficient matrix 
  if (D_levels>1){
    ran1=sigma_id[1] * to_vector(eta);
  }else{
    ran1=rep_vector(0.0, D_Nfields);
  }
  
  if ((D_levels>1)&&(D_intcpt==1)) {// accounting 4 awkard multi-yr (shared vs separ. intercept) vs single yr
    shared_beta2[1]=beta1[1];
    shared_beta2[2:D_K]=beta2;
    beta = [beta1',to_vector(shared_beta2)']; 
  }else if ((D_levels>1)&&(D_intcpt>1)){
    beta = [beta1',to_vector(beta2)']; 
  }else{
    beta = to_matrix(beta1,1,D_K);
  }

  for (i in 1:(D_Nfields*D_levels))
  yhat[i]=to_vector(beta*D_X)[i]+ran1[D_subj[i]];




  if (D_levels>1) { 
  to_vector(eta)~normal(0.0, 1);// individual field variation
  sigma_id[1] ~ cauchy(0.0, 1);// sigma_id bounded at 0, so half-cauchy
  }
  sigma_eps ~ cauchy(0.0, 1);     // sigma_eps bounded at 0, so half-cauchy
  to_vector(beta1) ~ normal(0.0, 1);  // standard normal priors on FE coefficients
  if (D_levels>1) {
  to_vector(beta2) ~ normal(0.0, 1);
  }
  to_vector(D_y) ~ normal(yhat, sigma_eps);
}

generated quantities {      // Generated quantities block. 

  vector[D_levels] compoundPam; // want to output posteriors for response turning point wrt distance

  compoundPam[1]=-beta1[2]/(2*beta1[3]);
  if (D_levels>1) {
    compoundPam[2]=-beta2[1]/(2*beta2[2]); 
  }
}
