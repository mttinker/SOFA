// SOFAfit
// This Stan code fits a model of foraging behavior to observational data 
//  on sea otter foraging  [Refer to "SOFA_summary" for details] 
//  
// NOTE: this version does not incorporate sub-groups...
// alternate version will incorporate user-defined sub-groups
// (could be individual otters and/or temporal or spatial groupings)
// by estimating separate eta/alpha simplex vectors for each group,
// and making relevant forage params hierarchical (psi1, phi1, muSZ)
// to allow for differences in group-specific HT, CR, & SZ values
//   (and thus also the group-specific PD and ER values)
//
// Section 1. Data inputs for model
data {
  int<lower=1> Nbouts;            // Number of feeding bouts (for effort allocation observations)
  int<lower=1> K;                 // Number of prey types (1 to (K-1) for ID prey, UN-ID prey = K)
  int<lower=1> Km1;               // "K - 1", number of prey types excluding UN-ID prey  
  int<lower=0> EffortP[Nbouts,K]; // Minutes allocated to each prey type (incl. un-id) by bout  
  int<lower=1> NSz;               // Number of prey Size obs
  int<lower=1> NHt;               // Number of prey Handling time (HT) obs  
  int<lower=1> NCR;               // Number of prey Consumption rate obs
  int<lower=1> NU[2];             // Number of un-ID prey obs (2 different param types, FOR NOW)
  real<lower=0> SZmnU[NU[1]];     // Array of mean prey size vaues (cm) for UN-ID prey, by bout  
  real<lower=0> HTmnU[NU[2]];     // Array of mean handling times (sec), for UN-ID prey, by bout  
  int<lower=1,upper=Km1> Sp[NSz]; // Array of prey id values for size obs
  int<lower=1,upper=Km1> Hp[NHt]; // Array of prey id values for HT obs  
  int<lower=1,upper=Km1> Cp[NCR]; // Array of prey id values for Cons Rate obs
  real<lower=0> SZmn[NSz];        // Array of mean prey size vaues (cm) by prey, by bout  
  real<lower=0> HTmn[NHt];        // Array of mean handling times (sec), by prey, by bout   
  real<lower=0> CRate[NCR];       // Array of prey-specific consumption rates (CR) by bout (g/min)
  vector[NCR] Csz;                // Array of log(prey size), covariate for CR obs
  vector[NHt] Hsz;                // Array of log(prey size), covariate for HT obs
  vector[NU[2]] Hsz_u;            // Array of log(prey size), covariate for UN-ID HT obs
  vector<lower=0>[Km1] Cal_dns_mn;// Vector of caloric density values (kcal/g) by prey type, mean
  vector<lower=0>[Km1] Cal_dns_sg;// Vector of caloric density values (kcal/g) by prey type, std derr
}
// Section 2. The parameters to be estimated 
parameters {
  simplex[Km1] eta ;             // mean proportional effort allocation by prey type (if all ID'd)
  // simplex [Km1] q ;               // mean proportional contribution of each prey type to Un-ID prey
  // real<lower=0,upper=1> PUID ;   // mean TOTAL Proportion of effort on Un-ID prey  
  vector<lower=0,upper=1>[Km1] Pid ;// prey-specific probability of positive ID
  real<lower=0> tau ;            // relative precision of diet composition across bouts  
  simplex[K] theta[Nbouts] ;     // vectors of bout-specific prey-allocation probs for multinomial  
  vector[Km1] muSZ ;             // mean of log(Size_cm) of each prey
  //real<lower=0> muSZ_u ;         // mean of log(Size_cm) of UN-ID
  vector<lower=0>[Km1] sigSZ ;   // stdev of log(Size_cm) of each prey (incl. un-id)
  vector<lower=0>[Km1] sigHT ;   // stdev of log(HT per item) of each prey (incl. un-id)
  vector<lower=0>[Km1] sigCR ;   // stdev of log(HT per item) of each prey (incl. un-id)
  // real<lower=0> sigSZ_u ;        // stdev of log(Size_cm) of UN-ID
  // real<lower=0> sigHT_u ;        // stdev of log(HT per item) of UN-ID    
  vector[Km1] phi1 ;             // intercept of log(mean cons rate) fxn by prey 
  vector[Km1] phi2 ;             // size effect on log(mean cons rate) fxn by prey 
  vector[Km1] psi1 ;             // intercept of log(HT) fxn by prey 
  vector[Km1] psi2 ;             // size effect on log(HT) fxn by prey   
  // real psi1_u ;                  // intercept of log(HT) fxn for UN-ID prey 
  // real psi2_u ;                  // size effect on log(HT) fxn for UN-ID prey     
  // real<lower=0,upper=1> maxPunid;// max prob un-id if size/HT dist.s overlap 100% with UN-ID
  vector<lower=0>[Km1] Cal_dens; // Caloric density (kcal/g) by prey type
}
// Section 3. Derived (transformed) parameters
transformed parameters {
  // vector<lower=0,upper=1>[Km1] q ;// prey-specific contribution to UN-ID prey
  simplex[Km1] q ;               // mean proportional contribution of each prey type to Un-ID prey
  // vector<lower=0,upper=1>[Km1] Pid ;// prey-specific probability of positive ID
  // vector<lower=0>[Km1] alphaP ;     // vector of dirichlet params: relative prey freq (w/o UN-ID)
  vector<lower=0>[K] alpha ;        // vector of dirichlet params: relative prey freq (inc. UN-ID)
  // real<lower=0> OV_SZ[Km1] ;    
  // real<lower=0> OV_HT[Km1] ;
  // real<lower=0> Dist_SZ[Km1] ;
  // real<lower=0> Dist_HT[Km1] ;   
  // vector[Km1] muHT ;             // mean of log(HT per item) of each prey  
  //real muHT_u ;         // mean of log(HT per item) of UN-ID     
  // Compute mean expected HT at the mean size of un-id prey:
  // muHT = psi1 + psi2 .* muSZ ;
  //muHT_u = psi1_u + psi2_u * muSZ_u ;  
  // Calculate alpha, the dirichlet params for bout-specific prey allocation probs  
  // Pid = 1 - ( (q * PUID) ./ eta );
  // alphaP = eta * (tau * Km1) ;
  alpha[1:Km1] = (eta * (tau * Km1)) .* Pid ;
  alpha[K] = sum((eta * (tau * Km1)) .* (1 - Pid)) ;
  // Calculate proporitonal contribution to un-ID 
  q = (eta .* (1 - Pid)) / sum(eta .* (1 - Pid)) ;
}
// Section 4. Estimating model parameters (drawing from probability distributions)
model {
  // A) Observed nodes:
  // Allocation of effort by prey type (proportional # minutes of each feeding bout)
  for (i in 1:Nbouts){
    theta[i] ~ dirichlet(alpha);  
    EffortP[i,] ~ multinomial(theta[i]) ;
  }
  // Observed consumption rate (g/min) by prey type, random samples
  CRate ~ lognormal(phi1[Cp] + phi2[Cp] .* Csz, sigCR[Cp]) ;
  // Mean prey size (cm), random samples, ID prey 
  SZmn ~ lognormal(muSZ[Sp],sigSZ[Sp]) ;
  // Handling time (sec), random samples, ID prey
  HTmn ~ lognormal(psi1[Hp] + psi2[Hp] .* Hsz, sigHT[Hp]) ;
  // Mean prey size (cm), random samples, UN-id prey
  //SZmnU ~ lognormal(muSZ_u,sigSZ_u) ;
  for (i in 1:NU[1]){
    real lmix[Km1];
    for(j in 1:Km1){
      lmix[j] = log(q[j]) + normal_lpdf(SZmnU[i] | muSZ[j], sigSZ[j]) ; 
    }
    target += log_sum_exp(lmix);
  }
  // Handling time (sec), random samples, UN-id prey
  //HTmnU ~ lognormal(psi1_u + psi2_u * Hsz_u[i], sigHT_u) ;
  for (i in 1:NU[2]){
    real lmix[Km1];
    for(j in 1:Km1){
      lmix[j] = log(q[j]) + normal_lpdf(HTmnU[i] |  psi1[j] + psi2[j] * Hsz_u[i], sigHT[j]) ; 
    }
    target += log_sum_exp(lmix);
  }
  //
  // B) Prior distributions for model parameters:
  Cal_dens ~ normal(Cal_dns_mn,Cal_dns_sg) ; // Caloric deinsity (allows uncertainty)
  // NOTE: also include normal prior for CR_uncert multiplier, with mean of 1,
  //   to incorporate uncertainty associated with Mass-size power functions
  tau ~ cauchy(0,2.5) ;  // precision param for diet comp variaiton across bouts
  phi1 ~ cauchy(0,2.5) ; // Intercept of log Cons Rate, by prey
  phi2 ~ cauchy(0,2.5) ; // effect of size on log Cons rate
  psi1 ~ cauchy(0,2.5) ; // Intercept of log HT, by prey
  psi2 ~ cauchy(0,2.5) ; // effect of size on log HT rate
  // psi1_u ~ cauchy(0,2.5) ; // Intercept of log HT, UN-ID prey
  // psi2_u ~ cauchy(0,2.5) ; // effect of size on log HTm UN-ID prey
  muSZ ~ cauchy(1.5,2.5);// log-mean size by prey type
  // muSZ_u ~ cauchy(1.5,2.5);// log-mean size by prey type
  // PUID ~ beta(1,10);     // mean TOTAL proportion of effort that is Un-ID
  Pid ~ beta(5,1);     // max prob un-ID, for prey overlapping UN-ID in size & HT
  sigSZ ~ cauchy(0,2.5); // variation in prey size across bouts, by prey type
  sigHT ~ cauchy(0,2.5); // variation in prey HT across bouts, by prey type
  sigCR ~ cauchy(0,2.5); // variation in prey HT across bouts, by prey type
  // sigSZ_u ~ cauchy(0,2.5); 
  // sigHT_u ~ cauchy(0,2.5); 
}
// Section 5. Derived parameters and statistics 
generated quantities {
  vector[Km1] SZ ;
  // vector[Km1] NI ;
  vector[Km1] HT ;
  // vector[Km1] LM ;
  vector[Km1] CR ;
  vector[Km1] PD ;
  vector[Km1] ER ;
  real CRmn ;
  real ERmn ;
  // Mean Size (cm) and HT/itm, by prey type (adjust for log-normal)
  SZ = exp(muSZ +  square(sigSZ)/2 );
  HT = exp( (psi1 + psi2 .* muSZ) +  square(sigHT)/2 );
  // Mean Cons rate (CR, g/min) by prey type, adjusted for mean prey size 
  CR = exp( (phi1 + phi2 .* muSZ) +  square(sigCR)/2 ) ;
  // Mean Energy intake rate (ER) by prey, incl. uncertainty in Caloric density
  ER = Cal_dens .* CR ;
  // Proportional contribution (biomass consumed) of each prey type to diet: 
  PD = (eta .* CR) / sum(eta .* CR) ;
  // Overall mean consumption rate (CRmn) given effort allocation to each prey: 
  CRmn = sum(eta .* CR) ;
  // Overall mean Energy Intake Rate (ERmn, kcal.min) given effort allocation: 
  ERmn = sum(eta .* CR .* Cal_dens) ;
  // 
}
