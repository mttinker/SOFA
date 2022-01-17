// SOFAfit
// This Stan code fits a model of foraging behavior to observational data 
//  on sea otter foraging  [Refer to "SOFA_summary" for details] 
//  - this is simple version for "Non-grouped" data
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
  int<lower=1> NLm;               // Number of dive success rate obs (lambda)
  int<lower=1> NU[2];             // Number of un-ID prey obs (2 different param types, FOR NOW)
  real<lower=0> SZmnU[NU[1]];     // Array of mean prey size vaues (mm) for UN-ID prey, by bout  
  real<lower=0> HTmnU[NU[2]];     // Array of mean handling times (sec), for UN-ID prey, by bout  
  int<lower=1,upper=Km1> Sp[NSz]; // Array of prey id values for size obs
  int<lower=1,upper=Km1> Hp[NHt]; // Array of prey id values for HT obs  
  int<lower=1,upper=Km1> Cp[NCR]; // Array of prey id values for Cons Rate obs
  real<lower=0> SZmn[NSz];        // Array of mean prey size vaues (mm) by prey, by bout  
  real<lower=0> HTmn[NHt];        // Array of mean handling times (sec), by prey, by bout   
  real<lower=0> CRate[NCR];       // Array of prey-specific consumption rates (CR) by bout (g/min)
  vector[NCR] Csz;                // Array of prey size (log pwr transform), covariate for CR obs
  vector[NCR] Css;                // Array of sample sizes for CR obs (N obs per bout)
  vector[NSz] Sss;                // Array of sample sizes (n obs per bout) for size observations 
  vector[NU[1]] Sss_u;            // Array of sample sizes (n obs per bout) for size observations, Un-ID   
  vector[NHt] Hsz;                // Array of prey size (log pwr transform), covariate for HT obs
  vector[NHt] Hss;                // Array of sample sizes for HT obs (n obs per bout)
  vector[NU[2]] Hss_u;            // Array of sample sizes for HT obs (n obs per bout), UN-ID
  vector[NU[2]] Hsz_u;            // Array of prey size (log pwr transform), covariate for UN-ID HT obs
  vector<lower=0>[Km1] Cal_dns_mn;// Vector of caloric density values (kcal/g) by prey type, mean
  vector<lower=0>[Km1] Cal_dns_sg;// Vector of caloric density values (kcal/g) by prey type, std err
  vector<lower=0>[Km1] logMass_sg;// Vector of std err vlues for log-log fxn of bomass vs size
  int<lower=1,upper=Km1> Lp[NLm]; // Array of prey id values for success rate obs
  real LMlg[NLm];                 // Array of logit(lambda), proportion of succesful dives by prey, by bout  
  vector[NLm] Lss;                // Array of sample sizes (n obs per bout) for logit lambda observations 
}
// Section 2. The parameters to be estimated 
parameters {
  simplex[Km1] eta ;             // mean proportional effort allocation by prey type (if all ID'd)
  real<lower=0> tauB ;           // relative precision of diet composition across bouts
  simplex[K] theta[Nbouts] ;     // vectors of bout-specific prey-allocation probs for multinomial  
  vector<lower=0>[Km1] muSZ ;    // mean of log(Size_mm) of each prey
  real<lower=0> muSZ_u ;         // mean of log(Size_mm) of UN-ID
  vector[Km1] lgtLM ;            // mean logit(lambda) (dive success rate) for each prey
  vector<lower=0>[Km1] sigSZ ;   // stdev of log(Size_cm) of each prey 
  vector<lower=0>[Km1] sigHT ;   // stdev of log(HT per item) of each prey 
  vector<lower=0>[Km1] sigCR ;   // stdev of log(HT per item) of each prey 
  vector<lower=0>[Km1] sigLM ;   // stdev of logit(lambda) of each prey 
  real<lower=0> sigSZ_u ;        // stdev of log(Size_cm) of UN-ID
  real<lower=0> sigHT_u ;        // stdev of log(HT per item) of UN-ID    
  vector<lower=0>[Km1] phi1 ;    // intercept of CRate v Size fxn by prey 
  vector<lower=0>[Km1] phi2 ;    // slope of CRate v Size fxn by prey 
  vector<lower=0>[Km1] psi1 ;    // intercept of HT v Size fxn by prey 
  vector<lower=0>[Km1] psi2 ;    // slope of HT v Size fxn by prey   
  real<lower=0> psi1_u ;         // intercept of HT v Size fxn for UN-ID prey 
  real<lower=0> psi2_u ;         //  slope of HT v Size fxn for UN-ID prey     
  real<lower=0,upper=1> maxPunid;// max prob un-id, if size/HT dist.s overlap 100% with UN-ID
  vector<lower=0>[Km1] Cal_dens; // Caloric density (kcal/g) by prey type
  vector[Km1] lgSz_adj; // Uncertainty adjustment for size-biomass fxn
}
// Section 3. Derived (transformed) parameters
transformed parameters {
  vector<lower=0,upper=1>[Km1] Omega ;// prey-specific probability of positive ID
  vector<lower=0>[K] alpha ;        // vector of dirichlet params: relative prey freq (inc. UN-ID)
  vector[Km1] muHT ;                // mean of log(HT per item) of each prey  
  real muHT_u ;                     // mean of log(HT per item) of UN-ID       
  //
  // Compute mean expected HT at avg size for un-id prey
  muHT_u = psi1_u + psi2_u * (2.5 * muSZ_u - 7);  
  // Loop through prey types, calculate prob of ID-ing prey (Pid) & contribution to un-ID 
  // NOTE: Bhattacharyya coefficient measures overlap between size/HT distributions of 
  // each prey type and UN-ID prey: if NO overlap then prey does not contribute to un-ID,
  // while if perfect overlap, then prob of Un-ID is at maximum (maxPunid)
  for (j in 1:Km1){
    real OV_SZ ;        // Overlap of size distribution with unknown prey, by prey type
    real OV_HT ;        // Overlap of HT distribution with unknown prey, by prey type
    real Dist_SZ ;      // Bhattacharyya distance between UnID andeach prey type, for size
    real Dist_HT ;      // Bhattacharyya distance between UnID andeach prey type, for HT
    // Compute mean expected HT at avg size for prey type j:
    muHT[j] = psi1[j] + psi2[j] * (2.5 * muSZ[j] - 7) ;  
    Dist_SZ = 0.25 * log(0.25 * (sigSZ[j]^2 / sigSZ_u^2 + sigSZ_u^2/sigSZ[j]^2 + 2))
              + 0.25 * (((muSZ[j] - muSZ_u)^2) / (sigSZ[j]^2 + sigSZ_u^2)) ;
    OV_SZ = exp(-Dist_SZ) ;
    Dist_HT = 0.25 * log(0.25 * (sigHT[j]^2 / sigHT_u^2 + sigHT_u^2/sigHT[j]^2 + 2))
              + 0.25 * (((muHT[j] - muHT_u)^2) / (sigHT[j]^2 + sigHT_u^2)) ;
    OV_HT = exp(-Dist_HT) ;    
    Omega[j] = 1 - maxPunid * (OV_SZ * OV_HT) ;
  }  
  // Calculate alpha, the dirichlet params for bout-specific prey allocation probs
  alpha[1:Km1] = eta .* Omega ;
  alpha[K] = sum(eta .* (1 - Omega)) ;
}
// Section 4. Estimating model parameters (drawing from probability distributions)
model {
  // A) Observed nodes:
  // Allocation of effort by prey type (proportional # minutes of each feeding bout)
  for (i in 1:Nbouts){
    theta[i] ~ dirichlet(alpha * tauB);  
    EffortP[i,] ~ multinomial(theta[i]) ;
  }
  // Observed consumption rate (g/min) by prey type, random samples
  // NOTE: Expected log attributes, given that the log median of means of lognormal samples of size n: 
  //           mu + (n*sg^2 - sg^2)/(2*n) 
  CRate ~ lognormal((phi1[Cp] + phi2[Cp] .* Csz) + ((Css .* square(sigCR[Cp])) - square(sigCR[Cp])) ./ (2 * Css),
          sigCR[Cp] ./ sqrt(Css) ) ;
  // Observed Handling time (sec), random samples, ID prey & UN-id prey
  HTmn ~ lognormal((psi1[Hp] + psi2[Hp] .* Hsz) + ((Hss .* square(sigHT[Hp])) - square(sigHT[Hp])) ./ (2 * Hss), 
          sigHT[Hp] ./ sqrt(Hss) ) ;
  HTmnU ~ lognormal((psi1_u + psi2_u * Hsz_u) + ((Hss_u * square(sigHT_u)) - square(sigHT_u)) ./ (2 * Hss_u), 
          sigHT_u ./ sqrt(Hss_u) ) ;  
  // Observed Mean prey size (mm), random samples, ID prey & UN-id prey
  SZmn ~ lognormal(muSZ[Sp] + ((Sss .* square(sigSZ[Sp])) - square(sigSZ[Sp])) ./ (2 * Sss), 
          sigSZ[Sp] ./ sqrt(Sss) ) ;
  SZmnU ~ lognormal(muSZ_u + ((Sss_u * square(sigSZ_u)) - square(sigSZ_u)) ./ (2 * Sss_u), 
          sigSZ_u ./ sqrt(Sss_u) ) ;
  // Observed logit(lambda), dive succes rate by bout/prey, random samples 
  LMlg ~ normal(lgtLM[Lp], sigLM[Lp] ./ sqrt(Lss) ) ;
  //
  // B) Prior distributions for model parameters:
  Cal_dens ~ normal(Cal_dns_mn,Cal_dns_sg) ; // Caloric density (incorporates est. uncertainty)
  lgSz_adj ~ normal(0,logMass_sg) ;          // CRate adjust (incorporates est. uncertainty)
  // NOTE: also include normal prior for CR_uncert multiplier, with mean of 1,
  //   to incorporate uncertainty associated with Mass-size power functions
  tauB ~ cauchy(0,2.5) ;    // precision param for diet comp variation across bouts
  phi1 ~ cauchy(0,1) ;   // Intercept of log Cons Rate, by prey
  phi2 ~ cauchy(0,1) ;   // effect of size on log Cons rate
  psi1 ~ cauchy(0,1) ;   // Intercept of log HT, by prey
  psi2 ~ cauchy(0,1) ;   // effect of size on log HT rate
  psi1_u ~ cauchy(0,1) ; // Intercept of log HT, UN-ID prey
  psi2_u ~ cauchy(0,1) ; // effect of size on log HTm UN-ID prey
  muSZ ~ normal(4,1);      // log-mean size by prey type
  muSZ_u ~ normal(4,1);    // log-mean size of UN-ID prey
  maxPunid ~ beta(3,1);    // max prob un-ID, for prey overlapping UN-ID in size & HT
  sigSZ ~ cauchy(0,2.5);   // variation in prey size across bouts, by prey type
  sigHT ~ cauchy(0,2.5);   // variation in prey HT across bouts, by prey type
  sigCR ~ cauchy(0,2.5);   // variation in prey CR across bouts, by prey type
  sigSZ_u ~ cauchy(0,2.5); // variation in prey size across bouts, unid prey
  sigHT_u ~ cauchy(0,2.5); // variation in prey HT across bouts, unid prey
  lgtLM ~ cauchy(0,2.5); // mean logit(lambda), dive success rate, by prey
  sigLM ~ cauchy(0,2.5);  // variation in logit(lambda), dive success rate
}
// Section 5. Derived parameters and statistics 
generated quantities {
  vector[Km1] SZ ;
  real SZ_u ;
  vector[Km1] HT ;
  real HT_u ;
  vector[Km1] CR ;
  vector[Km1] Pi ;
  vector[Km1] ER ;
  vector[Km1] LM ;
  real CRmn ;
  real ERmn ;
  real LMmn ;
  // Mean Size (mm) by prey type (adjust for log-normal)
  SZ = exp(muSZ +  square(sigSZ)/2 );
  SZ_u = exp(muSZ_u +  square(sigSZ_u)/2 );
  for (j in 1:Km1){
    // Mean Cons rate (CR, g/min) by prey type, adjusted for mean prey size and lognormal dist
    CR[j] = fmin(100, exp(phi1[j] + phi2[j] * (2.5*muSZ[j]-7) + square(sigCR[j])/2 + lgSz_adj[j]) );
    // Mean HT/itm, by prey type, adjusted for mean prey size and lognormal dist
    HT[j] = fmin(900, exp(psi1[j] + psi2[j] * (2.5*muSZ[j]-7) + square(sigHT[j])/2)) ;
    LM[j] = inv_logit(lgtLM[j]) ;
  }
  // Mean HT/item for Unid prey:
  HT_u = exp(psi1_u + psi2_u * (2.5*muSZ_u-7) + square(sigHT_u)/2);
  // Mean Energy intake rate (ER) by prey, incl. uncertainty in Caloric density
  ER = Cal_dens .* CR ;
  // Proportional contribution (biomass consumed) of each prey type to diet: 
  Pi = (eta .* CR) / sum(eta .* CR) ;
  // Overall mean consumption rate (CRmn) given effort allocation to each prey: 
  CRmn = sum(eta .* CR) ;
  // Overall mean Energy Intake Rate (ERmn, kcal.min) given effort allocation: 
  ERmn = sum(eta .* CR .* Cal_dens) ;
  // Overall mean dive success rate (Lambda)
  LMmn = sum(eta .* LM) ;
}
