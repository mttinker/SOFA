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
  array[Nbouts,K] int<lower=0> EffortP; // Minutes allocated to each prey type (incl. un-id) by bout  
  int<lower=1> NSz;               // Number of prey Size obs
  int<lower=1> NHt;               // Number of prey Handling time (HT) obs  
  int<lower=1> NCR;               // Number of prey Consumption rate obs
  int<lower=1> NLm;               // Number of dive success rate obs (lambda)
  array[2] int<lower=1> NU;             // Number of un-ID prey obs (2 different param types, FOR NOW)
  array[NU[1]] real<lower=0> SZmnU;     // Array of mean prey size vaues (mm) for UN-ID prey, by bout  
  array[NU[2]] real<lower=0> HTmnU;     // Array of mean handling times (sec), for UN-ID prey, by bout  
  array[NSz] int<lower=1,upper=Km1> Sp; // Array of prey id values for size obs
  array[NHt] int<lower=1,upper=Km1> Hp; // Array of prey id values for HT obs  
  array[NCR] int<lower=1,upper=Km1> Cp; // Array of prey id values for Cons Rate obs
  array[NSz] real<lower=0> SZmn;        // Array of mean prey size vaues (mm) by prey, by bout  
  array[NHt] real<lower=0> HTmn;        // Array of mean handling times (sec), by prey, by bout   
  array[NCR] real<lower=0> CRate;       // Array of prey-specific consumption rates (CR) by bout (g/min)
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
  array[NLm] int<lower=1,upper=Km1> Lp; // Array of prey id values for success rate obs
  array[NLm] real LMlg;           // Array of logit(lambda), proportion of succesful dives by prey, by bout  
  vector[NLm] Lss;                // Array of sample sizes (n obs per bout) for logit lambda observations 
  vector<lower=0>[Km1] CR_max ;   // maximum possible mean CR per prey type
  vector[Km1] eta_pri ;           // dirichlet prior for eta
  vector[Km1] MnMass ;            // Mean biomass by prey type
}
transformed data{
	vector[Km1] logMass_mn = rep_vector(0, Km1) ; // -1 * square(logMass_sg)/2 ;
}
// Section 2. The parameters to be estimated 
parameters {
  simplex[Km1] eta ;                    // mean proportional effort allocation by prey type (if all ID'd)
  real<lower=0,upper=20> tauB ;         // relative precision of diet composition across bouts
  array[Nbouts] simplex[K] theta ;      // vectors of bout-specific prey-allocation probs for multinomial  
  vector<lower=0,upper=6>[Km1] muSZ ;   // mean of log(Size_mm) of each prey
  real<lower=0,upper=6> muSZ_u ;        // mean of log(Size_mm) of UN-ID
  vector[Km1] lgtLM ;                   // mean logit(lambda) (dive success rate) for each prey
  vector<lower=0,upper=20>[Km1] sigSZ ; // stdev of log(Size_cm) of each prey 
  vector<lower=0,upper=20>[Km1] sigHT ; // stdev of log(HT per item) of each prey 
  vector<lower=0,upper=20>[Km1] sigCR ; // stdev of log(HT per item) of each prey 
  vector<lower=0,upper=20>[Km1] sigLM ; // stdev of logit(lambda) of each prey 
  real<lower=0,upper=20> sigSZ_u ;      // stdev of log(Size_cm) of UN-ID
  real<lower=0,upper=20> sigHT_u ;      // stdev of log(HT per item) of UN-ID    
  vector<lower=0>[Km1] phi1 ;           // intercept of CRate v Size fxn by prey 
  vector<lower=0>[Km1] phi2 ;           // slope of CRate v Size fxn by prey 
  vector<lower=0>[Km1] psi1 ;           // intercept of HT v Size fxn by prey 
  vector<lower=0>[Km1] psi2 ;           // slope of HT v Size fxn by prey   
  real<lower=0> psi1_u ;                // intercept of HT v Size fxn for UN-ID prey 
  real<lower=0> psi2_u ;                //  slope of HT v Size fxn for UN-ID prey     
  real<lower=0,upper=1> maxPunid;       // max prob un-id across prey types
  vector<lower=0>[Km1] Cal_dens;        // Caloric density (kcal/g) by prey type
  vector[Km1] lgSz_adj;                 // Uncertainty adjustment for size-biomass fxn
}
// Section 3. Derived (transformed) parameters
transformed parameters {
	vector<lower=0,upper=1>[Km1] UNID_sim ;// prey-specific similarity (overlap) with UNID
  vector<lower=0,upper=1>[Km1] Omega ;  // prey-specific probability of positive ID
  vector<lower=0>[K] alpha ;            // vector of dirichlet params: relative prey freq (inc. UN-ID)
  vector[Km1] muHT ;                    // mean of log(HT per item) of each prey  
  real muHT_u ;                         // mean of log(HT per item) of UN-ID       
  vector[2] mu1 ;                       // vector of mean log size and log HT for Un-ID
  matrix[2,2] Sig1 ;                    // covar matix of log size and log HT for Un-ID
  //
  // Compute mean expected HT at avg size for un-id prey
  muHT_u = psi1_u + psi2_u * (2.5 * muSZ_u - 7);  
  // Loop through prey types, calculate prob of ID-ing prey (Pid) & contribution to un-ID 
  // NOTE: Bhattacharyya coefficient measures overlap between size/HT distributions of 
  // each prey type and UN-ID prey: if NO overlap then prey does not contribute to un-ID,
  // while if maximum overlap, then prob of Un-ID is at maximum (maxPunid)
  mu1 = [muSZ_u, muHT_u]' ;
  Sig1 = [ [ square(sigSZ_u), .5 * sigSZ_u * sigHT_u ], [.5 * sigSZ_u * sigHT_u, square(sigHT_u) ] ] ;
  for (j in 1:Km1){    
    vector[2] mu2 ;                     // vector of mean log size and log HT for prey type j
    matrix[2,2] Sig2 ;                  // covar matix of log size and log HT for prey type j
    real BDist ;                        // Bhattacharyya distance (bivariate, SZ/HT) b/t UnID & prey type j
    muHT[j] = psi1[j] + psi2[j] * (2.5 * muSZ[j] - 7) ;  
    mu2 = [muSZ[j], muHT[j]]' ;
    Sig2 = [ [ square(sigSZ[j]), .5 * sigSZ[j] * sigHT[j] ] , [.5 * sigSZ[j] * sigHT[j], square(sigHT[j]) ] ] ;
    // NOTE: for efficiency, replace inverse of covariance matrix with matrix division 
    BDist = 0.125 * ( (mu1 - mu2)' / ((Sig1 + Sig2) ./ 2) * (mu1-mu2) ) + 
          0.5 * log(determinant((Sig1 + Sig2) ./ 2) / sqrt(determinant(Sig1) * determinant(Sig2))) ;
    UNID_sim[j] = exp(-BDist) ; // bivariate Bhattacharyya coefficient
  }  
  Omega = 1 - maxPunid * ( UNID_sim ./ max(UNID_sim) ) ; 
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
  // Bout attributes:
  // NOTE: For lognormal comparisons of observed to expected attributes, the number 
  //    of dives going into a mean per bout value must be accounted for... specifically,
  //    the log median of the means of lognormal samples of size "n" is calculated as: 
  //           mu + (n*sg^2 - sg^2)/(2*n)    
  // Observed consumption rate (g/min) by prey type, random samples
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
  lgSz_adj ~ normal(logMass_mn,logMass_sg) ; // CRate adjust (incorporates mass-lng est. uncertainty)
  // NOTE: also include normal prior for CR_uncert multiplier, with mean of 1,
  //   to incorporate uncertainty associated with Mass-size power functions
  eta ~ dirichlet(eta_pri) ; 
  tauB ~ cauchy(0,1) ;                  // precision param for diet comp variation across bouts
  phi1 ~ cauchy(0,.5) ;                 // Intercept of log Cons Rate, by prey
  phi2 ~ cauchy(0,.1) ;                 // effect of size on log Cons rate
  psi1 ~ cauchy(0,1) ;                  // Intercept of log HT, by prey
  psi2 ~ cauchy(0,.1) ;                 // effect of size on log HT rate
  psi1_u ~ cauchy(0,1) ;                // Intercept of log HT, UN-ID prey
  psi2_u ~ cauchy(0,.1) ;               // effect of size on log HTm UN-ID prey
  muSZ ~ normal(4,1);                   // log-mean size by prey type
  muSZ_u ~ normal(4,1);                 // log-mean size of UN-ID prey
  maxPunid ~ beta(1,1);                 // max prob un-ID, for prey overlapping UN-ID in size & HT
  sigSZ ~ cauchy(0,.5);                 // variation in prey size across bouts, by prey type
  sigHT ~ cauchy(0,.5);                 // variation in prey HT across bouts, by prey type
  sigCR ~ cauchy(0,.5);                 // variation in prey CR across bouts, by prey type
  sigSZ_u ~ cauchy(0,.5);               // variation in prey size across bouts, unid prey
  sigHT_u ~ cauchy(0,.5);               // variation in prey HT across bouts, unid prey
  lgtLM ~ normal(1,2);                  // mean logit(lambda), dive success rate, by prey
  sigLM ~ cauchy(0,1.5);                // variation in logit(lambda), dive success rate
}
// Section 5. Derived parameters and statistics 
generated quantities {
  vector[Km1] upsilon ;                 // relative contribution to unknown prey dives
  vector[Km1] SZ ;
  vector[Km1] logSZ ;
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
  logSZ = muSZ + square(sigSZ)/2 ;
  SZ = exp(muSZ + square(sigSZ)/2 ) ;
  SZ_u = exp(muSZ_u + square(sigSZ_u)/2 ) ;
  for (j in 1:Km1){
    // Mean Cons rate (CR, g/min) by prey type, adjusted for mean prey size and lognormal dist
    CR[j] = fmin(CR_max[j], exp(phi1[j] + phi2[j] * (2.5*logSZ[j]-7) + square(fmin(1.5,sigCR[j]))/2 + lgSz_adj[j]) ) ;
    // Mean HT/itm, by prey type, adjusted for mean prey size and lognormal dist
    HT[j] = fmin(600, exp(psi1[j] + psi2[j] * (2.5*logSZ[j]-7) + square(sigHT[j])/2)) ;
    LM[j] = inv_logit(lgtLM[j]) ;
  }
  // Mean HT/item for Unid prey:
  HT_u = exp(psi1_u + psi2_u * (2.5*muSZ_u-7) + square(sigHT_u)/2);
  // Mean Energy intake rate (ER) by prey, incl. uncertainty in Caloric density
  ER = Cal_dens .* CR ;
  // Proportional contribution (biomass consumed) of each prey type to diet: 
  Pi = (eta .* CR) / sum(eta .* CR) ;
  // Proportional contribution of each prey type to unkown prey items 
  upsilon = (1000 .* Pi ./ MnMass) .* (1 - Omega) ;
  upsilon = upsilon ./ sum(upsilon) ;
  // Overall mean consumption rate (CRmn) given effort allocation to each prey: 
  CRmn = sum(eta .* CR) ;
  // Overall mean Energy Intake Rate (ERmn, kcal.min) given effort allocation: 
  ERmn = sum(eta .* CR .* Cal_dens) ;
  // Overall mean dive success rate (Lambda)
  LMmn = sum(eta .* LM) ;
}
