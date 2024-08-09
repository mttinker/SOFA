// SOFAfit
// This Stan code fits a model of foraging behavior to observational data 
//  on sea otter foraging  [Refer to "SOFA_summary" for details] 
//  - this is simple version for "Non-grouped" data
//  - *** ECRate in this version = encounter/attack/capture rate (alpha in a type-II fctn response)
// Section 1. Data inputs for model
data {
  int<lower=1> Nbouts;                  // Number of feeding bouts (for effort allocation observations)
  int<lower=1> K;                       // Number of prey types (1 to (K-1) for ID prey, UN-ID prey = K)
  int<lower=1> Km1;                     // "K - 1", number of prey types excluding UN-ID prey  
  array[Nbouts,K] int<lower=0> EffortP; // Minutes allocated to each prey type (incl. un-id) by bout  
  int<lower=1> NSz;                     // Number of prey Size obs
  int<lower=1> NHt;                     // Number of prey Handling time (HT) obs  
  int<lower=1> NCR;                     // Number of prey Consumption rate obs
  int<lower=1> NLm;                     // Number of dive success rate obs (lambda)
  array[2] int<lower=1> NU;             // Number of un-ID prey obs (2 different param types, FOR NOW)
  array[NU[1]] real<lower=0> SZmnU;     // Array of mean prey size vaues (mm) for UN-ID prey, by bout  
  array[NU[2]] real<lower=0> HTmnU;     // Array of mean handling times (sec), for UN-ID prey, by bout  
  array[NSz] int<lower=1,upper=Km1> Sp; // Array of prey id values for size obs
  array[NHt] int<lower=1,upper=Km1> Hp; // Array of prey id values for HT obs  
  array[NCR] int<lower=1,upper=Km1> Cp; // Array of prey id values for Encounter/Capture Rate obs
  array[NSz] real<lower=0> SZmn;        // Array of mean prey size vaues (mm) by prey, by bout  
  array[NHt] real<lower=0> HTmn;        // Array of mean handling times (sec), by prey, by bout   
  array[NCR] real<lower=0> ECRate;      // Array of prey-specific encounter/capture rates by bout (#/min)
  vector[NCR] Csz;                      // Array of prey size (log pwr transform), covariate for encounter rates
  vector[NCR] Css;                      // Array of sample sizes for encounter rates (N obs per bout)
  vector[NSz] Sss;                      // Array of sample sizes (n obs per bout) for size observations 
  vector[NU[1]] Sss_u;                  // Array of sample sizes (n obs per bout) for size observations, Un-ID   
  vector[NHt] Hsz;                      // Array of prey size (log pwr transform), covariate for HT obs
  vector[NHt] Hss;                      // Array of sample sizes for HT obs (n obs per bout)
  vector[NU[2]] Hss_u;                  // Array of sample sizes for HT obs (n obs per bout), UN-ID
  vector[NU[2]] Hsz_u;                  // Array of prey size (log pwr transform), covariate for UN-ID HT obs
  vector[Km1] Cal_dns_mn;               // Vector of log caloric density values (kcal/g) by prey type, mean
  vector<lower=0>[Km1] Cal_dns_sg;      // Vector of sdt err of log caloric density values (kcal/g) by prey type
  vector<lower=0>[Km1] logMass_sg;      // Vector of std err vlues for log-log fxn of bomass vs size
  array[NLm] int<lower=1,upper=Km1> Lp; // Array of prey id values for success rate obs
  array[NLm] real LMlg;                 // Array of logit(lambda), proportion of succesful dives by prey, by bout  
  vector[NLm] Lss;                      // Array of sample sizes (n obs per bout) for logit lambda observations 
  real CR_max ;                         // average maximum possible mean Consumption rate 
  vector[Km1] eta_pri ;                 // dirichlet prior for eta
  vector[Km1] MnMass ;                  // Mean biomass by prey type
  vector[Km1] MLpar1 ;                  // Param 1 of mass vs biomass log-log fxn, by prey type
  vector[Km1] MLpar2 ;                  // Param 2 of mass vs biomass log-log fxn, by prey type
  array[K] corr_matrix[2] rho_HS ;      // Prey-specific Correlation matrices between log(HT) and log(SZ)
}
transformed data{
	// 
}
// Section 2. The parameters to be estimated 
parameters {
  simplex[Km1] eta ;                    // mean proportional effort allocation by prey type (if all ID'd)
  real<lower=0,upper=5> tauB ;          // relative precision of diet composition across bouts
  real<lower=0,upper=1> kappa;          // max prob un-id across prey types, average
  vector<lower=0,upper=6>[Km1] muSZ ;   // median of log(Size_mm) of each prey
  real<lower=0,upper=6> muSZ_u ;        // median of log(Size_mm) of UN-ID
  vector[Km1] lgtLM ;                   // median logit(lambda) (dive success rate) for each prey
  vector<lower=0,upper=2.5>[Km1] sigSZ ;// stdev of log(Size_cm) of each prey 
  vector<lower=0,upper=2.5>[Km1] sigHT ;// stdev of log(HT per item) of each prey 
  vector<lower=0,upper=5>[Km1] sigCR ;  // stdev of log(ECR per item) of each prey 
  vector<lower=0,upper=10>[Km1] sigLM ; // stdev of logit(lambda) of each prey 
  real<lower=0,upper=2.5> sigSZ_u ;     // stdev of log(Size_cm) of UN-ID
  real<lower=0,upper=2.5> sigHT_u ;     // stdev of log(HT per item) of UN-ID    
  vector<lower=-10,upper=10>[Km1] phi1 ;// intercept of ECRate v Size fxn by prey 
  vector<lower=0>[Km1] phi2 ;           // slope of ECRate v Size fxn by prey (assumed negative)
  vector<lower=0>[Km1] psi1 ;           // intercept of HT v Size fxn by prey 
  vector<lower=0>[Km1] psi2 ;           // slope of HT v Size fxn by prey   
  real<lower=0> psi1_u ;                // intercept of HT v Size fxn for UN-ID prey 
  real<lower=0> psi2_u ;                // slope of HT v Size fxn for UN-ID prey     
}
// Section 3. Derived (transformed) parameters
transformed parameters {
  vector<lower=0,upper=1>[Km1] Omega ;  // scaled similarity to UNID: relative prob of contributing to un-ID
  vector[K] theta ;                     // prey-specific effort by prey type, average, including Un-ID
  vector[Km1] muHT ;                    // median of log(HT per item) of each prey  
  real muHT_u ;                         // median of log(HT per item) of UN-ID       
  vector[2] mu1 ;                       // vector of median log size and log HT for Un-ID
  matrix[2,2] Sig1 ;                    // covar matrix of log size and log HT for Un-ID   
  //
  // Compute mean expected HT at avg size for un-id prey
  muHT_u = psi1_u + psi2_u * (2.5 * muSZ_u - 7);  
  // Loop through prey types, calculate prob of ID-ing prey (Pid) & contribution to un-ID 
  // NOTE: Bhattacharyya coefficient measures overlap between size/HT distributions of 
  // each prey type and UN-ID prey: if NO overlap then prey does not contribute to un-ID,
  // while if maximum overlap, then prob of Un-ID is at maximum (maxPunid)
  mu1 = [muSZ_u, muHT_u]' ;
  // "quad_form_diag" generates covariance matrix from correlation matrix and vector of std errors
  Sig1 = quad_form_diag(rho_HS[K], [sigSZ_u, sigHT_u]') ;
  for (j in 1:Km1){    
    vector[2] mu2 ;                     // vector of median log size and log HT for prey type j
    matrix[2,2] Sig2 ;                  // covar matrix of log size and log HT for prey type j
    real BDist ;                        // Bhattacharyya distance (bivariate, SZ/HT) b/t UnID & prey type j
    muHT[j] = psi1[j] + psi2[j] * (2.5 * muSZ[j] - 7) ;  
    mu2 = [muSZ[j], muHT[j]]' ;
    Sig2 = quad_form_diag(rho_HS[j], [sigSZ[j], sigHT[j]]') ;
    // NOTE: for efficiency, we replace inverse of covariance matrix with matrix division 
    BDist = 0.125 * ( (mu1 - mu2)' / ((Sig1 + Sig2) ./ 2) * (mu1-mu2) ) + 
          0.5 * log(determinant((Sig1 + Sig2) ./ 2) / sqrt(determinant(Sig1) * determinant(Sig2))) ;
    Omega[j] = exp(-BDist) ; // bivariate Bhattacharyya coefficient
  }  
  // scaled similarity to UNID: relative prob of contributing to un-ID
  Omega = ( Omega ./ max(Omega) ) * kappa ; 
  // mean expected effort allocation by prey, with UN-ID
  theta = append_row(eta .* (1 - Omega), sum(eta .* Omega ) ) ;
}
// Section 4. Estimating model parameters (drawing from probability distributions)
model {
  // A) Observed nodes:
  // Allocation of effort by prey type (proportional # minutes of each feeding bout)     append_row(vector x, real y)
  for (i in 1:Nbouts){
    EffortP[i] ~ dirichlet_multinomial( theta * tauB ) ;
  }
  // Bout attributes:
  // Observed encounter rate during dives (#/min) by prey type, random samples
  ECRate ~ lognormal((phi1[Cp] - phi2[Cp] .* Csz), sigCR[Cp]) ;
  // Observed Handling time (sec), random samples, ID prey & UN-id prey
  HTmn ~ lognormal((psi1[Hp] + psi2[Hp] .* Hsz), sigHT[Hp]) ;
  HTmnU ~ lognormal((psi1_u + psi2_u * Hsz_u), sigHT_u) ;  
  // Observed Mean prey size (mm), random samples, ID prey & UN-id prey
  SZmn ~ lognormal(muSZ[Sp],sigSZ[Sp]) ;
  SZmnU ~ lognormal(muSZ_u, sigSZ_u) ;
  // Observed logit(lambda), dive succes rate by bout/prey, random samples 
  LMlg ~ normal(lgtLM[Lp], sigLM[Lp]) ;
  //
  // B) Prior distributions for model parameters:
  eta ~ dirichlet(eta_pri) ; 
  kappa ~ beta(2,3);                    // prob un-ID by bout, scaled by similarity to UN-ID in size & HT
  tauB ~ cauchy(0,.5) ;                 // precision param for diet comp variation across bouts
  phi1 ~ normal(0,1.5) ;                // Intercept of log Encounter Rate, by prey
  phi2 ~ cauchy(0,.1) ;                 // effect of size on log Cons rate
  psi1 ~ lognormal(0.5,0.6) ;           // Intercept of log HT, by prey
  psi2 ~ cauchy(0,.1) ;                 // effect of size on log HT rate
  psi1_u ~ lognormal(0.5,0.6) ;         // Intercept of log HT, UN-ID prey
  psi2_u ~ cauchy(0,.1) ;               // effect of size on log HTm UN-ID prey
  muSZ ~ normal(4,1);                   // log-mean size by prey type
  muSZ_u ~ normal(4,1);                 // log-mean size of UN-ID prey 
  sigSZ ~ cauchy(0,.1);                 // variation in prey size across bouts, by prey type
  sigHT ~ cauchy(0,.1);                 // variation in prey HT across bouts, by prey type
  sigCR ~ cauchy(0,.1);                 // variation in prey ECR across bouts, by prey type
  sigSZ_u ~ cauchy(0,.1);               // variation in prey size across bouts, unid prey
  sigHT_u ~ cauchy(0,.1);               // variation in prey HT across bouts, unid prey
  lgtLM ~ normal(1,2);                  // mean logit(lambda), dive success rate, by prey
  sigLM ~ cauchy(0,1.5);                // variation in logit(lambda), dive success rate
}
// Section 5. Derived parameters and statistics 
generated quantities {
  vector[Km1] upsilon ;                 // relative contribution to unknown prey dives
  vector[Km1] SZ = exp( muSZ + square(sigSZ)/2 ); // Mean Size (mm) by prey type (adjust for log-normal)
  real SZ_u = exp(muSZ_u + square(sigSZ_u)/2 ) ; // Mean size (mm) of Un_ID prey
  real HT_u = exp(psi1_u + psi2_u * (2.5* log(SZ_u) - 7) + square(sigHT_u)/2) ; // Mean HT/item for Un_ID prey
  vector[Km1] HT = rep_vector(0,Km1) ;  // mean handling time by prey
  vector[Km1] CP = rep_vector(0,Km1) ;  // mean encounter/capture rate by prey
  vector[Km1] CR = rep_vector(0,Km1) ;  // mean consumption rate by prey
  vector[Km1] ER = rep_vector(0,Km1) ;  // mean energy rate by prey
  vector[Km1] Pi ;                      // mean dietary composition by prey  
  real CRmn = 0 ;                       // overall mean consumption rate
  real ERmn = 0 ;                       // overall mean energy intake rate
  vector[Km1] LM = inv_logit(lgtLM) ;   // mean dive success rate y prey
  real LMmn = sum(eta .* LM) ;          // overall mean dive success rate
  //
  // Loop through 1000 random bouts, to characterize uncertainty
  for (r in 1:1000){
    vector[Km1] s ;
    vector[Km1] m ;
    vector[Km1] a ;
    vector[Km1] h ;
    vector[Km1] Cal_dens ;
    real CR_mx = normal_rng(CR_max, .15 * CR_max) ;
    for (j in 1:Km1){
      real Csz_r ;        
      Cal_dens[j] = fmax(.15,fmin(2.25,lognormal_rng(Cal_dns_mn[j],Cal_dns_sg[j]))) ;
      s[j] = lognormal_rng(muSZ[j],sigSZ[j]) ;
      m[j] = exp(MLpar1[j] + MLpar2[j] * log(s[j]) + normal_rng(0,logMass_sg[j])) ;
      Csz_r = 2.5 * log(s[j]) - 7 ;      
      a[j] = lognormal_rng((phi1[j] - phi2[j] * Csz_r), sigCR[j]) ;
      h[j] = lognormal_rng((psi1[j] + psi2[j] * Csz_r), sigHT[j]) / 60 ;
      // Prey-specific consumption & energy intake rates: type-II functional response (Hollings disc eqn)
      CR[j] = CR[j] + fmin(CR_mx, (a[j] * m[j]) / (1 + h[j] * a[j]) ) ;  
      ER[j] = ER[j] + fmin(CR_mx, (a[j] * m[j]) / (1 + h[j] * a[j]) ) * Cal_dens[j] ;
      HT[j] = HT[j] + fmin(600, h[j]*60) ;
      CP[j] = CP[j] + a[j] ;
    }
  }
  CR = CR ./ 1000 ;
  ER = ER ./ 1000 ;
  HT = HT ./ 1000 ;
  CP = CP ./ 1000 ;
  CRmn = sum(eta .* CR) ;
  ERmn = sum(eta .* ER) ;
  // Proportional contribution (biomass consumed) of each prey type to diet: 
  Pi = (eta .* CR) / sum(eta .* CR) ;
  // Proportional contribution of items of each prey type to unknown prey items 
  upsilon = (100 .* Pi ./ MnMass) .* Omega  ;
  upsilon = upsilon ./ sum(upsilon) ;
}
