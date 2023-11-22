// SOFAfit
// This Stan code fits a model of foraging behavior to observational data 
//  on sea otter foraging  [Refer to "SOFA_summary" for details] 
//  - this version analzes data categorized by some grouping variable
//
// Section 1. Data inputs for model
data {
  int<lower=1> Nbouts;            // Number of feeding bouts (for effort allocation observations)
  int<lower=1> Ngrp ;             // Number of feeding bouts (for effort allocation observations)
  int<lower=1> K;                 // Number of prey types (1 to (K-1) for ID prey, UN-ID prey = K)
  int<lower=1> Km1;               // "K - 1", number of prey types excluding UN-ID prey  
  array[Nbouts,K] int<lower=0> EffortP; // Minutes allocated to each prey type (incl. un-id) by bout  
  array[Nbouts] int<lower=1> GrpE;      // Group ID for each bout (for effort allocations)  
  int<lower=1> NSz;               // Number of prey Size obs
  int<lower=1> NHt;               // Number of prey Handling time (HT) obs  
  int<lower=1> NCR;               // Number of prey Consumption rate obs
  int<lower=1> NLm;               // Number of dive success rate obs (lambda)
  array[2] int<lower=1> NU;             // Number of un-ID prey obs (2 different param types, FOR NOW)
  array[NU[1]] real<lower=0> SZmnU;     // Array of mean prey size vaues (mm) for UN-ID prey, by bout  
  array[NU[2]] real<lower=0> HTmnU;     // Array of mean handling times (sec), for UN-ID prey, by bout 
  array[NU[1]] int<lower=1,upper=Ngrp> Sg_u; // Array of group id values for size obs of UN-ID
  array[NU[2]] int<lower=1,upper=Ngrp> Hg_u; // Array of group id values for HT obs of UN-ID  
  array[NSz] int<lower=1,upper=Km1> Sp; // Array of prey id values for size obs
  array[NHt] int<lower=1,upper=Km1> Hp; // Array of prey id values for HT obs  
  array[NCR] int<lower=1,upper=Km1> Cp; // Array of prey id values for Cons Rate obs
  array[NSz] int<lower=1,upper=Ngrp> Sg; // Array of group id values for size obs
  array[NHt] int<lower=1,upper=Ngrp> Hg; // Array of group id values for HT obs  
  array[NCR] int<lower=1,upper=Ngrp> Cg; // Array of group id values for Cons Rate obs
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
  array[NLm] int<lower=1,upper=Ngrp> Lg;// Array of group id values for Dive success rate
  array[NLm] real LMlg;                 // Array of logit(lambda), proportion of succesful dives by prey, by bout  
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
  array[Ngrp] simplex[Km1] etaG ;       // mean proportional effort allocation by prey type, by group
  array[Ngrp] real<lower=0,upper=20> tauB ; // relative precision of diet composition across bouts
  real<lower=0,upper=100> tauG ;        // relative precision of diet composition across groups
  array[Nbouts] simplex[K] theta ;      // vectors of bout-specific prey-allocation probs for multinomial  
  vector<lower=0,upper=6>[Km1] muSZ ;   // mean of log(Size_mm) of each prey type
  array[Ngrp] vector<lower=0,upper=6>[Km1] muSZG ; // mean of log(Size_mm) of each prey type, by group 
  real<lower=0,upper=6> muSZ_u ;        // mean of log(Size_mm) of UN-ID
  array[Ngrp] real<lower=0,upper=6> muSZG_u ; // mean of log(Size_mm) of UN-ID, by group
  vector[Km1] lgtLM ;                   // mean logit(lambda) (dive success rate) for each prey
  array[Ngrp] vector[Km1] lgtLMG ;      // mean logit(lambda) (dive success rate) for each prey, by group
  vector<lower=0,upper=20>[Km1] sigLM ; // stdev of logit(lambda) of each prey 
  vector<lower=0,upper=20>[Km1] sigSZ ; // stdev of log(Size_cm) of each prey 
  vector<lower=0,upper=20>[Km1] sigHT ; // stdev of log(HT per item) of each prey 
  vector<lower=0,upper=20>[Km1] sigCR ; // stdev of log(HT per item) of each prey 
  real<lower=0,upper=20> sigSZ_u ;      // stdev of log(Size_mm) of UN-ID
  real<lower=0,upper=20> sigHT_u ;      // stdev of log(HT per item) of UN-ID    
  vector<lower=0>[Km1] phi1 ;           // intercept of CRate v Size fxn by prey 
  array[Ngrp] vector<lower=0>[Km1] phi1G ; // intercept of CRate v Size fxn by prey 
  vector<lower=0>[Km1] phi2 ;           // slope of CRate v Size fxn by prey 
  vector<lower=0>[Km1] psi1 ;           // intercept of HT v Size fxn by prey 
  array[Ngrp] vector<lower=0>[Km1] psi1G ;// intercept of HT v Size fxn by prey 
  vector<lower=0>[Km1] psi2 ;           // slope of HT v Size fxn by prey   
  real<lower=0> psi1_u ;                // intercept of HT v Size fxn for UN-ID prey 
  array[Ngrp] real<lower=0> psi1G_u ;   // intercept of HT v Size fxn for UN-ID prey 
  real<lower=0> psi2_u ;                //  slope of HT v Size fxn for UN-ID prey     
  vector<lower=0,upper=1>[Ngrp] maxPunid;// max prob un-id across prey types, by group
  vector<lower=0>[Km1] Cal_dens;        // Caloric density (kcal/g) by prey type
  vector[Km1] lgSz_adj; // Uncertainty adjustment for size-biomass fxn
  real<lower=0,upper=2.5> sg1 ;         // stdev of hierarchical param 
  real<lower=0,upper=2.5> sg2 ;         // stdev of hierarchical param 
  real<lower=0,upper=2.5> sg3 ;         // stdev of hierarchical param 
  real<lower=0,upper=2.5> sg4 ;         // stdev of hierarchical param 
  real<lower=0,upper=2.5> sg5 ;         // stdev of hierarchical param 
  real<lower=0,upper=2.5> sg6 ;         // stdev of hierarchical param 
}
// Section 3. Derived (transformed) parameters
transformed parameters {
  array[Ngrp] vector<lower=0,upper=1>[Km1] UNID_sim;//  prey-specific similarity (overlap) with UNID
  array[Ngrp] vector<lower=0,upper=1>[Km1] OmegaG ; // prey-specific probability of positive ID
  array[Ngrp] vector<lower=0>[K] alpha ;            // vector of dirichlet params: relative prey freq (inc. UN-ID)
  //
  for (g in 1:Ngrp){  
    vector[Km1] muHT ;                  // mean of log(HT per item) of each prey  
    real muHT_u ;                       // mean of log(HT per item) of UN-ID   
    vector[2] mu1 ;                     // vector of mean log size and log HT for Un-ID
    matrix[2,2] Sig1 ;                  // covar matix of log size and log HT for Un-ID
    // Compute mean expected HT at avg size for un-id prey
    muHT_u = psi1G_u[g] + psi2_u * (2.5 * muSZG_u[g] - 7);  
    // Loop through prey types, calculate prob of ID-ing prey (Pid) & contribution to un-ID 
    // NOTE: Bhattacharyya coefficient measures overlap between size/HT distributions of 
    // each prey type and UN-ID prey: if NO overlap then prey does not contribute to un-ID,
    // while if maximum overlap, then prob of Un-ID is at maximum (maxPunid)
    mu1 = [muSZG_u[g], muHT_u]' ;
    Sig1 = [ [ square(sigSZ_u), .5 * sigSZ_u * sigHT_u ], [.5 * sigSZ_u * sigHT_u, square(sigHT_u) ] ] ;
    for (j in 1:Km1){
      vector[2] mu2 ;     // vector of mean log size and log HT for prey type j
      matrix[2,2] Sig2 ;  // covar matix of log size and log HT for prey type j
      real BDist ;        // Bhattacharyya distance (bivariate, SZ/HT) b/t UnID & prey type j
      muHT[j] = psi1G[g][j] + psi2[j] * (2.5 * muSZG[g][j] - 7) ; 
      mu2 = [muSZG[g][j], muHT[j]]' ;
      Sig2 = [ [ square(sigSZ[j]), .5 * sigSZ[j] * sigHT[j] ] , [.5 * sigSZ[j] * sigHT[j], square(sigHT[j]) ] ] ;
      // NOTE: for efficiency, replace inverse of covariance matrix with matrix division 
      BDist = 0.125 * ( (mu1 - mu2)' / ((Sig1 + Sig2) ./ 2) * (mu1-mu2) ) + 
            0.5 * log(determinant((Sig1 + Sig2) ./ 2) / sqrt(determinant(Sig1) * determinant(Sig2))) ;
      UNID_sim[g][j] = exp(-BDist) ; // bivariate Bhattacharyya coefficient
    } 
    OmegaG[g] = 1 - maxPunid[g] * ( UNID_sim[g] ./ max(UNID_sim[g]) ) ;
    // Calculate alpha, the dirichlet params for bout-specific prey allocation probs
    alpha[g][1:Km1] = etaG[g] .* OmegaG[g] ;
    alpha[g][K] = sum(etaG[g] .* (1 - OmegaG[g])) ;
  }
}
// Section 4. Estimating model parameters (drawing from probability distributions)
model {
  // A) Observed nodes:
  // Allocation of effort by prey type (proportional # minutes of each feeding bout)
  for (i in 1:Nbouts){
    theta[i] ~ dirichlet(alpha[GrpE[i]] * tauB[GrpE[i]]);  
    EffortP[i,] ~ multinomial(theta[i]) ;
  }
  // Bout attributes:
  // NOTE: For lognormal comparisons of observed to expected attributes, the number 
  //    of dives going into a mean per bout value must be accounted for... specifically,
  //    the log median of the means of lognormal samples of size "n" is calculated as: 
  //           mu + (n*sg^2 - sg^2)/(2*n)    
  // Observed consumption rate (g/min) by prey type, random samples
  for(i in 1:NCR){
    CRate[i] ~ lognormal((phi1G[Cg[i]][Cp[i]] + phi2[Cp[i]] * Csz[i]) + ((Css[i] * square(sigCR[Cp[i]])) - square(sigCR[Cp[i]])) / (2 * Css[i]),
            sigCR[Cp[i]] / sqrt(Css[i]) ) ;
  }
  // Observed Handling time (sec), random samples, ID prey & UN-id prey
  for(i in 1:NHt){
    HTmn[i] ~ lognormal( (psi1G[Hg[i]][Hp[i]] + psi2[Hp[i]] * Hsz[i]) + ((Hss[i] * square(sigHT[Hp[i]])) - square(sigHT[Hp[i]])) / (2 * Hss[i]),
              sigHT[Hp[i]] / sqrt(Hss[i]) ) ;
  }
  for(i in 1:NU[2]){
    HTmnU ~ lognormal((psi1G_u[Hg_u[i]] + psi2_u * Hsz_u[i]) + ((Hss_u[i] * square(sigHT_u)) - square(sigHT_u)) / (2 * Hss_u[i]),
            sigHT_u / sqrt(Hss_u[i]) ) ;  
  }
  // Observed Mean prey size (mm), random samples, ID prey & UN-id prey
  for(i in 1:NSz){
    SZmn[i] ~ lognormal(muSZG[Sg[i]][Sp[i]] + ((Sss[i] * square(sigSZ[Sp[i]])) - square(sigSZ[Sp[i]])) / (2 * Sss[i]), 
              sigSZ[Sp[i]] / sqrt(Sss[i]) ) ;
  }  
  for(i in 1:NU[1]){
    SZmnU[i] ~ lognormal(muSZG_u[Sg_u[i]] + ((Sss_u[i] * square(sigSZ_u)) - square(sigSZ_u)) / (2 * Sss_u[i]), 
            sigSZ_u / sqrt(Sss_u[i]) ) ;
  }
  // Observed logit(lambda), dive succes rate by bout/prey, random samples 
  for(i in 1:NLm){
    LMlg[i] ~ normal(lgtLMG[Lg[i]][Lp[i]], sigLM[Lp[i]] / sqrt(Lss[i]) ) ; 
  }
  //  
  // B) Prior distributions for model parameters:
  Cal_dens ~ normal(Cal_dns_mn,Cal_dns_sg) ; // Caloric density (incorporates est. uncertainty)
  lgSz_adj ~ normal(logMass_mn,logMass_sg) ;          // CRate adjust (incorporates est. uncertainty)
  eta ~ dirichlet(eta_pri) ; 
  for (g in 1:Ngrp){
    etaG[g] ~ dirichlet(eta * 20 * tauG);
    phi1G[g] ~ normal(phi1,sg1);
    psi1G[g] ~ normal(psi1,sg2);
    psi1G_u[g] ~ normal(psi1_u,sg3);
    muSZG[g] ~ normal(muSZ,sg4);
    muSZG_u[g] ~ normal(muSZ_u,sg5);
    lgtLMG[g] ~ normal(lgtLM,sg6);
  }
  tauB ~ cauchy(0,1) ;      // precision param for diet comp variation across bouts
  tauG ~ cauchy(0,1) ;      // precision param for diet comp variation across groups
  phi1 ~ cauchy(0,.5) ;     // Intercept of log Cons Rate, by prey
  phi2 ~ cauchy(0,.1) ;     // effect of size on log Cons rate
  psi1 ~ cauchy(0,1) ;      // Intercept of log HT, by prey
  psi2 ~ cauchy(0,.1) ;     // effect of size on log HT rate
  psi1_u ~ cauchy(0,1) ;    // Intercept of log HT, UN-ID prey
  psi2_u ~ cauchy(0,.1) ;   // effect of size on log HTm UN-ID prey
  muSZ ~ normal(4,1);       // log-mean size by prey type
  muSZ_u ~ normal(4,1);     // log-mean size of UN-ID prey
  maxPunid ~ beta(1,1);     // max prob un-ID, for prey overlapping UN-ID in size & HT
  sigSZ ~ cauchy(0,.5);     // variation in prey size across bouts, by prey type
  sigHT ~ cauchy(0,.5);     // variation in prey HT across bouts, by prey type
  sigCR ~ cauchy(0,.5);     // variation in prey CR across bouts, by prey type
  sigSZ_u ~ cauchy(0,.5);   // variation in prey size across bouts, unid prey
  sigHT_u ~ cauchy(0,.5);   // variation in prey HT across bouts, unid prey
  lgtLM ~ normal(1,2);      // mean logit(lambda), dive success rate, by prey
  sigLM ~ cauchy(0,1.5);    // variation in logit(lambda), dive success rate  
  sg1 ~ cauchy(0,.05);      // variation across groups, phi1 (Crate intercept)
  sg2 ~ cauchy(0,.05);      // variation across groups, psi1 (HT intercept)
  sg3 ~ cauchy(0,.05);      // variation across groups, psi1 un-ID
  sg4 ~ cauchy(0,.05);      // variation across groups, mn log size
  sg5 ~ cauchy(0,.05);      // variation across groups, mn log size un-ID
  sg6 ~ cauchy(0,.05);      // variation across groups, dive success rate
}
// Section 5. Derived parameters and statistics   
generated quantities {
  vector[Km1] upsilon ;     // relative contribution to unknown prey dives
  array[Ngrp] vector[Km1] upsilonG ;
  vector<lower=0,upper=1>[Km1] Omega ;
  vector[Km1] SZ ;
  array[Ngrp] vector[Km1] SZg ;
  vector[Km1] logSZ ;
  array[Ngrp] vector[Km1] logSZg ;
  real SZ_u ;
  array[Ngrp] real SZg_u ;
  vector[Km1] HT ;
  array[Ngrp] vector[Km1] HTg ;
  real HT_u ;
  array[Ngrp] real HTg_u ;
  vector[Km1] CR ;
  array[Ngrp] vector[Km1] CRg ;
  vector[Km1] Pi ;
  array[Ngrp] vector[Km1] PiG ;
  vector[Km1] ER ;
  array[Ngrp] vector[Km1] ERg ;
  real CRmn ;
  array[Ngrp] real CRgmn ;
  real ERmn ;
  array[Ngrp] real ERgmn ;
  vector[Km1] LM ;
  array[Ngrp] vector[Km1] LMg ;
  real LMmn ;
  array[Ngrp] real LMgmn ;
  // Mean Size (mm) by prey type (adjust for log-normal)
  logSZ = muSZ + square(sigSZ)/2 ;
  SZ = exp(muSZ + square(sigSZ)/2 );
  SZ_u = exp(muSZ_u + square(sigSZ_u)/2 );
  for(g in 1:Ngrp){
  	logSZg[g] = muSZG[g] + square(sigSZ)/2 ;
    SZg[g] = exp(muSZG[g] + square(sigSZ)/2 );
    SZg_u[g] = exp(muSZG_u[g] + square(sigSZ_u)/2 );
  }
  for (j in 1:Km1){
    // Mean Cons rate (CR, g/min) by prey type, adjusted for mean prey size and lognormal dist
    CR[j] = fmin(CR_max[j], exp(phi1[j] + phi2[j] * (2.5*logSZ[j]-7) + square(fmin(1.5,sigCR[j]))/2 + lgSz_adj[j]) );
    // Mean HT/itm, by prey type, adjusted for mean prey size and lognormal dist
    HT[j] = fmin(600, exp(psi1[j] + psi2[j] * (2.5*logSZ[j]-7) + square(sigHT[j])/2)) ;
    LM[j] = inv_logit(lgtLM[j]) ;
    Omega[j] = 0 ;
    for(g in 1:Ngrp){
      Omega[j] = Omega[j] + ((1.0 / Ngrp) * OmegaG[g][j]) ;
      CRg[g][j] = fmin(CR_max[j], exp(phi1G[g][j] + phi2[j] * (2.5*logSZg[g][j]-7) + square(fmin(1.5,sigCR[j]))/2 + lgSz_adj[j]) );
      HTg[g][j] = fmin(600, exp(psi1G[g][j] + psi2[j] * (2.5*logSZg[g][j]-7) + square(sigHT[j])/2)) ;
      LMg[g][j] = inv_logit(lgtLMG[g][j]) ;
    }
  }
  upsilon = (eta .* (1 - Omega)) / sum(eta .* (1 - Omega));
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
  // Repeat above stats, by group
  for(g in 1:Ngrp){    
    HTg_u[g] = exp(psi1G_u[g] + psi2_u * (2.5*muSZG_u[g]-7) + square(sigHT_u)/2);
    ERg[g] = Cal_dens .* CRg[g] ;
    PiG[g] = (etaG[g] .* CRg[g]) / sum(etaG[g] .* CRg[g]) ;
    CRgmn[g] = sum(etaG[g] .* CRg[g]) ;
    ERgmn[g] = sum(etaG[g] .* CRg[g] .* Cal_dens) ;
    LMgmn[g] = sum(etaG[g] .* LMg[g]) ;
    upsilonG[g] = (1000 .* PiG[g] ./ MnMass) .* (1 - OmegaG[g]) ; 
    upsilonG[g] = upsilonG[g] ./ sum(upsilonG[g]) ;
  }
}
