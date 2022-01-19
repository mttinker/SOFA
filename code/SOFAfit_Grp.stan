// SOFAfit
// This Stan code fits a model of foraging behavior to observational data 
//  on sea otter foraging  [Refer to "SOFA_summary" for details] 
//  - this is complex version for "Grouped" data
//
// Section 1. Data inputs for model
data {
  int<lower=1> Nbouts;            // Number of feeding bouts (for effort allocation observations)
  int<lower=1> Ngrp ;             // Number of feeding bouts (for effort allocation observations)
  int<lower=1> K;                 // Number of prey types (1 to (K-1) for ID prey, UN-ID prey = K)
  int<lower=1> Km1;               // "K - 1", number of prey types excluding UN-ID prey  
  int<lower=0> EffortP[Nbouts,K]; // Minutes allocated to each prey type (incl. un-id) by bout  
  int<lower=1> GrpE[Nbouts];      // Group ID for each bout (for effort allocations)  
  int<lower=1> NSz;               // Number of prey Size obs
  int<lower=1> NHt;               // Number of prey Handling time (HT) obs  
  int<lower=1> NCR;               // Number of prey Consumption rate obs
  int<lower=1> NLm;               // Number of dive success rate obs (lambda)
  int<lower=1> NU[2];             // Number of un-ID prey obs (2 different param types, FOR NOW)
  real<lower=0> SZmnU[NU[1]];     // Array of mean prey size vaues (mm) for UN-ID prey, by bout  
  real<lower=0> HTmnU[NU[2]];     // Array of mean handling times (sec), for UN-ID prey, by bout 
  int<lower=1,upper=Ngrp> Sg_u[NU[1]]; // Array of group id values for size obs of UN-ID
  int<lower=1,upper=Ngrp> Hg_u[NU[2]]; // Array of group id values for HT obs of UN-ID  
  int<lower=1,upper=Km1> Sp[NSz]; // Array of prey id values for size obs
  int<lower=1,upper=Km1> Hp[NHt]; // Array of prey id values for HT obs  
  int<lower=1,upper=Km1> Cp[NCR]; // Array of prey id values for Cons Rate obs
  int<lower=1,upper=Ngrp> Sg[NSz]; // Array of group id values for size obs
  int<lower=1,upper=Ngrp> Hg[NHt]; // Array of group id values for HT obs  
  int<lower=1,upper=Ngrp> Cg[NCR]; // Array of group id values for Cons Rate obs
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
  int<lower=1,upper=Ngrp> Lg[NLm];// Array of group id values for Dive success rate
  real LMlg[NLm];                 // Array of logit(lambda), proportion of succesful dives by prey, by bout  
  vector[NLm] Lss;                // Array of sample sizes (n obs per bout) for logit lambda observations   
  vector<lower=0>[Km1] CR_max ;   // maximum possible mean CR per prey type
}
transformed data{
	vector[Km1] logMass_mn;
	logMass_mn = -1 * square(logMass_sg)/2 ;
}
// Section 2. The parameters to be estimated 
parameters {
  simplex[Km1] eta ;             // mean proportional effort allocation by prey type (if all ID'd)
  simplex[Km1] etaG[Ngrp] ;      // mean proportional effort allocation by prey type, by group
  real<lower=0,upper=20> tauB[Ngrp] ;     // relative precision of diet composition across bouts
  real<lower=0,upper=10> tauG ;           // relative precision of diet composition across groups
  simplex[K] theta[Nbouts] ;     // vectors of bout-specific prey-allocation probs for multinomial  
  vector<lower=0,upper=6>[Km1] muSZ ;    // mean of log(Size_mm) of each prey type
  vector<lower=0,upper=6>[Km1] muSZG[Ngrp] ;  // mean of log(Size_mm) of each prey type, by group 
  real<lower=0,upper=6> muSZ_u ;         // mean of log(Size_mm) of UN-ID
  real<lower=0,upper=6> muSZG_u[Ngrp] ;         // mean of log(Size_mm) of UN-ID, by group
  vector[Km1] lgtLM ;            // mean logit(lambda) (dive success rate) for each prey
  vector[Km1] lgtLMG[Ngrp] ;     // mean logit(lambda) (dive success rate) for each prey, by group
  vector<lower=0,upper=20>[Km1] sigLM ;   // stdev of logit(lambda) of each prey 
  vector<lower=0,upper=20>[Km1] sigSZ ;   // stdev of log(Size_cm) of each prey 
  vector<lower=0,upper=20>[Km1] sigHT ;   // stdev of log(HT per item) of each prey 
  vector<lower=0,upper=20>[Km1] sigCR ;   // stdev of log(HT per item) of each prey 
  real<lower=0,upper=20> sigSZ_u ;        // stdev of log(Size_mm) of UN-ID
  real<lower=0,upper=20> sigHT_u ;        // stdev of log(HT per item) of UN-ID    
  vector<lower=0>[Km1] phi1 ;    // intercept of CRate v Size fxn by prey 
  vector<lower=0>[Km1] phi1G[Ngrp] ;    // intercept of CRate v Size fxn by prey 
  vector<lower=0>[Km1] phi2 ;    // slope of CRate v Size fxn by prey 
  vector<lower=0>[Km1] psi1 ;    // intercept of HT v Size fxn by prey 
  vector<lower=0>[Km1] psi1G[Ngrp] ;    // intercept of HT v Size fxn by prey 
  vector<lower=0>[Km1] psi2 ;    // slope of HT v Size fxn by prey   
  real<lower=0> psi1_u ;         // intercept of HT v Size fxn for UN-ID prey 
  real<lower=0> psi1G_u[Ngrp] ;         // intercept of HT v Size fxn for UN-ID prey 
  real<lower=0> psi2_u ;         //  slope of HT v Size fxn for UN-ID prey     
  real<lower=0,upper=1> maxPunid;// max prob un-id, if size/HT dist.s overlap 100% with UN-ID
  vector<lower=0>[Km1] Cal_dens; // Caloric density (kcal/g) by prey type
  vector[Km1] lgSz_adj; // Uncertainty adjustment for size-biomass fxn
  real<lower=0,upper=5> sg1 ;        // stdev of hierarchical param 
  real<lower=0,upper=5> sg2 ;        // stdev of hierarchical param 
  real<lower=0,upper=5> sg3 ;        // stdev of hierarchical param 
  real<lower=0,upper=5> sg4 ;        // stdev of hierarchical param 
  real<lower=0,upper=5> sg5 ;        // stdev of hierarchical param 
  real<lower=0,upper=5> sg6 ;        // stdev of hierarchical param 
}
// Section 3. Derived (transformed) parameters
transformed parameters {
  vector<lower=0,upper=1>[Km1] OmegaG[Ngrp] ;// prey-specific probability of positive ID
  vector<lower=0>[K] alpha[Ngrp] ;        // vector of dirichlet params: relative prey freq (inc. UN-ID)
  //
  for (g in 1:Ngrp){  
    real muHT_u ;                     // mean of log(HT per item) of UN-ID   
    vector[Km1] muHT ;                // mean of log(HT per item) of each prey  
    // Compute mean expected HT at avg size for un-id prey
    muHT_u = psi1G_u[g] + psi2_u * (2.5 * muSZG_u[g] - 7);  
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
      muHT[j] = psi1G[g][j] + psi2[j] * (2.5 * muSZG[g][j] - 7) ;  
      Dist_SZ = 0.25 * log(0.25 * (sigSZ[j]^2 / sigSZ_u^2 + sigSZ_u^2/sigSZ[j]^2 + 2))
        + 0.25 * (((muSZ[j] - muSZ_u)^2) / (sigSZ[j]^2 + sigSZ_u^2)) ;
      OV_SZ = exp(-Dist_SZ) ;
      Dist_HT = 0.25 * log(0.25 * (sigHT[j]^2 / sigHT_u^2 + sigHT_u^2/sigHT[j]^2 + 2))
        + 0.25 * (((muHT[j] - muHT_u)^2) / (sigHT[j]^2 + sigHT_u^2)) ;
      OV_HT = exp(-Dist_HT) ;    
      OmegaG[g][j] = 1 - maxPunid * (OV_SZ * OV_HT) ;
    } 
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
  // Observed consumption rate (g/min) by prey type, random samples
  // NOTE: Expected log attributes, given that the log median of means of lognormal samples of size n: 
  //           mu + (n*sg^2 - sg^2)/(2*n) 
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
  for (g in 1:Ngrp){
    etaG[g] ~ dirichlet(eta * 20 * tauG);
    phi1G[g] ~ normal(phi1,sg1);
    psi1G[g] ~ normal(psi1,sg2);
    psi1G_u[g] ~ normal(psi1_u,sg3);
    muSZG[g] ~ normal(muSZ,sg4);
    muSZG_u[g] ~ normal(muSZ_u,sg5);
    lgtLMG[g] ~ normal(lgtLM,sg6);
  }
  tauB ~ cauchy(0,1) ;  // precision param for diet comp variation across bouts
  tauG ~ cauchy(0,.5) ;  // precision param for diet comp variation across groups
  phi1 ~ cauchy(0,.5) ;   // Intercept of log Cons Rate, by prey
  phi2 ~ cauchy(0,.5) ;   // effect of size on log Cons rate
  psi1 ~ cauchy(0,1) ;   // Intercept of log HT, by prey
  psi2 ~ cauchy(0,.5) ;   // effect of size on log HT rate
  psi1_u ~ cauchy(0,1) ; // Intercept of log HT, UN-ID prey
  psi2_u ~ cauchy(0,.5) ; // effect of size on log HTm UN-ID prey
  muSZ ~ normal(4,1);      // log-mean size by prey type
  muSZ_u ~ normal(4,1);    // log-mean size of UN-ID prey
  maxPunid ~ beta(3,1);    // max prob un-ID, for prey overlapping UN-ID in size & HT
  sigSZ ~ cauchy(0,.5);    // variation in prey size across bouts, by prey type
  sigHT ~ cauchy(0,.5);    // variation in prey HT across bouts, by prey type
  sigCR ~ cauchy(0,.5);    // variation in prey CR across bouts, by prey type
  sigSZ_u ~ cauchy(0,.5);  // variation in prey size across bouts, unid prey
  sigHT_u ~ cauchy(0,.5);  // variation in prey HT across bouts, unid prey
  lgtLM ~ normal(1,2); // mean logit(lambda), dive success rate, by prey
  sigLM ~ cauchy(0,1.5);   // variation in logit(lambda), dive success rate  
  sg1 ~ cauchy(0,.25);
  sg2 ~ cauchy(0,.25);
  sg3 ~ cauchy(0,.25);
  sg4 ~ cauchy(0,.25);
  sg5 ~ cauchy(0,.25);
  sg6 ~ cauchy(0,.25);
}
// Section 5. Derived parameters and statistics 
generated quantities {
  vector<lower=0,upper=1>[Km1] Omega ;
  vector[Km1] SZ ;
  vector[Km1] SZg[Ngrp] ;
  vector[Km1] logSZ ;
  vector[Km1] logSZg[Ngrp] ;
  real SZ_u ;
  real SZg_u[Ngrp] ;
  vector[Km1] HT ;
  vector[Km1] HTg[Ngrp] ;
  real HT_u ;
  real HTg_u[Ngrp] ;
  vector[Km1] CR ;
  vector[Km1] CRg[Ngrp] ;
  vector[Km1] Pi ;
  vector[Km1] PiG[Ngrp] ;
  vector[Km1] ER ;
  vector[Km1] ERg[Ngrp] ;
  real CRmn ;
  real CRgmn[Ngrp] ;
  real ERmn ;
  real ERgmn[Ngrp] ;
  vector[Km1] LM ;
  vector[Km1] LMg[Ngrp] ;
  real LMmn ;
  real LMgmn[Ngrp] ;
  // Mean Size (mm) by prey type (adjust for log-normal)
  logSZ = muSZ + square(sigSZ)/2 ;
  SZ = exp(muSZ + square(sigSZ)/2 );
  SZ_u = exp(muSZ_u +  square(sigSZ_u)/2 );
  for(g in 1:Ngrp){
  	logSZg[g] = muSZG[g] + square(sigSZ)/2  ;
    SZg[g] = exp(muSZG[g] +  square(sigSZ)/2 );
    SZg_u[g] = exp(muSZG_u[g] +  square(sigSZ_u)/2 );
  }
  for (j in 1:Km1){
    // Mean Cons rate (CR, g/min) by prey type, adjusted for mean prey size and lognormal dist
    CR[j] = fmin(CR_max[j], exp(phi1[j] + phi2[j] * (2.5*logSZ[j]-7) + square(sigCR[j])/2 + lgSz_adj[j]) );
    // Mean HT/itm, by prey type, adjusted for mean prey size and lognormal dist
    HT[j] = fmin(600, exp(psi1[j] + psi2[j] * (2.5*logSZ[j]-7) + square(sigHT[j])/2)) ;
    LM[j] = inv_logit(lgtLM[j]) ;
    Omega[j] = 0 ;
    for(g in 1:Ngrp){
      Omega[j] = Omega[j] + ((1.0 / Ngrp) * OmegaG[g][j]) ;
      CRg[g][j] = fmin(CR_max[j], exp(phi1G[g][j] + phi2[j] * (2.5*logSZg[g][j]-7) + square(sigCR[j])/2 + lgSz_adj[j]) );
      HTg[g][j] = fmin(600, exp(psi1G[g][j] + psi2[j] * (2.5*logSZg[g][j]-7) + square(sigHT[j])/2)) ;
      LMg[g][j] = inv_logit(lgtLMG[g][j]) ;
    }
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
  // Repeat above stats, by group
  for(g in 1:Ngrp){
    HTg_u[g] = exp(psi1G_u[g] + psi2_u * (2.5*muSZG_u[g]-7) + square(sigHT_u)/2);
    ERg[g] = Cal_dens .* CRg[g] ;
    PiG[g] = (etaG[g] .* CRg[g]) / sum(etaG[g] .* CRg[g]) ;
    CRgmn[g] = sum(etaG[g] .* CRg[g]) ;
    ERgmn[g] = sum(etaG[g] .* CRg[g] .* Cal_dens) ;
    LMgmn[g] = sum(etaG[g] .* LMg[g]) ;
  }
}
