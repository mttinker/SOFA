// SOFAfit
// This Stan code fits a model of foraging behavior
//  to observational data on sea otter foraging
//  [Refer to "SOFA_summary" for details] 
//
// Section 1. Data inputs for model
data {
  int<lower=1> Nbouts;            // Number of feeding bouts (for effort allocation observations)
  int<lower=1> K;                 // Number of prey types (1 to (K-1) for ID prey, UN-ID prey = K)
  int<lower=1> Km1;               // "K - 1", number of prey types excluding UN-ID prey  
  int<lower=0> EffortP[Nbouts,K]; // Minutes allocated to each prey type (incl. un-id) by bout  
  int<lower=1> NSz;               // Number of prey Size obs
  int<lower=1> NNi;               // Number of prey Nitem obs
  int<lower=1> NHt;               // Number of prey Handling time (HT) obs
  int<lower=1> NLm;               // Number of prey success rate (Lambda) obs
  int<lower=1> NCR;               // Number of prey Consumption rate obs
  int<lower=1> NU[4];             // Number of un-ID prey obs (4 different param types)
  real<lower=0> SZmnU[NU[1]];     // Array of mean prey size vaues (cm) for UN-ID prey, by bout  
  real<lower=0> NImnU[NU[2]];     // Array of mean number items per dive, for UN-ID prey, by bout  
  real<lower=0> HTmnU[NU[3]];     // Array of mean handling times (sec), for UN-ID prey, by bout   
  real LMlgU[NU[4]];              // Array of mean logit success rates, for UN-ID prey, by bout    
  int<lower=1,upper=Km1> Sp[NSz]; // Array of prey id values for size obs
  int<lower=1,upper=Km1> Np[NNi]; // Array of prey id values for Nitm obs 
  int<lower=1,upper=Km1> Hp[NHt]; // Array of prey id values for HT obs
  int<lower=1,upper=Km1> Lp[NLm]; // Array of prey id values for Lambda obs
  int<lower=1,upper=Km1> Cp[NCR]; // Array of prey id values for Cons Rate obs
  real<lower=0> SZmn[NSz];        // Array of mean prey size vaues (cm) by prey, by bout  
  real<lower=0> NImn[NNi];        // Array of mean number items per dive, by prey, by bout  
  real<lower=0> HTmn[NHt];        // Array of mean handling times (sec), by prey, by bout   
  real LMlg[NLm];                 // Array of mean logit success rates, by prey, by bout    
  real<lower=0> CRate[NCR];       // Array of prey-specific consumption rates (CR) by bout (g/min)
  real<lower=0> Csz[NCR];         // Array of prey size (cm) covariate for CR obs
  real<lower=0> Css[NCR];         // Array of sample size covariate for CR obs
  vector<lower=0>[Km1] Cal_dens;  // Array of caloric density values (kcal/g) by prey type      
}
// Section 2. The parameters to be estimated 
parameters {
  // Annual proportions of pups in areas 1 and 2, predicted
  simplex[Km1] eta ;           // mean proportional effort allocation by prey type (if all ID'd)
  real<lower=0> tau ;          // relative precision of diet composition across bouts
  simplex[K] theta[Nbouts] ;   // vectors of bout-specific prey-allocation probs for multinomial  
  vector[Km1] muSZ ;           // mean of log(Size_cm) of each prey
  vector[Km1] muNI ;           // mean of log(N_items) of each prey
  vector[Km1] muHT ;           // mean of log(HT per item) of each prey
  vector[Km1] mnlgLM ;         // mean of logit(Lambda) of each prey
  real<lower=0> sigSZ[K] ;     // stdev of log(Size_cm) of each prey (incl. un-id)
  real<lower=0> sigNI[K] ;     // stdev of log(N_items) of each prey (incl. un-id)
  real<lower=0> sigHT[K] ;     // stdev of log(HT per item) of each prey (incl. un-id)
  real<lower=0> sigLM[K] ;     // stdev of logit(Lambda) of each prey (incl. un-id)
  real<lower=0> phi1[Km1] ;    // intercept of mean cons rate fxn by prey 
  real<lower=0> phi2[Km1] ;    // size effect on mean cons rate fxn by prey 
  real<lower=0> psi1[Km1] ;    // base inv. scale param fxn for gamma dist CR, by prey
  real<lower=0> psi2 ;         // expon for inv.scal param fxn for gamma dist CR
  real beta0 ;                 // logit param intercept for P (prob of pos id by prey type)
  real beta1 ;                 // logit param effect of pr. size on P (prob of pos id by prey type)
  real beta2 ;                 // logit param effect of pr. HT on P (prob of pos id by prey type)
  real epsP[Km1] ;             // Unexplained deviations in prob of ID by prey type (logit)
  real sigP  ;                 // Variance across prey types in logit-P (prob of pos ID)
}
// Section 3. Derived (transformed) parameters
transformed parameters {
  vector<lower=0>[Km1] Pid ;        // prey-specific probability of positive ID
  vector<lower=0>[Km1] q ;        // vector of prey-specific contributions to un-ID, by bout  
  vector<lower=0>[Km1] alphaP ;   // vector of dirichlet params: relative prey freq (w/o UN-ID)
  vector<lower=0>[K] alpha ;      // vector of dirichlet params: relative prey freq (inc. UN-ID)
  real muSZ_u ;                   // mean of log(Size_cm) of UN-ID
  real muNI_u ;                   // mean of log(N_items) of UN-ID 
  real muHT_u ;                   // mean of log(HT per item) of UN-ID  
  real mnlgLM_u ;                 // mean of logit(Lambda) of UN-ID 
  //
  // Calculate alpha, dirichlet params for bout-specific prey allocation probs  
  alphaP = eta * (tau * Km1) ;  
  // Loop through prey to calculate prob of positive ID (P) and contribution to un-ID (q)
  for (j in 1:Km1){
    Pid[j] = inv_logit(beta0 + beta1 * muSZ[j] + beta2 * muHT[j] + epsP[j]) ;  
    alpha[j] = alphaP[j] * (Pid[j]) ;
  }  
  alpha[K] = sum(alphaP .* (1 - Pid)) ;
  q = (alphaP .* (1 - Pid) ) / sum(alphaP .* (1 - Pid)) ;
  // Calculate un-ID prey params based on contributions of prey types to un-ID class
  muSZ_u = log(sum(q .* exp(muSZ) )) ;
  muNI_u = log(sum(q .* exp(muNI) )) ;
  muHT_u = log(sum(q .* exp(muHT) )) ;
  mnlgLM_u = logit(sum(q .* inv_logit(mnlgLM) )) ;
}
// Section 4. Estimating model parameters (drawing from probability distributions)
model {
  // A) Observed nodes:
  // Allocation of effort by prey type (proportional and # minutes of each feeding bout)
  for (i in 1:Nbouts){
    theta[i] ~ dirichlet(alpha);  
    EffortP[i,] ~ multinomial(theta[i]) ;
  }
  // Observed consumption rates by prey type, random samples
  for (i in 1:NCR){
    CRate[i] ~ gamma( (phi1[Cp[i]]+phi2[Cp[i]]*Csz[Cp[i]])*(psi1[Cp[i]]*Css[i]^psi2), psi1[Cp[i]]*Css[i]^psi2);
  }
  // Mean prey size, random samples, ID prey
  for (i in 1:NSz){
    SZmn[i] ~ lognormal(muSZ[Sp[i]],sigSZ[Sp[i]]) ;
  }
  // Mean prey Nitems, random samples, ID prey
  for (i in 1:NNi){
    NImn[i] ~ lognormal(muNI[Np[i]],sigNI[Np[i]]) ;
  }
  // Handling time, random samples, ID prey
  for (i in 1:NHt){
    HTmn[i] ~ lognormal(muHT[Hp[i]],sigHT[Hp[i]]) ;
  }  
  // Lamda (success rate), logit-transformed samples, ID prey
  for (i in 1:NLm){
    LMlg[i] ~ normal(mnlgLM[Lp[i]],sigLM[Lp[i]]) ;
  }  
  // Observed params for UN-ID prey, random samples
  SZmnU ~ lognormal(muSZ_u,sigSZ[K]) ;
  NImnU ~ lognormal(muNI_u,sigNI[K]) ;
  HTmnU ~ lognormal(muHT_u,sigHT[K]) ;
  LMlgU ~ normal(mnlgLM_u,sigLM[K]) ;
  //
  // B) Prior distributions for model parameters:
  tau ~ cauchy(0,2.5) ;  // precision param for diet comp variaiton across bouts
  epsP ~ normal(0, sigP);// Unexplained deviation in prob of ID by prey type (logit)
  phi1 ~ gamma(1,.05) ;  // Intercept for Cons Rate, by prey: positive real ~ 0-50
  phi2 ~ cauchy(0,2.5) ; // effect of size on Cons rate
  psi1 ~ cauchy(0,2.5) ; // intercept for inv.scale param, gamma dist Cons Rate
  psi2 ~ cauchy(0,2.5) ; // sample size effect for inv.scale param, gamma dist Cons Rate
  beta0 ~ cauchy(0,2.5); // intercept, logit fxn for prey-specific prob of ID
  beta1 ~ cauchy(0,2.5); // effect sz, logit fxn for prey-specific prob of ID
  beta2 ~ cauchy(0,2.5); // effect ht, logit fxn for prey-specific prob of ID
  muSZ ~ cauchy(1.5,2.5);// log-mean size by prey type
  muNI ~ cauchy(1,2.5);  // log-mean N items by prey type
  muHT ~ cauchy(3,2.5);  // log-mean HT by prey type
  mnlgLM ~ normal(1,3);  // logit of prey-specific success rates
  sigSZ ~ cauchy(0,2.5); // variation in prey size across bouts, by prey type
  sigNI ~ cauchy(0,2.5); // variation in prey Nitem across bouts, by prey type
  sigHT ~ cauchy(0,2.5); // variation in prey HT across bouts, by prey type
  sigLM ~ cauchy(0,2.5); // variation in prey success rate across bouts, by prey type
  sigP ~ cauchy(0,2.5);  // variation in prob of ID accross prey types (logit)
}
// Section 5. Derived parameters and statistics 
generated quantities {
  vector[Km1] SZ ;
  vector[Km1] NI ;
  vector[Km1] HT ;
  vector[Km1] LM ;
  vector[Km1] CR ;
  vector[Km1] PD ;
  vector[Km1] ER ;
  real CRmn ;
  real ERmn ;
  // Mean Size (cm), Nitm/dv, HT/itm and Lambda, by prey type (adjust for log-normal)
  SZ = exp(muSZ +  square(to_vector(sigSZ[1:Km1]))/2 );
  NI = exp(muNI +  square(to_vector(sigNI[1:Km1]))/2 );
  HT = exp(muHT +  square(to_vector(sigHT[1:Km1]))/2 );
  LM = inv_logit(mnlgLM) ;
  // Mean Cons rate (CR) by prey (g/min), adjusted for mean prey size (in cm)
  CR = to_vector(phi1) + to_vector(phi2) .* SZ  ;
  // Mean Energy intake rate (ER) by prey (kcal/min)
  ER = Cal_dens .* CR ;
  // Proportional contribution of each prey type to diet (PD, biomass consumed): (eta.*CR)/sum(eta .* CR)
  PD = (eta .* CR) / sum(eta .* CR) ;
  // Overall mean consumption rate (CRmn) given effort allocation to each prey: sum(eta .* CR)
  CRmn = sum(eta .* CR) ;
  // Overall mean Energy Intake Rate (ERmn) given effort allocation: sum(eta .* CR .* Cal_dens)
  ERmn = sum(eta .* CR .* Cal_dens) ;
  // 
}
