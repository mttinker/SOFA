# Simulation of foraging data
# 
require(gtools)
require(stats)
require(fitdistrplus)
require(parallel)
require(rstan)
require(bayesplot)
require(lattice)
require(coda)
require(readxl)
# Number bouts
Nbouts = 500
# Number of prey captures per bout
Npr_p_bt = 20
# Set up some variables ------------------------------------------------------
# Prey size key
PrSzC = factor(c("1a","1b","1c","1",
                "2a","2b","2c","2",
                "3a","3b","3c","3",
                "4a","4b","4c","4"))
PrszR = c(0.25,.5,.75,.5); PrszR1 = PrszR
PrszR = c(PrszR,PrszR1+1);PrszR = c(PrszR,PrszR1+2);PrszR = c(PrszR,PrszR1+3)
PreySzLst = data.frame(Szcode = PrSzC,SzRatio = PrszR)
PreySzLst$SzClass = rep(c(1,2,3,4),each=4)
PreySzLst$SzQual = rep(c("a","b","c",""),4)
PreySzLst = PreySzLst[order(PreySzLst$SzRatio),]; rownames(PreySzLst) = seq(1,16)
Pawsz = 47.5 # in mm, for prey size to mass conversions 
# Prey params  
NPtypes = 7
prlist = c("can","urc","mus","kcr","sna","was","clm")
# Biomass params (for power function): Mass = PFa[i]*(Pawsz*PreySzLst$SzRatio[j])^PFb[i]
# NOTE : replace this with function that fits these params from prey data
PFa = c(.0053, 0.0005, 0.0034, 0.00882, 0.00187, 0.00023, 0.00046)
PFb = c(2.195, 2.9035, 2.0464, 2.15248, 2.49975, 2.82750, 2.70797)
Cal_dns_mn = c(1.22, 0.393, 0.497, 0.733, 1.021, 0.505, .535)
Cal_dns_sg = c(0.44, 0.05, 0.02, 0.188, 0.15, 0.05, 0.05)
Preytypes = data.frame(Prey=factor(prlist,levels =prlist))
# Relative freq of capture for simulations:
eta_CAP = c(.05,.2,.23,.12,.18,.08,.14)
tau_CAP = 2
Drchvec = eta_CAP*(tau_CAP*NPtypes)
# To set random prey type capture probs per bout: Pprey = rdirichlet(1,Drchvec)
#   (then use multinomial to assigh prey)
# Size class params (gamma fxn): ceiling(rgamma(1,shape=Spar1*Spar2, rate=Spar2))
Spar1 = c(9,3.5,3,4,1.5,7.5,3.2)
Spar2 = c(2,3.0,3,2,5.0,2.0,2)
# N-items params (gamma fxn): ceiling(rgamma(1,shape=Npar1*Npar2, rate=npar2))
Npar1 = c(0.6,2.5,3,1,9,0.8,2)
Npar2 = c(7.0,2,5,2,1,5.0,2)
# HT params for adjusting baseline fxn
Hadj = c(1.4,.9,.3, 1.2, .4,.8,.3)
# Note: HT(s) = rpois((55*PreySzLst$SzRatio^1.5)*Hadj)
# Dive succes rates (determines number of unsucc dives prior to successful dive)
Lmbda = c(.6,.93,.97,.88,.99,.55,.8)
# Note: use neg binomial to get number of preceeding unsucc dives given prob:
# NumUnsucc = rnbinom(1,1,Lmbda[i])
# Second prey type can be captured on single dive (must be another with MultP==1)?
MultP = c(0,1,1,1,0,0,1)
PrbSc = 0.35 # for prey types allowing "second types", prob of a second prey type
# Prob of positive ID prey, logit fxn:
#   probID_Prey = inv.logit(B0 + B1*log(Sz_cm) + B2*log(HT) + B3*log(Sz_cm)*log(HT))
B0 = -1.1 # Intercept
B1 = .2 # Effect of size (cm)  
B2 = .05 # Effect of HT (in decimal minutes)
B3 = .4 # Effect of size*HT interaction
# B3 = c(1,0,0,0.5,0.5,1,0) # Prey-specific modifiers of ID prob
#
# Prob of ID Nitems: probID_N = inv.logit(6 - 3*Nitm^0.25)
# Prob of ID Sz:   inv.logit(1 + .2*Sz_cm)
#
# Prob of carry over
CryOv_P = c(.2,0,0,0,0,.15,0)
CryOv_Psz = c(0,.1,.5,1)
# Proportion of prey type lost, possible values (only for sizecl >2):
prplst = c(0,.5,.25) # Nlst = prplst[pmin(3,ceiling(rgamma(1,shape=1, rate=1.5)))]
prplstC = c(0, 0, 1, 1) # Size class-based modifiers for proportion lost vals
#
# Initialize variables for foraging data
Bout = character()
Divenum = numeric()
DT = numeric()
ST = numeric()
Success = character()
SuccessV = numeric()
Prey = character()
PreyV = numeric()
SzClss = numeric()
Size = numeric()
Qual = numeric()
Sz_cm = numeric()
Sz_rt = numeric()
Mss_est = numeric()
Nitem = numeric()
Ppnlost = numeric()
Ncrct = numeric()
HT = numeric()
#
# Simulate bouts----------------------------------------------------------------------
rc = 0
for (b in 1:Nbouts){
  # Prefprey = rmultinom(1,1,Drchvec)*8
  Pprey = rdirichlet(1,Drchvec)
  Preycaps = rmultinom(Npr_p_bt,1,Pprey)
  dvn = 0
  btid = paste0("Bout_",b)
  prn = 0
  while(prn < Npr_p_bt){
    rc = rc + 1
    dvn = dvn + 1
    Bout[rc] = btid
    Divenum[rc] = dvn
    DT[rc] = rpois(1,60)
    prn = prn + 1
    p = which(Preycaps[,prn]==1)
    prey = as.character(Preytypes[p,1])
    NumUnsucc = rnbinom(1,1,Lmbda[p])
    if(NumUnsucc==1){
      Success[rc] = "n"
      SuccessV[rc] = 0
      ST[rc] = rpois(1,5)+5
      rc = rc + 1
      dvn = dvn + 1
      Bout[rc] = btid
      Divenum[rc] = dvn
      DT[rc] = rpois(1,60)
    }else if(NumUnsucc>1){
      ST[rc] = rpois(1,5)+5
      Success[rc] = "n"
      SuccessV[rc] = 0
      for(i in 1:(NumUnsucc-1)){
        rc = rc + 1
        dvn = dvn + 1
        Bout[rc] = btid
        Divenum[rc] = dvn
        DT[rc] = rpois(1,60)
        ST[rc] = rpois(1,5)+5
        Success[rc] = "n"
        SuccessV[rc] = 0
      }
      rc = rc + 1
      dvn = dvn + 1
      Bout[rc] = btid
      Divenum[rc] = dvn
      DT[rc] = rpois(1,60)
    }
    Prey[rc] = prey
    PreyV[rc] = p
    S = min(16,ceiling(rgamma(1,shape=Spar1[p]*Spar2[p], rate=Spar2[p])))
    Size[rc] = PreySzLst$SzClass[S]
    Qual[rc] = PreySzLst$SzQual[S]   
    SzClss[rc] = S
    Sz_cm[rc] = 0.1*Pawsz*PreySzLst$SzRatio[S]
    Sz_rt[rc] = PreySzLst$SzRatio[S]
    Nitem[rc] = ceiling(rgamma(1,shape=Npar1[p]*Npar2[p], rate=Npar2[p]))
    CrryOv = rbinom(1,1,CryOv_P[p]*CryOv_Psz[Size[rc]])
    if (CrryOv == 1 & Nitem[rc]==1){
      nco = rpois(1,1)+1
      Success[rc] = "y"
      SuccessV[rc] = 1        
      Ppnlost[rc] = 0
      Ncrct[rc] = Nitem[rc] - (Ppnlost[rc]*Nitem[rc])
      Mss_est[rc] = (PFa[p]*(Pawsz*Sz_rt[rc])^PFb[p])*Ncrct[rc]  
      HT[rc] = round((max(5,rpois(1,(55*PreySzLst$SzRatio[S]^1.5)*Hadj[p]))*Ncrct[rc])/(nco+1))
      ST[rc] = HT[rc]
      for(i in 1:nco){
        rc = rc + 1
        dvn = dvn + 1
        Bout[rc] = Bout[rc-1] 
        Divenum[rc] = dvn
        DT[rc] = DT[rc-1]
        Success[rc] = "c"
        SuccessV[rc] = 0.5   
        Prey[rc] = Prey[rc-1]
        PreyV[rc] = PreyV[rc-1]
        Size[rc] = Size[rc-1]
        Qual[rc] = Qual[rc-1]
        SzClss[rc] = SzClss[rc-1]
        Sz_cm[rc] = Sz_cm[rc-1]
        Sz_rt[rc] = Sz_rt[rc-1]
        Nitem[rc] = Nitem[rc-1]
        Ppnlost[rc] = 0
        Ncrct[rc] = Ncrct[rc-1]
        Mss_est[rc] =  Mss_est[rc-1] 
        HT[rc] = HT[rc-1] + round(runif(1,-10,10))
        ST[rc] = HT[rc]  
      }
    }else{
      Success[rc] = "y"
      SuccessV[rc] = 1        
      Ppnlost[rc] = prplstC[Size[rc]]*prplst[pmin(3,ceiling(rgamma(1,shape=1, rate=1.5)))]
      Ncrct[rc] = Nitem[rc] - (Ppnlost[rc]*Nitem[rc])
      Mss_est[rc] = (PFa[p]*(Pawsz*Sz_rt[rc])^PFb[p])*Ncrct[rc] 
      HT[rc] = round(max(5,rpois(1,(55*Sz_rt[rc]^1.5)*Hadj[p]))*Ncrct[rc])
      SecndPrey = rbinom(1,1,PrbSc*MultP[p])
      if(SecndPrey==1){
        p1 = p
        S1 = S
        rc = rc + 1
        Bout[rc] = Bout[rc-1]
        Divenum[rc] = Divenum[rc-1]
        DT[rc] = DT[rc-1]
        Success[rc] = "y"
        SuccessV[rc] = 1
        tmp = rmultinom(1,1,MultP*Pprey)
        p = which(tmp==1)
        prey = as.character(Preytypes[p,1])
        Prey[rc] = prey
        PreyV[rc] = p
        if(p==p1){
          if (S<4){
            S = S1+1
          }else{
            S = S1-1
          }
        }else{
          S = min(16,ceiling(rgamma(1,shape=Spar1[p]*Spar2[p], rate=Spar2[p])))
        }
        Size[rc] = PreySzLst$SzClass[S]
        Qual[rc] = PreySzLst$SzQual[S]   
        Sz_cm[rc] = 0.1*Pawsz*PreySzLst$SzRatio[S]
        Sz_rt[rc] = PreySzLst$SzRatio[S]
        SzClss[rc] = S
        Nitem[rc] = ceiling(rgamma(1,shape=Npar1[p]*Npar2[p], rate=Npar2[p]))
        Ppnlost[rc] = prplstC[Size[rc]]*prplst[pmin(3,ceiling(rgamma(1,shape=1, rate=1.5)))]
        Ncrct[rc] = Nitem[rc] - (Ppnlost[rc]*Nitem[rc])
        Mss_est[rc] = (PFa[p]*(Pawsz*Sz_rt[rc])^PFb[p])*Ncrct[rc]
        HT[rc] = round(max(5,rpois(1,(55*Sz_rt[rc]^1.5)*Hadj[p]))*Ncrct[rc])
        ST[rc] = HT[rc] + HT[rc-1]
        ST[rc-1] = HT[rc] + HT[rc-1]
      }else{
        ST[rc] = HT[rc]
      }
    }
  }
}
#
# Make "observed" data (simulate missed records) ----------------------------------------
Nrec = length(Bout)
Preynum = rep(1,Nrec)
PreyObs = Prey
PreyVObs = PreyV
SizeObs = Size
QualObs = Qual
Sz_cmObs = Sz_cm
Sz_rtObs = Sz_rt
Mss_estObs = Mss_est
NitemObs = Nitem
PpnlostObs = Ppnlost
NcrctObs = Ncrct
HTObs = HT
for (i in 1:(Nrec)){
  if(i>1){
    if(Bout[i-1]==Bout[i] & Divenum[i-1]==Divenum[i]){
      Preynum[i] = Preynum[i-1]+1
    }else if(is.na(Prey[i])){
      Preynum[i] = NA
    }
  }
  if(i<(Nrec-3)){
    if(SuccessV[i]==1 & SuccessV[i+1]!=0.5){
      prbID = inv.logit(B0 + B1*log(Sz_cm[i]) + B2*log(HT[i]) + B3*log(Sz_cm[i])*log(HT[i]))
      probID_Prey = rbinom(1,1,prbID)
      if(probID_Prey==0){
        PreyObs[i] = "uni"
        PreyVObs[i] = 0
        Mss_estObs[i] = NA
      }
    }
    if(SuccessV[i]==1 & SuccessV[i+1]!=0.5 & PreyObs[i] == "uni"){
      probID_Size = rbinom(1,1,.85)
      if(probID_Size==0){
        SizeObs[i] = NA
        QualObs[i] = NA
        Sz_cmObs[i] = NA
        Sz_rtObs[i] = NA
        Mss_estObs[i] = NA
      }
    }else if(SuccessV[i]==1 & SuccessV[i+1]!=0.5 & PreyObs[i] != "uni"){
      probID_Size = rbinom(1,1,.98)
      if(probID_Size==0){
        SizeObs[i] = NA
        QualObs[i] = NA
        Sz_cmObs[i] = NA
        Sz_rtObs[i] = NA
        Mss_estObs[i] = NA
      }      
    }
    if(SuccessV[i]==1 & SuccessV[i+1]!=0.5){
      probID_N = rbinom(1,1,inv.logit(6 - 3*Nitem[i]^0.25))
      if(probID_N==0){
        NitemObs[i] = NA
        Mss_estObs[i] = NA
        PpnlostObs[i] = NA
        NcrctObs[i] = NA
      }
    }
    if(SuccessV[i]==1 & SuccessV[i+1]!=0.5 & Preynum[i]==1){
      probID_HT = rbinom(1,1,.65)
      if(probID_HT==0){
        HTObs[i] = NA
        if(i<Nrec & Divenum[i+1]==Divenum[i] & Preynum[i+1]==Preynum[i]+1){
          nxt = 1
          while(nxt==1){
            i = i+1
            HTObs[i] = NA
            if(i<Nrec & Divenum[i+1]==Divenum[i] & Preynum[i+1]==Preynum[i]+1){
              nxt = 1
            }else{
              nxt = 0
            }
          }
        }
      }
    }  
  }
}
#
# Create data frame of forage data ---------------------------------------------------------
Fdat = data.frame(ID = seq(1,length(Bout)), Bout = Bout, Divenum = Divenum, DT=DT, ST=ST,
                  Success=Success,Preynum=Preynum,Prey=Prey,PreyV=PreyV,
                  Size=Size,Qual=Qual,Nitem=Nitem,Ppnlost=Ppnlost,
                  HT=HT,Ncrct=Ncrct,Sz_cm=Sz_cm,Sz_rt=Sz_rt,Pawsz=rep(Pawsz,length(Bout)),
                  SuccessV=SuccessV,Mss_est=Mss_est,SzClss=SzClss)
FdatObs = data.frame(ID = seq(1,length(Bout)), Bout = Bout, Divenum = Divenum, DT=DT, ST=ST,
                  Success=Success,Preynum=Preynum,Prey=PreyObs,PreyV=PreyVObs,
                  Size=SizeObs,Qual=QualObs,Nitem=NitemObs,Ppnlost=PpnlostObs,
                  HT=HTObs,Ncrct=NcrctObs,Sz_cm=Sz_cmObs,Sz_rt=Sz_rtObs,Pawsz=rep(Pawsz,length(Bout)),
                  SuccessV=SuccessV,Mss_est=Mss_estObs)
# How to assign size ratio and size_cm normally? Would need pup paw sz by bout
tmp = paste0(FdatObs$Size,FdatObs$Qual)
szes = unique(tmp); szes = szes[-which(szes=="NANA")]
FdatObs$Sz_rt = numeric(length = length(tmp))
FdatObs$Sz_cm = numeric(length = length(tmp))
for (i in 1:length(szes)){
  ii = which(tmp==szes[i])
  FdatObs$Sz_rt[ii] = as.numeric(PreySzLst$SzRatio[PreySzLst$Szcode==szes[i]])
  FdatObs$Sz_cm[ii] = as.numeric(0.1*FdatObs$Pawsz[ii]*PreySzLst$SzRatio[PreySzLst$Szcode==szes[i]])
}
# Mas estimation for records with known size, known prey ID and positive Ncrct values
FdatObs$Mass_est = numeric(length = nrow(FdatObs))
ii = which(FdatObs$Sz_cm>0 & FdatObs$PreyV>0 & FdatObs$Ncrct>0)
FdatObs$Mass_est[ii] = (PFa[FdatObs$PreyV[ii]]*(10*FdatObs$Sz_cm[ii])^PFb[FdatObs$PreyV[ii]])*FdatObs$Ncrct[ii]
FdatObs$Mass_est[FdatObs$Mass_est==0] = NA
FdatObs$Sz_rt[FdatObs$Sz_rt==0] = NA
FdatObs$Sz_cm[FdatObs$Sz_cm==0] = NA
# Un-id prey: gets index for 1 greater than number of ID'd prey types
FdatObs$PreyV[FdatObs$PreyV==0] = NPtypes+1

NpreyDv = aggregate(Fdat$Preynum, by = list(Fdat$Bout,Fdat$Divenum), max)
colnames(NpreyDv) = c("Bout","Divenum","NPreyTp")
Fdat2 = merge(Fdat, NpreyDv, by=c("Bout","Divenum"))
Fdat = Fdat2[order(Fdat2$ID),]
Fdat$Nszunits = Fdat$Size*Fdat$Ncrct
rm(Fdat2,NpreyDv)

NpreyDv = aggregate(FdatObs$Preynum, by = list(FdatObs$Bout,FdatObs$Divenum), max)
colnames(NpreyDv) = c("Bout","Divenum","NPreyTp")
Fdat2 = merge(FdatObs, NpreyDv, by=c("Bout","Divenum"))
FdatObs = Fdat2[order(Fdat2$ID),]
FdatObs$Nszunits = FdatObs$Size*FdatObs$Ncrct
rm(Fdat2,NpreyDv)

# Calc TRUE params and Allocation of time by prey, per bout -------------------
#  plus log(mean) of Sz_cm, Nitm, HT and lambda
source("Fprocess.r")
NPtypes = length(unique(Fdat$PreyV[!is.na(Fdat$PreyV)]))
rslt = Fprocess(Fdat,Nbouts,NPtypes); 
TotMin = rslt$TotMin
TotMinP = rslt$TotMinP
AlloctP = rslt$AlloctP
SmnP = rslt$SmnP
NmnP = rslt$NmnP
HTmnP = rslt$HTmnP
LMDmnP = rslt$LMDmnP
CRmnP = rslt$CRmnP
SmnP_n = rslt$SmnP_n
NmnP_n = rslt$NmnP_n
HTmnP_n = rslt$HTmnP_n
LMDmnP_n = rslt$LMDmnP_n
CRmnP_n = rslt$CRmnP_n

AlloctP_TRUE = colMeans(AlloctP)/sum(colMeans(AlloctP))
SZmn = numeric(); Sp = numeric()
NImn = numeric(); Np = numeric()
HTmn = numeric(); Hp = numeric()
LMlg = numeric(); Lp = numeric()
CRate = numeric(); Cp = numeric(); Csz = numeric(); Css = numeric()
Smn_TRUE = numeric()
Nmn_TRUE = numeric()
Hmn_TRUE = numeric()
Lmn_TRUE = numeric(); 
Cmn_TRUE = numeric();
MnN = 5
for (p in 1:NPtypes){
  ii = which(SmnP_n[,p] >= MnN)
  Smn_TRUE[p] = weighted.mean(SmnP[ii,p],SmnP_n[ii,p])
  SZmn = c(SZmn,(SmnP[ii,p])); Sp = c(Sp,rep(p,length(ii)))
  ii = which(NmnP_n[,p] >= MnN)
  Nmn_TRUE[p] = weighted.mean(NmnP[ii,p],NmnP_n[ii,p])
  NImn = c(NImn,(NmnP[ii,p])); Np = c(Np,rep(p,length(ii)))
  ii = which(HTmnP_n[,p] >= MnN)
  Hmn_TRUE[p] = weighted.mean(HTmnP[ii,p],HTmnP_n[ii,p])
  HTmn = c(HTmn,(HTmnP[ii,p])); Hp = c(Hp,rep(p,length(ii)))
  ii = which(LMDmnP_n[,p] >= MnN)
  Lmn_TRUE[p] = weighted.mean(LMDmnP[ii,p],LMDmnP_n[ii,p])
  LMlg = c(LMlg,logit(pmin(0.99,LMDmnP[ii,p]))); Lp = c(Lp,rep(p,length(ii)))   
  ii = which(CRmnP_n[,p] >= MnN)
  ft = lm(CRmnP[ii,p] ~ SmnP[ii,p])
  Cmn_TRUE[p] = coef(ft)[1] + coef(ft)[2]* Smn_TRUE[p]
  CRate = c(CRate,CRmnP[ii,p]); Cp = c(Cp,rep(p,length(ii)))      
  Csz = c(Csz,log(SmnP[ii,p])); Css = c(Css,CRmnP_n[ii,p])
}
# 
# Calc OBSERVED params and Allocation of time by prey, per bout -------------------
#  plus log(mean) of Sz_cm, HT, Nitm and lambda
UseObs = 1
if(UseObs==1){
  NPtypes = length(unique(FdatObs$PreyV[!is.na(FdatObs$PreyV)]))
  rslt = Fprocess(FdatObs,Nbouts,NPtypes)
  TotMin = rslt$TotMin
  TotMinP = rslt$TotMinP
  AlloctP = rslt$AlloctP
  SmnP = rslt$SmnP
  NmnP = rslt$NmnP
  HTmnP = rslt$HTmnP
  LMDmnP = rslt$LMDmnP
  CRmnP = rslt$CRmnP
  SmnP_n = rslt$SmnP_n
  NmnP_n = rslt$NmnP_n
  HTmnP_n = rslt$HTmnP_n
  LMDmnP_n = rslt$LMDmnP_n
  CRmnP_n = rslt$CRmnP_n
  SZmn = numeric(); Sp = numeric()
  NImn = numeric(); Np = numeric()
  HTmn = numeric(); Hp = numeric(); Hsz = numeric()
  LMlg = numeric(); Lp = numeric()
  CRate = numeric(); Cp = numeric(); Csz = numeric(); Css = numeric()
  eta_prior = numeric()
  MnN = 3
  MaxSS = 100
  Nprcaps = length(which(FdatObs$PreyV>0 & FdatObs$PreyV<NPtypes))
  for (p in 1:(NPtypes-1)){
    ii = which(FdatObs$PreyV==p)
    eta_prior[p] = (length(ii)/Nprcaps)*(NPtypes-1)*2.5
    ii = which(SmnP_n[,p] >= MnN)
    ii = ii[order(SmnP_n[ii,p], decreasing = T)][1:min(MaxSS,length(ii))]
    SZmn = c(SZmn,(SmnP[ii,p])); Sp = c(Sp,rep(p,length(ii)))
    ii = which(NmnP_n[,p] >= MnN)
    ii = ii[order(NmnP_n[ii,p], decreasing = T)][1:min(MaxSS,length(ii))]
    NImn = c(NImn,(NmnP[ii,p])); Np = c(Np,rep(p,length(ii)))
    ii = which(HTmnP_n[,p] >= MnN)
    ii = ii[order(HTmnP_n[ii,p], decreasing = T)][1:min(MaxSS,length(ii))]
    HTmn = c(HTmn,(HTmnP[ii,p])); Hp = c(Hp,rep(p,length(ii)))
    Hsz = c(Hsz,log(SmnP[ii,p]))
    ii = which(LMDmnP_n[,p] >= MnN)
    ii = ii[order(LMDmnP_n[ii,p], decreasing = T)][1:min(MaxSS,length(ii))]
    LMlg = c(LMlg,logit(pmin(0.99,LMDmnP[ii,p]))); Lp = c(Lp,rep(p,length(ii)))  
    ii = which(CRmnP_n[,p] >= MnN)
    ii = ii[order(CRmnP_n[ii,p], decreasing = T)][1:min(MaxSS,length(ii))]
    CRate = c(CRate,CRmnP[ii,p]); Cp = c(Cp,rep(p,length(ii)))      
    Csz = c(Csz,log(SmnP[ii,p])) ; Css = c(Css,CRmnP_n[ii,p])
  }
  # Repeat for un-id prey
  MaxSS = 100
  p = NPtypes
  ii = which(SmnP_n[,p] >= MnN)
  ii = ii[order(SmnP_n[ii,p], decreasing = T)][1:min(MaxSS,length(ii))]
  SZmnU = (SmnP[ii,p]); 
  ii = which(NmnP_n[,p] >= MnN)
  ii = ii[order(NmnP_n[ii,p], decreasing = T)][1:min(MaxSS,length(ii))]
  NImnU = (NmnP[ii,p]); 
  ii = which(HTmnP_n[,p] >= MnN)
  ii = ii[order(HTmnP_n[ii,p], decreasing = T)][1:min(MaxSS,length(ii))]
  HTmnU = (HTmnP[ii,p]); 
  Hsz_u = log(SmnP[ii,p])
  ii = which(LMDmnP_n[,p] >= MnN)
  ii = ii[order(LMDmnP_n[ii,p], decreasing = T)][1:min(MaxSS,length(ii))]
  LMlgU = logit(pmin(0.99,LMDmnP[ii,p])); 
}
NU = numeric()
NU[1] = length(SZmnU);  NU[2] = length(HTmnU); # NU[3] = length(NImnU);NU[4] = length(LMlgU);
NSz = length(Sp)
NNi = length(Np)
NHt = length(Hp)
NLm = length(Lp)
NCR = length(Cp)  

# Set up and run Stan model to fit to data --------------------------------------------
# Set options for STAN
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
#
fitmodel = c("SOFAfit.stan")
#
stan.data <- list(Nbouts=Nbouts,K=NPtypes,Km1=NPtypes-1,EffortP=TotMinP,
                  NSz=NSz,NHt=NHt,NCR=NCR,NU=NU,SZmnU=SZmnU,HTmnU=HTmnU,
                  Sp=Sp,Hp=Hp,Cp=Cp,SZmn=SZmn,HTmn=HTmn,Hsz=Hsz,Hsz_u=Hsz_u,
                  CRate=CRate,Csz=Csz,Cal_dns_mn=Cal_dns_mn,Cal_dns_sg=Cal_dns_sg)
#
params <- c("eta","SZ","HT","CR","PD","ER","CRmn","ERmn","Pid","tau",
            "phi1","phi2","psi1","psi2","muSZ","sigSZ","sigHT","sigCR",
            "muSZ_u","psi1_u","psi2_u","sigSZ_u","sigHT_u","maxPunid") # NI, LM,
#
nsamples <- 500
nburnin <- 250
cores = detectCores()
ncore = min(20,cores-1)
#
out <- stan(
  file = fitmodel,         # Stan program
  data = stan.data,        # named list of data
  pars = params,           # list of params to monitor
  init= "random",          # initial values    "random"            
  chains = ncore,          # number of Markov chains
  warmup = nburnin,        # number of warmup iterations per chain
  iter = nburnin+nsamples, # total number of iterations per chain
  cores = ncore,           # number of available cores 
  refresh = 100           # show progress every 'refresh' iterations
)
#
mcmc <- as.matrix(out)
vn = colnames(mcmc)
Nsims = nrow(mcmc)
sumstats = summary(out)$summary
vns = row.names(sumstats)
#
rstan::traceplot(out, pars=c("eta","Pid"), inc_warmup = F, nrow = 7)
rstan::traceplot(out, pars=c("maxPunid"), inc_warmup = F, nrow = 1)
plot(out,pars="eta")
plot(out,pars="CR")

