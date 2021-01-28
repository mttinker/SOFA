Boutprocess <- function(Fdat,Nbouts,NPtypes,Boutlist,MnN1,MnN2,GrpOpt,Ngrp){
  # Function to process foraging data
  require(dplyr)
  TotMin = matrix(0,nrow = Nbouts,ncol = 1)
  TotMinP = matrix(0,nrow = Nbouts,ncol = NPtypes)
  AlloctP = matrix(0,nrow = Nbouts,ncol = NPtypes)
  CrctF_Un = matrix(0,nrow = Nbouts,ncol = NPtypes)
  CRmnP = matrix(NA,nrow = Nbouts,ncol = NPtypes)
  SmnP = matrix(NA,nrow = Nbouts,ncol = NPtypes)
  NmnP = matrix(NA,nrow = Nbouts,ncol = NPtypes)
  HTmnP = matrix(NA,nrow = Nbouts,ncol = NPtypes)
  LMDmnP = matrix(NA,nrow = Nbouts,ncol = NPtypes)
  
  CRmnP_n = matrix(NA,nrow = Nbouts,ncol = NPtypes)
  SmnP_n = matrix(NA,nrow = Nbouts,ncol = NPtypes)
  NmnP_n = matrix(NA,nrow = Nbouts,ncol = NPtypes)
  HTmnP_n = matrix(NA,nrow = Nbouts,ncol = NPtypes)
  LMDmnP_n = matrix(NA,nrow = Nbouts,ncol = NPtypes)
  
  # Note: to allocate time, need to account for unsucc, carry over and mult items per dive
  for (b in 1:Nbouts){
    ii = which(Fdat$BoutN==b) 
    FD = Fdat[ii,]
    FD$HTT = numeric(length = nrow(FD))
    Ndv = max(FD$DiveN)
    STpr = matrix(0,nrow = Ndv,ncol = NPtypes)
    DTpr = matrix(0,nrow = Ndv,ncol = NPtypes)
    STun = matrix(0,nrow = Ndv,ncol = 1)
    DTun = matrix(0,nrow = Ndv,ncol = 1)
    NDvUn = 0 
    NDvUnP = matrix(0,nrow = Ndv,ncol = NPtypes)
    for (i in 1:Ndv){
      ii = which(FD$DiveN==i)
      if(FD$SuccessV[ii][1]==1){
        if(length(ii)>1){
          HTmlt = FD$HT[ii]
          NSU = FD$Nszunits[ii]
          if(length(which(!is.na(HTmlt)))==length(ii)){
            NSU = NA
            for(j in 1:length(ii)){
              FD$HT[ii[j]] = FD$ST[ii[j]]*(HTmlt[j]/sum(HTmlt))
              FD$HTT[ii[j]] = FD$HT[ii[j]] + FD$DT[ii[j]]/length(ii)
              STpr[i,FD$PreyV[ii[j]]] = STpr[i,FD$PreyV[ii[j]]]+FD$HT[ii[j]] 
              DTpr[i,FD$PreyV[ii[j]]] = DTpr[i,FD$PreyV[ii[j]]]+FD$DT[ii[j]]/length(ii)
            }
          }else if(length(which(!is.na(NSU)))==length(ii)){
            for(j in 1:length(ii)){
              FD$HT[ii[j]] = NA
              FD$HTT[ii[j]] = NA
              STpr[i,FD$PreyV[ii[j]]] = STpr[i,FD$PreyV[ii[j]]]+FD$ST[ii[j]]*(NSU[j]/sum(NSU))
              DTpr[i,FD$PreyV[ii[j]]] = DTpr[i,FD$PreyV[ii[j]]]+FD$DT[ii[j]]/length(ii)
            }
          }else{
            for(j in 1:length(ii)){
              FD$HT[ii[j]] = NA
              FD$HTT[ii[j]] = NA
              STpr[i,FD$PreyV[ii[j]]] = STpr[i,FD$PreyV[ii[j]]]+FD$ST[ii[j]]/length(ii)
              DTpr[i,FD$PreyV[ii[j]]] = DTpr[i,FD$PreyV[ii[j]]]+FD$DT[ii[j]]/length(ii)
            }
          }
        }else{
          FD$HT[ii[1]] = FD$ST[ii[1]]
          FD$HTT[ii[1]] = FD$HT[ii[1]] + FD$DT[ii[1]]
          STpr[i,FD$PreyV[ii[1]]] = STpr[i,FD$PreyV[ii[1]]] + FD$ST[ii[1]]
          DTpr[i,FD$PreyV[ii[1]]] = DTpr[i,FD$PreyV[ii[1]]] + FD$DT[ii[1]]
        }
      }else if(FD$SuccessV[ii][1]==0){
        # If reasonable, apply unsuc dives to specific prey
        diveThis = FD$DiveN[ii[1]]
        if (diveThis == Ndv){
          STun[i] = FD$ST[ii[1]]
          DTun[i] = FD$DT[ii[1]]
          NDvUn = NDvUn + 1
          succNext = 1
        }else{
          succNext = 0
          gofor = 0          
        }
        while(succNext<1){
          gofor = gofor + 1
          diveNext = FD$DiveN[ii[1]+gofor]
          succNext = FD$SuccessV[ii[1]+gofor]
          if (diveNext == Ndv & succNext < 1){
            STun[i] = FD$ST[ii[1]]
            DTun[i] = FD$DT[ii[1]]
            NDvUn = NDvUn + 1
            succNext = 1
          }else if(diveNext < Ndv & succNext == 0){
            succNext = 0
          }else if(diveNext < Ndv & succNext == 0.5){
            succNext = 1
            Prtyp = FD$PreyV[ii[1]+gofor]
            NDvUnP[i, Prtyp] = 1
            STpr[i, Prtyp] = FD$ST[ii[1]]
            DTpr[i, Prtyp] = FD$DT[ii[1]]          
          }else if(succNext == 1 & (!is.na(FD$Sz_mm[ii[1]+gofor]) & FD$Sz_mm[ii[1]+gofor] > 90) & FD$Npr[ii[1]+gofor]==1){          
            succNext = 1
            Prtyp = FD$PreyV[ii[1]+gofor]
            NDvUnP[i, Prtyp] = 1
            STpr[i, Prtyp] = FD$ST[ii[1]]
            DTpr[i, Prtyp] = FD$DT[ii[1]]
          }else{
            STun[i] = FD$ST[ii[1]]
            DTun[i] = FD$DT[ii[1]]
            NDvUn = NDvUn + 1
            succNext = 1
          }
        }
      }else if(FD$SuccessV[ii][1]==0.5){
        if (FD$DiveN[ii][1]==1){
          STun[i] = FD$ST[ii[1]]
          DTun[i] = FD$DT[ii[1]]
          succPrev = 1
          FD$HT[ii[1]] = NA
          FD$HTT[ii[1]] = NA
        }else{
          succPrev = 0
          gobak = 0
        }
        while(succPrev<1){
          gobak = gobak + 1
          succPrev = FD$SuccessV[ii[1]-gobak]
          if (succPrev == 1){
            FD$HT[ii[1]-gobak] = FD$HT[ii[1]-gobak] + FD$ST[ii[1]]
            FD$HTT[ii[1]-gobak] = FD$HTT[ii[1]-gobak] + FD$ST[ii[1]] + FD$DT[ii[1]]
            FD$HT[ii[1]] = NA
            FD$HTT[ii[1]] = NA
            FD$N_items[ii[1]-gobak] = ceiling(FD$N_items[ii[1]-gobak])
            FD$Ncrct[ii[1]-gobak] = ceiling(FD$Ncrct[ii[1]-gobak])
            Prtyp = FD$PreyV[ii[1]-gobak]
            STpr[FD$DiveN[ii[1]-gobak], Prtyp] = STpr[FD$DiveN[ii[1]-gobak], Prtyp] + FD$ST[ii[1]]
            DTpr[FD$DiveN[ii[1]-gobak], Prtyp] = DTpr[FD$DiveN[ii[1]-gobak], Prtyp] + FD$DT[ii[1]]
          }else if(succPrev < 1 & FD$DiveN[ii[1]-gobak]==1){
            STun[i] = FD$ST[ii[1]]
            DTun[i] = FD$DT[ii[1]]
            FD$HT[ii[1]] = NA
            FD$HTT[ii[1]] = NA
            succPrev = 1
          }
        }
      }
    }
    TotMinutes = (sum(FD$ST[which(is.na(FD$Preynum) | FD$Preynum == 1)]) +
                    sum(FD$DT[which(is.na(FD$Preynum) | FD$Preynum == 1)]))/60
    TimeAllocP = colSums(STpr) + colSums(DTpr) 
    TimeAllocUnP = colSums(STpr*NDvUnP) + colSums(DTpr*NDvUnP) 
    AlloctP[b,] = TimeAllocP/sum(TimeAllocP)
    TimeUnsP = (sum(STun) + sum(DTun))*AlloctP[b,]
    # To correct CR for unsuccesful dives and unallocated time: 
    CrctF_Un[b,] = (TimeAllocP - TimeAllocUnP) / (TimeAllocP + TimeUnsP) 
    TotMinP[b,] = round((TimeAllocP + TimeUnsP)/60)
    TotMin[b] = sum(TotMinP[b,])
    for (p in 1:NPtypes){
      # Mean Success rate by prey (and sample size)
      ndvSp = length(which(STpr[,p]>0 & NDvUnP[,p]==0))
      ndvUp = length(which(NDvUnP[,p]>0))
      LMDmnP_n[b,p] = ndvSp
      if(ndvSp>0){
        LMDmnP[b,p] = ndvSp/(ndvSp+ndvUp+NDvUn*AlloctP[b,p]) 
      }
      # Mean Sz_mm by prey (and sample size)
      ii = which(FD$SuccessV==1 & FD$PreyV==p & FD$Sz_mm>0)
      SmnP_n[b,p] = length(ii)
      if(length(ii)>0){
        SmnP[b,p] = mean(FD$Sz_mm[ii])
      }
      # Mean N_items by prey (and sample size)
      ii = which(FD$SuccessV==1 & FD$PreyV==p & FD$Ncrct>0)
      NmnP_n[b,p] = length(ii)
      if(length(ii)>0){
        NmnP[b,p] = mean(FD$Ncrct[ii])
      }
      # Mean HT by prey (and sample size) **need to also have known size
      ii = which(FD$SuccessV==1 & FD$PreyV==p & FD$Ncrct>0 & FD$HT>0 & FD$Sz_mm>0)
      HTmnP_n[b,p] = length(ii)
      if(length(ii)>0){
        HTmnP[b,p] = mean(FD$HT[ii]/FD$Ncrct[ii])
      }
      # Mean CR by prey (and sample size) ** NOTE: also select on Eflag==0
      ii = which(FD$SuccessV==1 & FD$PreyV==p & FD$Mass_est>0 & FD$HTT>0 & FD$Tmtag ==0)
      CRmnP_n[b,p] = length(ii)
      if(length(ii)>0){
        CRmnP[b,p] = mean(FD$Mass_est[ii]/(FD$HTT[ii]/60))*CrctF_Un[b,p]
        if(is.nan(CRmnP[b,p]) | is.na(CRmnP[b,p])){
          CRmnP_n[b,p] = 0
          CRmnP[b,p] = 0
        }
      }
    }
  }
  # # Remove any records with NaN for time allocation
  GrpBt =  Boutlist$Grp; GrpE = GrpBt
  ix = which(is.nan(TotMin) | is.na(TotMin) | Boutlist$NdvSucc< MnN2 )
  TotMin = TotMin[-ix]; TotMinP = TotMinP[-ix,]; AlloctP = AlloctP[-ix,]; GrpE = GrpE[-ix]
  NboutsE = length(TotMin)
  # Result = list(TotMin=TotMin,TotMinP=TotMinP,AlloctP=AlloctP,
  #               CRmnP=CRmnP,SmnP=SmnP,NmnP=NmnP,HTmnP=HTmnP,LMDmnP=LMDmnP,
  #               CRmnP_n=CRmnP_n,SmnP_n=SmnP_n,NmnP_n=NmnP_n,HTmnP_n=HTmnP_n,
  #               LMDmnP_n=LMDmnP_n)
  # return(Result)    
  SZmn = numeric(); Sp = numeric(); Sss = numeric(); Sg = numeric()
  NImn = numeric(); Np = numeric(); Nss = numeric(); Ng = numeric()
  HTmn = numeric(); Hp = numeric(); Hsz = numeric(); Hss = numeric(); Hg = numeric()
  LMlg = numeric(); Lp = numeric(); Lss = numeric(); Lg = numeric()
  CRate = numeric(); Cp = numeric(); Csz = numeric(); Css = numeric(); Cg = numeric()
  eta_prior = numeric()
  MaxSS = 100
  Nprcaps = length(which(Fdat$PreyV>0 & Fdat$PreyV<NPtypes))
  # Cap Crate at 250
  for (p in 1:(NPtypes-1)){
    #
    ii = which(Fdat$PreyV==p)
    eta_prior[p] = (length(ii)/Nprcaps)*(NPtypes-1)*2.5
    #
    ii = which(SmnP_n[,p] >= MnN1)
    if(length(ii)<6){ii = which(SmnP_n[,p] > 0)}
    if(length(ii)>0){
      ii = ii[order(SmnP_n[ii,p], decreasing = T)][1:min(MaxSS,length(ii))]
      SZmn = c(SZmn,(SmnP[ii,p])); Sp = c(Sp,rep(p,length(ii)));  Sss = c(Sss,SmnP_n[ii,p])
      Sg = c(Sg,GrpBt[ii])
    }
    #
    ii = which(HTmnP_n[,p] >= MnN1 & SmnP_n[,p] > 0)
    if(length(ii)<6){ii = which(HTmnP_n[,p] > 0 & SmnP_n[,p] > 0)}
    if(length(ii)>0){
      ii = ii[order(HTmnP_n[ii,p], decreasing = T)][1:min(MaxSS,length(ii))]
      HTmn = c(HTmn,(HTmnP[ii,p])); Hp = c(Hp,rep(p,length(ii)))
      Hsz = c(Hsz,2.5*log(SmnP[ii,p])-7); Hss = c(Hss,HTmnP_n[ii,p])
      Hg = c(Hg,GrpBt[ii])
    }
    #
    ii = which(NmnP_n[,p] >= MnN1)
    if(length(ii)<6){ii = which(NmnP_n[,p] > 0)}
    if(length(ii)>0){
      ii = ii[order(NmnP_n[ii,p], decreasing = T)][1:min(MaxSS,length(ii))]
      NImn = c(NImn,(NmnP[ii,p])); Np = c(Np,rep(p,length(ii))); Nss = c(Nss,NmnP_n[ii,p])
      Ng = c(Ng,GrpBt[ii]) 
    }
    #
    ii = which(LMDmnP_n[,p] >= MnN1)
    if(length(ii)==0){ii = which(LMDmnP_n[,p] > 0)}
    if(length(ii)>0){    
      ii = ii[order(LMDmnP_n[ii,p], decreasing = T)][1:min(MaxSS,length(ii))]
      LMlg = c(LMlg,logit(pmin(0.99,LMDmnP[ii,p]))); Lp = c(Lp,rep(p,length(ii))) ; Lss = c(Lss,LMDmnP_n[ii,p]) 
      Lg = c(Lg,GrpBt[ii]) 
    }
    #
    ii = which(CRmnP_n[,p] >= MnN1)
    if(length(ii)<6){ii = which(CRmnP_n[,p] > 0)}
    if(length(ii)>0){    
      ii = ii[order(CRmnP_n[ii,p], decreasing = T)][1:min(MaxSS,length(ii))]
      CRate = c(CRate,pmin(250,CRmnP[ii,p])); Cp = c(Cp,rep(p,length(ii)))      
      Csz = c(Csz,2.5*log(SmnP[ii,p])-7) ; Css = c(Css,CRmnP_n[ii,p])
      Cg = c(Cg,GrpBt[ii]) 
    }
  }
  # Repeat for un-id prey
  p = NPtypes
  ii = which(SmnP_n[,p] >= MnN1)
  if(length(ii)<6){ii = which(SmnP_n[,p] > 0)}
  if(length(ii)>0){ 
  ii = ii[order(SmnP_n[ii,p], decreasing = T)][1:min(5*MaxSS,length(ii))]
  SZmnU = (SmnP[ii,p]); Sss_u = SmnP_n[ii,p]; 
  Sg_u = GrpBt[ii]
  }
  ii = which(NmnP_n[,p] >= MnN1 )
  if(length(ii)<6){ii = which(NmnP_n[,p] > 0)}
  if(length(ii)>0){   
  ii = ii[order(NmnP_n[ii,p], decreasing = T)][1:min(5*MaxSS,length(ii))]
  NImnU = (NmnP[ii,p]); Nss_u = NmnP_n[ii,p]
  Ng_u = GrpBt[ii]
  }
  ii = which(HTmnP_n[,p] >= MnN1 & SmnP_n[,p] > 0)
  if(length(ii)<6){ii = which(HTmnP_n[,p] > 0 & SmnP_n[,p] > 0)}
  if(length(ii)>0){     
  ii = ii[order(HTmnP_n[ii,p], decreasing = T)][1:min(5*MaxSS,length(ii))]
  HTmnU = (HTmnP[ii,p]); Hss_u = HTmnP_n[ii,p]
  Hsz_u = 2.5*log(SmnP[ii,p])-7
  Hg_u = GrpBt[ii]
  }
  ii = which(LMDmnP_n[,p] >= MnN1)
  if(length(ii)<6){ii = which(LMDmnP_n[,p] > 0)}
  if(length(ii)>0){   
  ii = ii[order(LMDmnP_n[ii,p], decreasing = T)][1:min(MaxSS,length(ii))]
  LMlgU = logit(pmin(0.99,LMDmnP[ii,p])); Lss_u = LMDmnP_n[ii,p]
  Lg_u = GrpBt[ii]
  }
  NU = numeric()
  NU[1] = length(SZmnU);  NU[2] = length(HTmnU); # NU[3] = length(NImnU);NU[4] = length(LMlgU);
  NSz = length(Sp)
  NNi = length(Np)
  NHt = length(Hp)
  NLm = length(Lp)
  NCR = length(Cp)  
  #
  if(GrpOpt==0){
    Result = list(Nbouts=NboutsE,K=NPtypes,Km1=NPtypes-1,EffortP=TotMinP,
                  NSz=NSz,NHt=NHt,NCR=NCR,NLm=NLm,NU=NU,SZmnU=SZmnU,HTmnU=HTmnU,
                  Sp=Sp,Hp=Hp,Cp=Cp,SZmn=SZmn,HTmn=HTmn,Hsz=Hsz,Hsz_u=Hsz_u,
                  LMlg=LMlg,Lp=Lp,Lss=Lss,LMlgU=LMlgU,Lss_u=Lss_u,
                  Sss=Sss,Sss_u=Sss_u,Hss=Hss,Hss_u=Hss_u,CRate=CRate,Css=Css,Csz=Csz,
                  Cal_dns_mn=Cal_dns_mn,Cal_dns_sg=Cal_dns_sg,logMass_sg=logMass_sg)
  }else{
    Result = list(Nbouts=NboutsE,K=NPtypes,Km1=NPtypes-1,EffortP=TotMinP,
                  NSz=NSz,NHt=NHt,NCR=NCR,NLm=NLm,NU=NU,SZmnU=SZmnU,HTmnU=HTmnU,
                  Sp=Sp,Hp=Hp,Cp=Cp,SZmn=SZmn,HTmn=HTmn,Hsz=Hsz,Hsz_u=Hsz_u,
                  LMlg=LMlg,Lp=Lp,Lss=Lss,LMlgU=LMlgU,Lss_u=Lss_u,Lg=Lg,Lg_u=Lg_u,
                  Sss=Sss,Sss_u=Sss_u,Hss=Hss,Hss_u=Hss_u,CRate=CRate,Css=Css,Csz=Csz,
                  Cal_dns_mn=Cal_dns_mn,Cal_dns_sg=Cal_dns_sg,logMass_sg=logMass_sg,
                  Ngrp=Ngrp,GrpE=GrpE,Sg=Sg,Hg=Hg,Cg=Cg,Sg_u=Sg_u,Hg_u=Hg_u)
  }
  return(Result)  
}
