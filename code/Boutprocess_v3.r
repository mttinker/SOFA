Boutprocess <- function(Fdat,Nbouts,NPtypes,Boutlist,MnN1,MnN2,GrpOpt,Ngrp){
  # Function to process foraging data
  require(dplyr)
  require(MASS)
  require(gtools)
	# Average max prey consumption rate (g/min of feeding time) based on empirical observations
  #  *NOTE: includes dive times as well as surface times
	CR_max = 35
	# Min size (mm) for which a previous series of uncsuccesful dives associated with prey capture
	Sz_Lgry = 75
	Ratio_ST = sum(Fdat$ST,na.rm=T) / (sum(Fdat$ST,na.rm=T) + sum(Fdat$DT,na.rm=T))
  TotMin = matrix(0,nrow = Nbouts,ncol = 1)
  TotMinP = matrix(0,nrow = Nbouts,ncol = NPtypes)
  AlloctP = matrix(0,nrow = Nbouts,ncol = NPtypes)
  CrctF_Un = matrix(0,nrow = Nbouts,ncol = NPtypes)
  ECRmnP = matrix(NA,nrow = Nbouts,ncol = NPtypes)
  # ARmnP = matrix(NA,nrow = Nbouts,ncol = NPtypes)
  SmnP = matrix(NA,nrow = Nbouts,ncol = NPtypes)
  GmnP = matrix(NA,nrow = Nbouts,ncol = NPtypes)
  NmnP = matrix(NA,nrow = Nbouts,ncol = NPtypes)
  HTmnP = matrix(NA,nrow = Nbouts,ncol = NPtypes)
  LMDmnP = matrix(NA,nrow = Nbouts,ncol = NPtypes)
  ECRmnP_n = matrix(NA,nrow = Nbouts,ncol = NPtypes)
  SmnP_n = matrix(NA,nrow = Nbouts,ncol = NPtypes)
  GmnP_n = matrix(NA,nrow = Nbouts,ncol = NPtypes)
  NmnP_n = matrix(NA,nrow = Nbouts,ncol = NPtypes)
  HTmnP_n = matrix(NA,nrow = Nbouts,ncol = NPtypes)
  LMDmnP_n = matrix(NA,nrow = Nbouts,ncol = NPtypes)
  # Note: to allocate foraging effort, need to account for unsucc, carry over and mult items per dive
  # We assume unsuccessful dives are allocated to prey types in ppn to freq with which
  # each prey observed in bout, adjusted for prey size (larger prey >> unsuccesful dives)
  SzMean = numeric()
  for (p in 1:NPtypes){
   ii = which(Fdat$SuccessV==1 & Fdat$PreyV==p & Fdat$Sz_mm>0 )
   SzMean[p] = mean(Fdat$Sz_mm[ii],na.rm = T)  
  }
  for (b in 1:Nbouts){
    ii = which(Fdat$BoutN==b) 
    FD = Fdat[ii,]
    FD$HTT = numeric(length = nrow(FD))
    FD$Preynum = rep(NA,nrow(FD))
    FD$Preynum[which(FD$SuccessV>0)] = 0
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
        FD$Preynum[ii][1]=1
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
          }else if(succNext == 1 & (!is.na(FD$Sz_mm[ii[1]+gofor]) & FD$Sz_mm[ii[1]+gofor] >= Sz_Lgry) & FD$Npr[ii[1]+gofor]==1){          
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
        FD$Preynum[ii][1]=1
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
          iipr = which(FD$DiveN==(i-gobak))
          succPrev = FD$SuccessV[iipr[1]]
          if (succPrev == 1){
            for(j in 1:length(ii)){
              matchprey = 0
              if(is.na(FD$PreyV[ii[j]])){
              	iiprv = which(!is.na(FD$PreyV[iipr]))
              	if (length(iiprv)>0){
              		FD$PreyV[ii[j]] = FD$PreyV[iipr[iiprv[1]]]
              	}
              }
              k = 0
              while(k<length(iipr)){
              	k = k + 1
                if (!is.na(FD$PreyV[ii[j]]) & !is.na(FD$PreyV[iipr[k]]) & 
                		FD$PreyV[ii[j]]==FD$PreyV[iipr[k]]){
                  matchprey = 1
                  Prtyp = FD$PreyV[iipr[k]]
                  FD$HT[iipr[k]] = FD$HT[iipr[k]] + FD$ST[ii[j]]/length(ii)
                  FD$HTT[iipr[k]] = FD$HTT[iipr[k]] + FD$ST[ii[j]] + FD$DT[ii[j]]/length(ii)
                  FD$HT[ii[j]] = NA
                  FD$HTT[ii[j]] = NA
                  FD$N_items[iipr[k]] = ceiling(FD$N_items[iipr[k]])
                  FD$Ncrct[iipr[k]] = ceiling(FD$Ncrct[iipr[k]])
                  STpr[(i-gobak), Prtyp] = STpr[(i-gobak), Prtyp] + FD$ST[ii[j]]/length(ii)
                  DTpr[(i-gobak), Prtyp] = DTpr[(i-gobak), Prtyp] + FD$DT[ii[j]]/length(ii)
                  k = length(iipr)
                }
              }
              if (matchprey == 0){
              	if(!is.na(FD$PreyV[ii[j]])){
              		Prtyp = FD$PreyV[ii[j]]
              		FD$HT[ii[j]] = NA
              		FD$HTT[ii[j]] = NA
              		STpr[i, Prtyp] = STpr[i, Prtyp] + FD$ST[ii[j]]/length(ii)
              		DTpr[i, Prtyp] = DTpr[i, Prtyp] + FD$DT[ii[j]]/length(ii)
              	}else{
              		FD$HT[ii[j]] = NA
              		FD$HTT[ii[j]] = NA              		
              	}
              }
            }
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
    ppn_alloc_raw = TimeAllocP/sum(TimeAllocP)
    ppn_alloc_adj = ppn_alloc_raw * (SzMean/max(SzMean))^2
    AlloctP[b,] = ppn_alloc_adj/sum(ppn_alloc_adj)
    TimeUnsP = (sum(STun) + sum(DTun)) * AlloctP[b,]
    # Need to correct CR for allocted and unallocated unsuccesful dives: 
    CrctF_Un[b,] = (TimeAllocP - TimeAllocUnP) / (TimeAllocP + TimeUnsP) 
    ii = which(is.infinite(CrctF_Un[b,]))
    CrctF_Un[b,ii] = 1
    TotMinP[b,] = round((TimeAllocP + TimeUnsP)/60)
    TotMin[b] = sum(TotMinP[b,])
    for (p in 1:NPtypes){
      # Mean Success rate by prey (and sample size)
      ndvSp = length(which(STpr[,p]>0 & NDvUnP[,p]==0))
      ndvUp = length(which(NDvUnP[,p]>0))
      LMDmnP_n[b,p] = ndvSp
      if(ndvSp >0){
        LMDmnP[b,p] = ndvSp/(ndvSp+ndvUp+NDvUn*AlloctP[b,p]) 
      }
      # Mean Sz_mm by prey (and sample size)
      ii = which(FD$SuccessV==1 & FD$PreyV==p & FD$Sz_mm>0 )
      SmnP_n[b,p] = length(ii)
      if(length(ii)>0){
        SmnP[b,p] = mean(FD$Sz_mm[ii])
        GmnP[b,p] = mean(FD$Mass_est[ii])
      }
      # Mean edible biomass by prey 
      ii = which(FD$SuccessV==1 & FD$PreyV==p & FD$Sz_mm>0 & FD$Mass_est > 0)
      GmnP_n[b,p] = length(ii)
      if(length(ii)>0){
        GmnP[b,p] = mean(FD$Mass_est[ii])
      }      
      # Mean N_items by prey (and sample size)
      ii = which(FD$SuccessV==1 & FD$PreyV==p & FD$Ncrct>0.5)
      NmnP_n[b,p] = length(ii)
      if(length(ii)>0){
        NmnP[b,p] = mean(FD$Ncrct[ii])
      }
      # Mean HT by prey (and sample size) **need to also have known size
      ii = which(FD$SuccessV==1 & FD$PreyV==p & FD$Ncrct>0.5 & FD$HT>0 & FD$Sz_mm>0)
      HTmnP_n[b,p] = length(ii)
      if(length(ii)>0){
        HTmnP[b,p] = mean(FD$HT[ii]/FD$Ncrct[ii])
      }
      # Mean encounter rate CRmn by prey (and sample size) ** NOTE: also select on Tmtag==0
      if(p < NPtypes){
      	ii = which(FD$SuccessV==1 & FD$PreyV==p & FD$Sz_mm>0 & FD$HTT>0 & (FD$HTT[ii] - FD$HT[ii]) > 30 & FD$Ncrct > 0.9 & FD$Tmtag ==0)
      	ECRmnP_n[b,p] = length(ii)
      	if(length(ii)>0){
      		ECRmnP[b,p] = mean(FD$Ncrct[ii] / ((FD$HTT[ii] - FD$HT[ii])/60 )) * CrctF_Un[b,p]
      		if(is.nan(ECRmnP[b,p]) | is.na(ECRmnP[b,p] | is.infinite(ECRmnP[b,p]) | ECRmnP[b,p] > 50 | ECRmnP[b,p] < 0.01 )){
      			ECRmnP_n[b,p] = 0
      			ECRmnP[b,p] = NA
      			# ARmnP[b,p] = 0
      		}
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
  #               ECRmnP=ECRmnP,SmnP=SmnP,NmnP=NmnP,HTmnP=HTmnP,LMDmnP=LMDmnP,
  #               ECRmnP_n=ECRmnP_n,SmnP_n=SmnP_n,NmnP_n=NmnP_n,HTmnP_n=HTmnP_n,
  #               LMDmnP_n=LMDmnP_n)
  # return(Result)    
  SZmn = numeric(); Sp = numeric(); Sss = numeric(); Sg = numeric()
  NImn = numeric(); Np = numeric(); Nss = numeric(); Ng = numeric()
  HTmn = numeric(); Hp = numeric(); Hsz = numeric(); Hss = numeric(); Hg = numeric()
  LMlg = numeric(); Lp = numeric(); Lss = numeric(); Lg = numeric()
  ECRate = numeric(); Cp = numeric(); Csz = numeric(); Css = numeric(); Cg = numeric()
  MnMass = numeric()
  MLpar1 = numeric()
  MLpar2 = numeric()
  # Replace rho by array of correlation matrices
  rho_HS = array(1,dim=c(NPtypes,2,2))
  MaxSST = 150
  MaxSS = max(25,round(MaxSST/Ngrp))
  Nprcaps = length(which(Fdat$PreyV>0 & Fdat$PreyV<NPtypes))
  set.seed(123)
  for (p in 1:(NPtypes-1)){
    ii = which(!is.na(SmnP[,p]) & !is.na(GmnP[,p]))
    wts = SmnP_n[ii,p] 
    x = log(SmnP[ii,p]); y = log(GmnP[ii,p])
    ft = rlm(y ~ x,weights = wts, method = "M",maxit = 50)
    # summary(ft)
    # xx = seq(min(x),max(x),by=.05)
    # plot(y ~ x)
    # lines(xx,predict(ft,newdata = data.frame(x=xx)),col="red")
    # 
    MLpar1[p] = coef(ft)[1]
    MLpar2[p] = coef(ft)[2]
    # Add variance of size-mass fxn reflecting multiple species in prey type
    logMass_sg[p] = sqrt(logMass_sg[p]^2 + (predict(ft,se.fit = T)$residual.scale)^2)
    # Correlation between prey size and HT
    ii = which(!is.na(SmnP[,p]) & !is.na(HTmnP[,p]) )
    rho_HS[p,,] = c(1, rep(max(0.05,cor(log(SmnP[ii,p]),log(HTmnP[ii,p]))),2),1)
    #
    ii = which(Fdat$PreyV==p & !is.na(Fdat$Mass_est) )
    if(length(ii) == 0){
      MnMass[p] = 10
    }else{
      MnMass[p] = mean(Fdat$Mass_est[ii],na.rm=T)
    }
    #
    ii = which(SmnP_n[,p] >= MnN1)
    if(length(ii)<6){ii = which(SmnP_n[,p] > 0)}
    if(length(ii)>0){
      if(length(ii)>MaxSST){
        iii = numeric()
        for (g in 1:Ngrp){
        	if(length(which(GrpBt[ii]==g)) < 25 ){
        		iiis = ii[GrpBt[ii]==g]
        	}else{
        		iiis = sample(ii[GrpBt[ii]==g],min(MaxSS,length(which(GrpBt[ii]==g))),replace = FALSE)
        	}
        	iii = c(iii,iiis)
        }
        ii = iii
      }
      SZmn = c(SZmn,(SmnP[ii,p])); Sp = c(Sp,rep(p,length(ii)));  Sss = c(Sss,SmnP_n[ii,p])
      Sg = c(Sg,GrpBt[ii])
    }
    #
    ii = which(HTmnP_n[,p] >= MnN1 & HTmnP[,p] >0 & SmnP_n[,p] > 0)
    if(length(ii)<6){ii = which(HTmnP_n[,p] > 0 & HTmnP[,p] >0 & SmnP_n[,p] > 0)}
    if(length(ii)>0){
      if(length(ii)>MaxSST){
        iii = numeric()
        for (g in 1:Ngrp){
        	if(length(which(GrpBt[ii]==g)) < 25 ){
        		iiis = ii[GrpBt[ii]==g]
        	}else{
        		iiis = sample(ii[GrpBt[ii]==g],min(MaxSS,length(which(GrpBt[ii]==g))),replace = FALSE)
        	}
        	iii = c(iii,iiis)
        }
        ii = iii
      }
      HTmn = c(HTmn,(HTmnP[ii,p])); Hp = c(Hp,rep(p,length(ii)))
      Hsz = c(Hsz,2.5*log(SmnP[ii,p])-7); Hss = c(Hss,HTmnP_n[ii,p])
      Hg = c(Hg,GrpBt[ii])
    }
    #
    ii = which(NmnP_n[,p] >= MnN1)
    if(length(ii)<6){ii = which(NmnP_n[,p] > 0)}
    if(length(ii)>0){
      if(length(ii)>MaxSST){
        iii = numeric()
        for (g in 1:Ngrp){
        	if(length(which(GrpBt[ii]==g)) < 25 ){
        		iiis = ii[GrpBt[ii]==g]
        	}else{
        		iiis = sample(ii[GrpBt[ii]==g],min(MaxSS,length(which(GrpBt[ii]==g))),replace = FALSE)
        	}
        	iii = c(iii,iiis)
        }
        ii = iii
      }
      NImn = c(NImn,(NmnP[ii,p])); Np = c(Np,rep(p,length(ii))); Nss = c(Nss,NmnP_n[ii,p])
      Ng = c(Ng,GrpBt[ii]) 
    }
    #
    ii = which(!is.na(LMDmnP[,p]) & LMDmnP_n[,p] >= 5)
    # if(length(ii)==0){ii = which(LMDmnP_n[,p] > 0)}
    if(length(ii)>0){  
      if(length(ii)>MaxSST){
        iii = numeric()
        for (g in 1:Ngrp){
        	if(length(which(GrpBt[ii]==g)) < 25 ){
        		iiis = ii[GrpBt[ii]==g]
        	}else{
        		iiis = sample(ii[GrpBt[ii]==g],min(MaxSS,length(which(GrpBt[ii]==g))),replace = FALSE)
        	}
        	iii = c(iii,iiis)
        }
        ii = iii
      }
      LMlg = c(LMlg,logit(pmin(0.99, pmax(.2,LMDmnP[ii,p]) ))); 
      Lp = c(Lp,rep(p,length(ii))) ; Lss = c(Lss,LMDmnP_n[ii,p]) 
      Lg = c(Lg,GrpBt[ii]) 
    }
    #
    ii = which(ECRmnP_n[,p] >= MnN1)
    if(length(ii)<6){ii = which(ECRmnP_n[,p] > 0)}
    if(length(ii)>0){    
      if(length(ii)>MaxSST){
        iii = numeric()
        for (g in 1:Ngrp){
        	if(length(which(GrpBt[ii]==g)) < 25 ){
        		iiis = ii[GrpBt[ii]==g]
        	}else{
        		iiis = sample(ii[GrpBt[ii]==g],min(MaxSS,length(which(GrpBt[ii]==g))),replace = FALSE)
        	}
        	iii = c(iii,iiis)
        }
        ii = iii
      }
      ECRate = c(ECRate,ECRmnP[ii,p]); Cp = c(Cp,rep(p,length(ii)))      
      Csz = c(Csz,2.5*log(SmnP[ii,p])-7) ; Css = c(Css,ECRmnP_n[ii,p])
      Cg = c(Cg,GrpBt[ii]) 
    }else{ # If no reliable encounter rates for this prey type...
    	ECRate = c(ECRate,NA); Cp = c(Cp,rep(p,1))      
    	Csz = c(Csz,2.5*log(mean(SmnP[,p],na.rm=T))-7) ; Css = c(Css,1)
    	Cg = c(Cg,1) 
    }
  }
  #
  # *** NOTE: for prey with ni reliable capture rates, use lower quartile for all other prey  
  ii = which(is.infinite(ECRate) | ECRate < 0.01)
  ECRate = ECRate[-ii]
  Cp = Cp[-ii]
  Csz = Csz[-ii]
  Css = Css[-ii]
  Cg = Cg[-ii]
  ii = which(is.na(ECRate))
  ECRate[ii] = quantile(ECRate, prob = 0.25, na.rm = T)
  # **** 
  #
  # Repeat for un-id prey
  p = NPtypes
  ii = which(!is.na(SmnP[,p]) & !is.na(HTmnP[,p]) )
  rho_HS[p,,] = c(1, rep(max(0.05,cor(log(SmnP[ii,p]),log(HTmnP[ii,p]))),2),1)
  ii = which(SmnP_n[,p] >= MnN1)
  if(length(ii)<6){ii = which(SmnP_n[,p] > 0)}
  if(length(ii)>0){ 
    if(length(ii)>(2*MaxSST)){
      iii = numeric()
      for (g in 1:Ngrp){
      	if(length(which(GrpBt[ii]==g)) < 50 ){
      		iiis = ii[GrpBt[ii]==g]
      	}else{
      		iiis = sample(ii[GrpBt[ii]==g],min(MaxSS*2,length(which(GrpBt[ii]==g))),replace = FALSE)
      	}
      	iii = c(iii,iiis)
      }
      ii = iii
    }
    SZmnU = (SmnP[ii,p]); Sss_u = SmnP_n[ii,p]; 
    Sg_u = GrpBt[ii]
  }
  ii = which(NmnP_n[,p] >= MnN1 )
  if(length(ii)<6){ii = which(NmnP_n[,p] > 0)}
  if(length(ii)>0){   
    if(length(ii)>(2*MaxSST)){
      iii = numeric()
      for (g in 1:Ngrp){
      	if(length(which(GrpBt[ii]==g)) < 50 ){
      		iiis = ii[GrpBt[ii]==g]
      	}else{
      		iiis = sample(ii[GrpBt[ii]==g],min(MaxSS*2,length(which(GrpBt[ii]==g))),replace = FALSE)
      	}
      	iii = c(iii,iiis)
      }
      ii = iii
    }
    NImnU = (NmnP[ii,p]); Nss_u = NmnP_n[ii,p]
    Ng_u = GrpBt[ii]
  }
  ii = which(HTmnP_n[,p] >= MnN1 & HTmnP[,p] >0 & SmnP_n[,p] > 0)
  if(length(ii)<6){ii = which(HTmnP_n[,p] > 0 & HTmnP[,p] >0 & SmnP_n[,p] > 0)}
  if(length(ii)>0){     
    if(length(ii)>(2*MaxSST)){
      iii = numeric()
      for (g in 1:Ngrp){
      	if(length(which(GrpBt[ii]==g)) < 50 ){
      		iiis = ii[GrpBt[ii]==g]
      	}else{
      		iiis = sample(ii[GrpBt[ii]==g],min(MaxSS*2,length(which(GrpBt[ii]==g))),replace = FALSE)
      	}
      	iii = c(iii,iiis)
      }
      ii = iii
    }
    HTmnU = (HTmnP[ii,p]); Hss_u = HTmnP_n[ii,p]
    Hsz_u = 2.5*log(SmnP[ii,p])-7
    Hg_u = GrpBt[ii]
  }
  ii = which(LMDmnP_n[,p] >= MnN1)
  if(length(ii)<6){ii = which(LMDmnP_n[,p] > 0)}
  if(length(ii)>0){   
    if(length(ii)>(2*MaxSST)){
      iii = numeric()
      for (g in 1:Ngrp){
      	if(length(which(GrpBt[ii]==g)) < 50 ){
      		iiis = ii[GrpBt[ii]==g]
      	}else{
      		iiis = sample(ii[GrpBt[ii]==g],min(MaxSS*2,length(which(GrpBt[ii]==g))),replace = FALSE)
      	}
      	iii = c(iii,iiis)
      }
      ii = iii
    }
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
  EffortP = TotMinP
  # eta_pri = colSums(EffortP[,-NPtypes]) + 10; eta_pri = 20*(eta_pri/sum(eta_pri))+1
  eta_pri = rep(1,NPtypes-1)
  #
  if(GrpOpt==0){
    Result = list(Nbouts=NboutsE,K=NPtypes,Km1=NPtypes-1,EffortP=EffortP,
                  NSz=NSz,NHt=NHt,NCR=NCR,NLm=NLm,NU=NU,SZmnU=SZmnU,HTmnU=HTmnU,
                  Sp=Sp,Hp=Hp,Cp=Cp,SZmn=SZmn,HTmn=HTmn,Hsz=Hsz,Hsz_u=Hsz_u,
                  LMlg=LMlg,Lp=Lp,Lss=Lss,CR_max=CR_max,eta_pri=eta_pri,MnMass=MnMass,
                  Sss=Sss,Sss_u=Sss_u,Hss=Hss,Hss_u=Hss_u,ECRate=ECRate,Css=Css,Csz=Csz,
                  Cal_dns_mn=Cal_dns_mn,Cal_dns_sg=Cal_dns_sg,logMass_sg=logMass_sg,
                  MLpar1=MLpar1,MLpar2=MLpar2,rho_HS=rho_HS)
  }else{
    Result = list(Nbouts=NboutsE,K=NPtypes,Km1=NPtypes-1,EffortP=EffortP,
                  NSz=NSz,NHt=NHt,NCR=NCR,NLm=NLm,NU=NU,SZmnU=SZmnU,HTmnU=HTmnU,
                  Sp=Sp,Hp=Hp,Cp=Cp,SZmn=SZmn,HTmn=HTmn,Hsz=Hsz,Hsz_u=Hsz_u,
                  LMlg=LMlg,Lp=Lp,Lss=Lss,Lg=Lg,CR_max=CR_max,eta_pri=eta_pri,MnMass=MnMass,
                  Sss=Sss,Sss_u=Sss_u,Hss=Hss,Hss_u=Hss_u,ECRate=ECRate,Css=Css,Csz=Csz,
                  Cal_dns_mn=Cal_dns_mn,Cal_dns_sg=Cal_dns_sg,logMass_sg=logMass_sg,
                  Ngrp=Ngrp,GrpE=GrpE,Sg=Sg,Hg=Hg,Cg=Cg,Sg_u=Sg_u,Hg_u=Hg_u,
                  MLpar1=MLpar1,MLpar2=MLpar2,rho_HS=rho_HS)
  }
  return(Result)  
}
