Fprocess <- function(Fdat,Nbouts,NPtypes){
  # Function to process foraging data
  Boutlist = unique(Fdat$Bout)
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
    ii = which(Fdat$Bout==Boutlist[b]) 
    FD = Fdat[ii,]
    FD$HTT = numeric(length = nrow(FD))
    Ndv = max(FD$Divenum)
    STpr = matrix(0,nrow = Ndv,ncol = NPtypes)
    DTpr = matrix(0,nrow = Ndv,ncol = NPtypes)
    STun = matrix(0,nrow = Ndv,ncol = 1)
    DTun = matrix(0,nrow = Ndv,ncol = 1)
    NDvUn = 0 
    NDvUnP = matrix(0,nrow = Ndv,ncol = NPtypes)
    for (i in 1:Ndv){
      ii = which(FD$Divenum==i)
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
        diveThis = FD$Divenum[ii[1]]
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
          diveNext = FD$Divenum[ii[1]+gofor]
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
          }else if(succNext == 1 & (!is.na(FD$Size[ii[1]+gofor]) & FD$Size[ii[1]+gofor] > 1) & FD$NPreyTp[ii[1]+gofor]==1){          
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
        if (FD$Divenum[ii][1]==1){
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
            FD$Nitem[ii[1]-gobak] = ceiling(FD$Nitem[ii[1]-gobak])
            FD$Ncrct[ii[1]-gobak] = ceiling(FD$Ncrct[ii[1]-gobak])
            Prtyp = FD$PreyV[ii[1]-gobak]
            STpr[FD$Divenum[ii[1]-gobak], Prtyp] = STpr[FD$Divenum[ii[1]-gobak], Prtyp] + FD$ST[ii[1]]
            DTpr[FD$Divenum[ii[1]-gobak], Prtyp] = DTpr[FD$Divenum[ii[1]-gobak], Prtyp] + FD$DT[ii[1]]
          }else if(succPrev < 1 & FD$Divenum[ii[1]-gobak]==1){
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
      # Mean Sz_cm by prey (and sample size)
      ii = which(FD$SuccessV==1 & FD$PreyV==p & FD$Sz_cm>0)
      SmnP_n[b,p] = length(ii)
      if(length(ii)>0){
        SmnP[b,p] = mean(FD$Sz_cm[ii])
      }
      # Mean Nitem by prey (and sample size)
      ii = which(FD$SuccessV==1 & FD$PreyV==p & FD$Ncrct>0)
      NmnP_n[b,p] = length(ii)
      if(length(ii)>0){
        NmnP[b,p] = mean(FD$Ncrct[ii])
      }
      # Mean HT by prey (and sample size)
      ii = which(FD$SuccessV==1 & FD$PreyV==p & FD$Ncrct>0 & FD$HT>0)
      HTmnP_n[b,p] = length(ii)
      if(length(ii)>0){
        HTmnP[b,p] = mean(FD$HT[ii]/FD$Ncrct[ii])
      }
      # Mean CR by prey (and sample size)
      ii = which(FD$SuccessV==1 & FD$PreyV==p & FD$Mss_est>0 & FD$HTT>0)
      CRmnP_n[b,p] = length(ii)
      if(length(ii)>0){
        CRmnP[b,p] = mean(FD$Mss_est[ii]/(FD$HTT[ii]/60))*CrctF_Un[b,p]
      }
    }
  }
  Result = list(TotMin=TotMin,TotMinP=TotMinP,AlloctP=AlloctP,
                CRmnP=CRmnP,SmnP=SmnP,NmnP=NmnP,HTmnP=HTmnP,LMDmnP=LMDmnP,
                CRmnP_n=CRmnP_n,SmnP_n=SmnP_n,NmnP_n=NmnP_n,HTmnP_n=HTmnP_n,
                LMDmnP_n=LMDmnP_n)
  return(Result)    
}
