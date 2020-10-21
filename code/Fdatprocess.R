Fdatprocess <- function(dat,dfPr,dfPtp,dfPcl,dfSz,dfPaw,dfM,dfE,dfElst){
# Data pre-processing prior to SOFA forage analysis:
require(dplyr)
require(stats)
require(MASS)
require(gtools)
Reg = dat$Region[1]
# Drop dives with unknown or other Success, create sequential dive numbers 
ix = which(dat$Success!="y" & dat$Success!="n" & dat$Success!="c")
dat = dat[-ix,]
dat$SuccessV = numeric(length = nrow(dat))
dat$SuccessV[which(dat$Success=="y")] = 1
dat$SuccessV[which(dat$Success=="c")] = 0.5
# Remove bouts with no observed successful dives
NoDiveSucc = 1
while(NoDiveSucc>0){
  ix = numeric()
  Bouts = unique(dat$Bout)
  Nbouts = length(Bouts)
  for (b in 1:Nbouts){
    ii = which(dat$Bout==Bouts[b])
    Nsucc = sum(dat$SuccessV[ii])
    if(Nsucc==0){
      ix = c(ix,ii)
    }
  }
  if(length(ix)>0){
    dat = dat[-ix,]
  }else{
    NoDiveSucc = 0
  }
}
dvlst = dat %>% group_by(Date,Bout,Subbout,Divenum) %>% 
  summarize(Npr = length(Prey))
dvlst$BoutN = as.numeric(as.factor(dvlst$Bout))
Nbouts = max(dvlst$BoutN)
Ndives = nrow(dvlst)
Nobs = nrow(dat)
dvlst$DiveN = numeric(length = Ndives)
for(b in 1:Nbouts){
  ii = which(dvlst$BoutN==b)
  dvlst$DiveN[ii] = seq(1,length(ii))
}
dat = merge(dat,dvlst,by=c("Date","Bout","Subbout","Divenum"),sort=F)
#
Boutlist = dat %>% group_by(BoutN,Bout,Date,Area,Site,Period,
                            Sex,Ageclass,Ottername,Pup) %>% 
            summarize(Ndv = max(DiveN),
            NdvSucc = sum(SuccessV/Npr))
# Check for and account for duplicated bouts in Boutlist table
# caused by errors in the various attribute fields (Sex, Ageclass, etc.)
id = which(duplicated(Boutlist$BoutN)==T)
if(length(id)>0){
  for(i in 1:length(id)){
    ii = id[i]
    chk = 0
    while(chk==0){
      ip = ii-1
      ic = which(id == ip)
      if(length(ic)==0){
        chk=1
      }
    }
    Boutlist$Ndv[ip] =  Boutlist$Ndv[ip] + Boutlist$Ndv[ii]
    Boutlist$NdvSucc[ip] =  Boutlist$NdvSucc[ip] + Boutlist$NdvSucc[ii]
  }
}
Boutlist = Boutlist[-id,]
# Create Sz_cm field, accounting for prey sizeclass and Qualifier and paw Size
# 1. Sizequal field for Size class
Sizequal = paste0( as.character(dat$Size),
                   replace(dat$Qualifier, is.na(dat$Qualifier), ""))
Sizequal[is.na(dat$Size)] = NA
Pawratio = merge(data.frame(ord = seq(1,length(Sizequal)), Szqual = Sizequal),
                 dfSz[,1:2],by=c("Szqual"),all.x=T,sort = F)
Pawratio = Pawratio$Pawratio[order(Pawratio$ord)]
Paw = merge(cbind(dat$Sex,dat$Ageclass,ord = seq(1,length(Sizequal))),
            dfPaw,by.x=c(1,2),by.y=c(1,2),sort=F)
Paw = Paw[order(as.numeric(Paw$ord)),]; 
dat$Paw_wd = Paw$Pawidth; dat$Pawratio = Pawratio; rm(Paw,Pawratio)
dat$Sz_mm = dat$Paw_wd*dat$Pawratio
# Override Sz_mm for prey caps where "est_cm" was recorded 
ii = which(!is.na(dat$est_cm) & dat$est_cm>0)
if(length(ii)>0){
  dat$Sz_mm[ii] = 10*dat$est_cm[ii]
}  
#
# Assign prey types (accounting for large/small classes where necessary)
#  and associated number codes (PreyV)
PreyT = character(length = Nobs)
NPtypes = nrow(dfPtp) 
NPcodes = nrow(dfPr)
for (i in 1:NPcodes){
  if (i==1){
    if(!is.na(dfPr$SzBrkmm[i])){
      ii = which(dat$Prey==dfPr$PreyCode[i] & dat$Sz_mm <= dfPr$SzBrkmm[i])
      PreyT[ii] = dfPr$PreyType[i]
    }else{
      ii = which(dat$Prey==dfPr$PreyCode[i])
      PreyT[ii] = dfPr$PreyType[i]
    }
  }else{
    if(!is.na(dfPr$SzBrkmm[i])){
      ii = which(dat$Prey==dfPr$PreyCode[i] & dat$Sz_mm <= dfPr$SzBrkmm[i])
      PreyT[ii] = dfPr$PreyType[i]
    }else if(!is.na(dfPr$SzBrkmm[i-1])){
      ii = which(dat$Prey==dfPr$PreyCode[i-1] & dat$Sz_mm > dfPr$SzBrkmm[i-1])
      dat$Prey[ii] = paste0(dat$Prey[ii],c("2"))
      PreyT[ii] = dfPr$PreyType[i]
    }else{
      ii = which(dat$Prey==dfPr$PreyCode[i])
      PreyT[ii] = dfPr$PreyType[i]
    }
  }
}
dat$PreyT = PreyT
tmp = merge(cbind(dat$PreyT,seq(1,Nobs)),dfPtp[,c(2,1)],by=c(1,1),sort=F,all.x=T)
tmp = tmp[order(as.numeric(tmp$V2)),]
dat$PreyV = tmp$TypeN
# fill in missing DT and ST with appropriate mean values
Tmtag = rep(0,Nobs)
zz = which(is.na(dat$DT))
Tmtag[zz] = 1
for(i in 1:length(zz)){
  ii = which(dat$Success==dat$Success[zz[i]] & dat$PreyT==dat$PreyT[zz[i]])
  dat$DT[zz[i]] = round(mean(dat$DT[ii],na.rm = T))
}
zz = which(is.na(dat$ST))
Tmtag[zz] = 1
for(i in 1:length(zz)){
  ii = which(dat$Success==dat$Success[zz[i]] & dat$PreyT==dat$PreyT[zz[i]])
  dat$ST[zz[i]] = round(mean(dat$ST[ii],na.rm = T))
}
dat$Tmtag = Tmtag
# For each prey code, fit mass-lng fxns, and get mean and SE of kcal/g
#  (use delta method to get combined SE for multiple proxy species)
#  then generate for each prey cap the expected mass, expected Edns, SE for each,
#  and then finally get mean SE for log mass predictions and SE_kcal by prey type
#
prx = which(startsWith(colnames(dfPr),"Prox"))
apar = numeric()
bpar = numeric()
Ped = numeric()
Edns = numeric()
SE_lgmss = numeric()
SE_Edens = numeric()
MassLngFits = list()
for(i in 1:NPcodes){
  if(dfPr$PreyType[i] != "UNID"){
    Nprox = length(which(!is.na(dfPr[i,prx])))
    jj = numeric()
    P_ed = numeric()
    wts = numeric()
    E_dns = numeric()
    Vr_E_dns = numeric()
    for(j in 1:Nprox){ # 
      kk = which(dfM$SppCode== as.character(dfPr[i,prx[j]]))
      wts = c(wts,rep(1/length(kk),length(kk)))
      jj = c(jj,kk)
      P_ed = c(P_ed,dfElst$AvgOfPrpn_Edible[which(dfElst$SppCode==as.character(dfPr[i,prx[j]]))])
      ll = which(dfE$SppCode== as.character(dfPr[i,prx[j]]))
      E_dns = c(E_dns,mean(dfE$kcal_g_edblwet[ll]))
      if(length(ll)<3){
        Vr_E_dns = c(Vr_E_dns,(mean(dfE$kcal_g_edblwet_SD[ll]))^2)
      }else{
        Vr_E_dns = c(Vr_E_dns,(sd(dfE$kcal_g_edblwet[ll])/sqrt(length(ll)))^2)
      }
    }
    wts = wts/max(wts)
    Ped[i] = mean(P_ed)
    x = log(dfM$MaxLinearDim_mm[jj]); y = log(dfM$TotWetMass[jj])
    ft = rlm(y ~ x,weights = wts, method = "MM")
    xx = seq(min(x),max(x),by=.05)
    # plot(x,y,main=dfPr$PreyCode[i])
    # lines(xx,predict(ft,newdata = data.frame(x=xx)),col="red")
    prdct = data.frame(x=xx,ypred=predict(ft,newdata = data.frame(x=xx)))
    ft$prdct = prdct
    MassLngFits[[i]] = ft
    apar[i] = as.numeric(coef(ft)[1])
    bpar[i] = as.numeric(coef(ft)[2])
    SE_lgmss[i] = mean(predict(ft,se.fit = T)$se.fit)  
    Edns[i] = mean(E_dns)
    SE_Edens[i] = sqrt(sum(Vr_E_dns)/ (Nprox^2) )
  }
}
dat$Mass_est = numeric(length = Nobs)*NA
dat$Edns_est = numeric(length = Nobs)*NA
dat$SE_lgmss = numeric(length = Nobs)*NA
dat$SE_Edens = numeric(length = Nobs)*NA
for(i in 1:NPcodes){
  if(dfPr$PreyType[i] != "UNID"){
    jj = which(dat$Prey==dfPr$PreyCode[i] & !is.na(dat$Sz_mm) )
    if(length(jj)>0){
      dat$Mass_est[jj] = exp( apar[i] + bpar[i] * log(dat$Sz_mm[jj]) ) * Ped[i]
      dat$Edns_est[jj] = Edns[i]
      dat$SE_lgmss[jj] = SE_lgmss[i]
      dat$SE_Edens[jj] = SE_Edens[i]
    }
  }
}
# Override Mass_est for prey caps where "est_kg" was recorded 
ii = which(!is.na(dat$est_kg) & dat$est_kg>0)
if(length(ii)>0){
  dat$Mass_est[ii] = pmin(5000,1000*dat$est_kg[ii])
}  
# Get mean sigma values for SE values by prey type
Cal_dns_mn = numeric(length = NPtypes-1)
Cal_dns_sg = numeric(length = NPtypes-1)
logMass_sg = numeric(length = NPtypes-1)
for(i in 1:(NPtypes-1)){
  jj = which(dat$PreyT==dfPtp$PreyType[i])
  Cal_dns_mn[i] = mean(dat$Edns_est[jj],na.rm = T)
  Cal_dns_sg[i] = mean(dat$SE_Edens[jj],na.rm = T)
  logMass_sg[i] = mean(dat$SE_lgmss[jj],na.rm = T)
}
PPnlost = dat$Prop_lost
PPnlost[is.na(PPnlost)] = 0
dat$Ncrct = dat$N_items - PPnlost*dat$N_items
dat$Nszunits = dat$Size*dat$Ncrct

Result = list(Fdat=dat,Cal_dns_mn=Cal_dns_mn,Cal_dns_sg=Cal_dns_sg,logMass_sg=logMass_sg,
              Nbouts=Nbouts,Ndives=Ndives,Nobs=Nobs,NPtypes=NPtypes,Boutlist=Boutlist,
              MassLngFits=MassLngFits)
return(Result) 
}

