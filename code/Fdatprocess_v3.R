Fdatprocess <- function(dat,dfPr,dfPtp,dfPcl,dfSz,dfPaw,dfM,dfE,dfElst,Adj_Sz_grp,dfSzAd){
# Data pre-processing prior to SOFA forage analysis:
require(dplyr)
require(stats)
require(MASS)
require(gtools)
Reg = dat$Region[1]
# account for prey codes that are "non-prey" such as rock or shell...
ii = which(is.na(dfPr$PreyType))
if(length(ii)>0){
  for(i in 1:length(ii)){
    iii = ii[i]
    ix = which(dat$Prey==dfPr$PreyCode[iii])
    if(length(ix)>0){
      for(j in 1:length(ix)){
        iz = ix[j]
        dnii = dat$Divenum[iz]
        btii = dat$Bout[iz]
        sbii = dat$Subbout[iz]
        ndii = which(dat$Bout==btii & dat$Subbout==sbii & dat$Divenum==dnii)
        if(length(ndii)==1){
          dat$Success[iz] = "n"
          dat$Prey[iz] = NA
          dat$N_items[iz] = NA
          dat$Size[iz] = NA
          dat$Qualifier[iz] = NA
          dat$HT[iz] = NA
          dat$Prop_lost[iz] = NA
        }else{
          dat = dat[-iz,]
        }
      }
    }
  }
}
# Remove dives with unknown or other Success codes, create sequential dive numbers 
dat$Success[dat$Success=="co"] = "c"
dat$Success[dat$Success=="m"] = "y"
ix = which(dat$Success!="y" & dat$Success!="n" & dat$Success!="c")
dat = dat[-ix,]
dat$SuccessV = numeric(length = nrow(dat))
dat$SuccessV[which(dat$Success=="y")] = 1
dat$SuccessV[which(dat$Success=="c")] = 0.5
# Clean-up Sex and Ageclass fields as necessary
ii = which(dat$Sex!="f" & dat$Sex!="m" & dat$Sex!="u")
dat$Sex[ii] = "u"
ii = which(is.na(dat$Sex))
dat$Sex[ii] = "u"
dat$Ageclass[dat$Ageclass=="s"] = "sa"
dat$Ageclass[dat$Ageclass=="o"] = "aa"
dat$Ageclass[dat$Ageclass=="p"] = "j"
ii = which(is.na(dat$Ageclass))
dat$Ageclass[ii] = "u"
ii = which(dat$Ageclass!="j" & dat$Ageclass!="sa" & dat$Ageclass!="a" & dat$Ageclass!="aa" & dat$Ageclass!="u")
dat$Ageclass[ii] = "u"

# ** for size codes of 9 or other "illegal" values, set to NA
dat$Size[dat$Size>4] = NA
# **

# Remove bouts with no observed successful dives (or <25% successful) 
NoDiveSucc = 1
while(NoDiveSucc>0){
  ix = numeric()
  Bouts = unique(dat$Bout)
  Nbouts = length(Bouts)
  for (b in 1:Nbouts){
    ii = which(dat$Bout==Bouts[b])
    Nsucc = length(which(dat$SuccessV[ii]>0))
    Nusucc = length(which(dat$SuccessV[ii]==0))
    if((Nsucc/(Nsucc+Nusucc))<0.25){
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
dvlst$BoutN = paste0(dvlst$Date,"_",dvlst$Bout)
dvlst$BoutN = as.numeric(as.factor(dvlst$BoutN))
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
Boutlist = dat %>% group_by(BoutN,Bout) %>% 
            summarize(Date = min(Date),
                      Area = first(Area),
                      Site = first(Site),
                      Period = first(Period),
                      Sex = first(Sex),
                      Ageclass = first(Ageclass),
                      Ottername = first(Ottername),
                      Pup = first(Pup),
                      Ndv = max(DiveN),
                      NdvSucc = sum(SuccessV/Npr))
# Check for and account for duplicated bouts in Boutlist table
# caused by errors in the various attribute fields (Sex, Ageclass, etc.)
# id = unique(c(which(duplicated(Boutlist$BoutN,fromLast = TRUE)==T),
#             which(duplicated(Boutlist$BoutN)==T)))
# if(length(id)>0){
#   for(i in 1:length(id)){
#     ii = id[i]
#     chk = 0
#     while(chk==0){
#       ip = ii-1
#       ic = which(id == ip)
#       if(length(ic)==0){
#         chk=1
#       }
#     }
#     Boutlist$Ndv[ip] =  Boutlist$Ndv[ip] + Boutlist$Ndv[ii]
#     Boutlist$NdvSucc[ip] =  Boutlist$NdvSucc[ip] + Boutlist$NdvSucc[ii]
#   }
#   Boutlist = Boutlist[-id,]
# }
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
# Adjust prey size to account for minimum reasonable size
if(ncol(dfPtp)==5){
  for(i in 1:(NPtypes-1)){
    if(!is.na(dfPtp[i,5]) & is.numeric(as.numeric(dfPtp[i,5])) & dfPtp[i,5]>0 ){
      mnszmm = as.numeric(dfPtp[i,5])
      ii = which(dat$PreyV==i & is.numeric(dat$Sz_mm) & dat$Sz_mm>0 )
      dat$Sz_mm[ii] = pmax(dat$Sz_mm[ii],mnszmm)
    }
  }
}
# NOTE: allow for user over-ride adjustment of size for particular prey and group
if(Adj_Sz_grp==1){
  n_sz_ch = nrow(dfSzAd)
  for(i in 1:n_sz_ch){
    iicol = which(colnames(dat)==dfSzAd$Groupvar[i])
    ii = which(dat[,iicol]==as.character(dfSzAd$Value[i]) & dat$PreyT==dfSzAd$PreyType[i])
    if(length(ii)>0){
      dat$Sz_mm[ii] = dfSzAd$Adjust_fact[i] * dat$Sz_mm[ii]
    }
  }
}
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
iic = which(dat$SuccessV==0.5 & dat$Divenum==1)
dat$SuccessV[iic] = 1
dat$Tmtag[iic]=1
iic = which(dat$SuccessV==0.5)
# If any carry-over dives, make sure prey type is assigned to correct (original) prey
if(length(iic)>0){
  for (i in 1:length(iic)){
    ii = iic[i]
    iicc = which(dat$BoutN == dat$BoutN[ii] & dat$DiveN==dat$DiveN[ii])
    successchk = 0
    gobak = 0
    while(successchk==0){
      gobak = gobak + 1
      iipr = which(dat$BoutN == dat$BoutN[ii] & dat$DiveN==(dat$DiveN[ii]-gobak))
      # If beginning of carry over sequence is success = N, fix it
      prvsuc = dat$SuccessV[iipr[1]]
      if(is.na(prvsuc)){
        if(is.na(dat$Prey[ii])){
          dat$SuccessV[iicc] = rep(0,length(iicc))
          successchk = -1
        }else{
          dat$SuccessV[iicc] = rep(1,length(iicc))
          successchk = -1
        }
      }else if(prvsuc==0){
        if(is.na(dat$Prey[ii])){
          dat$SuccessV[iicc] = rep(0,length(iicc))
          successchk = -1
        }else{
          dat$SuccessV[iicc] = rep(1,length(iicc))
          successchk = -1
        }
      }else if(prvsuc==1){
        successchk = 1
      }
    }
    if(length(iicc)==1 & length(iipr)==1 & successchk >0){
      if(is.na(dat$PreyV[ii]) | dat$PreyV[ii]==NPtypes){
        if(is.na(dat$PreyV[iipr[1]]) | dat$PreyV[iipr[1]]==NPtypes){
          dat$PreyV[ii] = NPtypes
        }else{
          dat$PreyV[ii] = dat$PreyV[iipr[1]]
        }
      }
      if(is.na(dat$PreyV[iipr[1]]) | dat$PreyV[iipr[1]]==NPtypes){
        dat$PreyV[iipr[1]] = dat$PreyV[ii]
      }
      if( dat$PreyV[ii] != dat$PreyV[iipr[1]] ){
        dat$SuccessV[ii] = 1
      }
    }
  }
}
# For each prey code, fit mass-lng fxns, and get mean and SE of kcal/g
#  (use delta method to get combined SE for multiple proxy species)
#  then generate for each prey cap the expected mass, expected Edns, SE for each,
#  and then finally get mean SE for log mass predictions and SE_kcal by prey type
#
prx = which(startsWith(colnames(dfPr),"Prox"))
apar = numeric()
bpar = numeric()
Edns = numeric()
SE_lgmss = numeric()
SE_Edens = numeric()
MassLngFits = list()
for(i in 1:NPcodes){
  if(!is.na(dfPr$PreyType[i])){
    if(dfPr$PreyType[i] != "UNID"){
      Nprox = length(which(!is.na(dfPr[i,prx])))
      jj = numeric()
      wts = numeric()
      E_dns = numeric()
      E_dns_sd = numeric()
      E_dns_wts = numeric()
      # Vr_E_dns = numeric()
      for(j in 1:Nprox){ # 
        kk = which(dfM$SppCode== as.character(dfPr[i,prx[j]]))
        wts = c(wts,rep(1/length(kk),length(kk)))
        jj = c(jj,kk)
        cc = which(dfElst$SppCode == as.character(dfPr[i,prx[j]]))
        p_cat = dfElst$PreyCat1[cc]
        dd = which(dfElst$PreyCat1 == p_cat)
        for(q in 1:length(dd)){
          pryprx = dfElst$SppCode[dd[q]]
          ll = which(dfE$SppCode== pryprx)
          if(length(ll)>0){
            for(m in 1:length(ll)){
              E_dns = c(E_dns, log(max(.1,dfE$kcal_g_edblwet[ll[m]])))
              E_dns_wts = c(E_dns_wts,dfE$N[ll[m]])
              if(!is.na(dfE$kcal_g_edblwet_SD[ll[m]])){
                M = log(max(.1,dfE$kcal_g_edblwet[ll[m]]))
                V = (dfE$kcal_g_edblwet_SD[ll[m]])^2
                E_dns_sd = c(E_dns_sd, sqrt(log(1 + V/M^2)) )
              }else{
                E_dns_sd = c(E_dns_sd, NA)
              }
            }
          }
        }
      }
      E_dns_sd = E_dns_sd / sqrt(E_dns_wts)
      E_dns_wts = E_dns_wts / sum(E_dns_wts)
      wts = wts/max(wts)
      # NOTE: 
      x = log(dfM$MaxLinearDim_mm[jj]); y = log(dfM$EdblWetMass[jj])
      ft = rlm(y ~ x,weights = wts, method = "M",maxit = 50)
      xx = seq(min(x),max(x),by=.05)
      # plot(x,y,main=dfPr$PreyCode[i])
      # lines(xx,predict(ft,newdata = data.frame(x=xx)),col="red")
      prdct = data.frame(x=xx,ypred=predict(ft,newdata = data.frame(x=xx)))
      ft$prdct = prdct
      MassLngFits[[i]] = ft
      apar[i] = as.numeric(coef(ft)[1])
      bpar[i] = as.numeric(coef(ft)[2])
      SE_lgmss[i] = mean(predict(ft,se.fit = T)$se.fit) 
      if(length(E_dns)==1){
        Edns[i] = E_dns
        if(!is.na(E_dns_sd)){
          SE_Edens[i] = E_dns_sd
        }else{
          SE_Edens[i] = 0.2
        }
      }else if(length(E_dns)>1 & length(E_dns)<5){
        Edns[i] = weighted.mean(E_dns,E_dns_wts)
        if(!is.na(mean(E_dns_sd, na.rm=T))){
          SE_Edens[i] = max(0.05,sqrt(mean(E_dns_sd^2, na.rm=T)))
        }else{
          SE_Edens[i] = 0.2
        }
      }else{
        Edns[i] = weighted.mean(E_dns,E_dns_wts)
        if(!is.na(mean(E_dns,na.rm=T))){
          SE_Edens[i] = max(sd(Edns),sqrt(mean(E_dns_sd^2, na.rm=T)))
        }else{
          SE_Edens[i] = sd(Edns)
        }
      }
    }
  }
}
dat$Mass_est = numeric(length = Nobs)*NA
dat$Edns_est = numeric(length = Nobs)*NA
dat$SE_lgmss = numeric(length = Nobs)*NA
dat$SE_Edens = numeric(length = Nobs)*NA
for(i in 1:NPcodes){
  if(!is.na(dfPr$PreyType[i])){
    if(dfPr$PreyType[i] != "UNID"){
      jj = which(dat$Prey==dfPr$PreyCode[i] & !is.na(dat$Sz_mm) )
      if(length(jj)>0){
      	# NOTE: Changed to using edibible wet biomass, so do not mult by Ppn Edible 
        dat$Mass_est[jj] = exp( apar[i] + bpar[i] * log(dat$Sz_mm[jj]) ) # * Ped[i]
        dat$Edns_est[jj] = Edns[i]
        dat$SE_lgmss[jj] = SE_lgmss[i]
        dat$SE_Edens[jj] = SE_Edens[i]
      }
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
  jj = which(dat$PreyT==dfPtp$PreyType[i] & !is.na(dat$Edns_est) & !is.na(dat$SE_Edens))
  if(length(jj)==0){
    Cal_dns_mn[i] = log(0.75)
    Cal_dns_sg[i] = 0.25   
  }else{
    Cal_dns_mn[i] = mean(dat$Edns_est[jj],na.rm = T)
    Cal_dns_sg[i] = min(0.25,sqrt(mean(dat$SE_Edens[jj]^2,na.rm = T)))
  }
  jj = which(dat$PreyT==dfPtp$PreyType[i] & !is.na(dat$SE_lgmss) )
  if(length(jj)==0){
    logMass_sg[i] = 0.05 
  }else{
    logMass_sg[i] = sqrt(mean((dat$SE_lgmss[jj])^2,na.rm = T))
  }
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

