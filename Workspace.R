# Work area for exploring functions and parameter effects
require(gtools)
require(fitdistrplus)
# Function for consumption rate variation among bouts-------------------------------------- 
# Assumes greater variance if Crate based on fewer dives
Ndv = seq(1,20)
reps = 10000
Mrt = 10
p1 = .3
p2 = .7
b = p1*Ndv^(p2)
Cse = numeric()
Cmn = numeric()
tmp = matrix(0,nrow = reps,ncol = length(Ndv))
for (i in 1:length(Ndv)){
  tmp[,i] = rgamma(reps,shape=Mrt*b[i],rate = b[i])
  Cmn[i] = mean(tmp[,i])
  Cse[i] = sd(tmp[,i])
}
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(Ndv,Cse,type="l",main=paste0("Std Err of Crt vs Ndives, p1 = ",p1,", p2 = ",p2))
hist(pmin(32,tmp[,1]),
     main="dist of Crt for 1 dive",xlim = c(0,32),breaks = seq(0,32,by=2))
hist(pmin(32,tmp[,20]),
     main="dist of Crt for 20 dives",xlim = c(0,32),breaks = seq(0,32,by=2))
layout(matrix(c(1), 1, 1, byrow = TRUE))

# Function to create N items per dive--------------------------------------
# low number items (1 or 2), with some proportion lost: M = 0.5 and b = 7
# mod number items (2 to 6): M = 3 and b = 2
# larger number items (5 to 15): M = 9 and b = 1
prplst = c(0,.5,.25); Ndvs = 20
M = 3.5
b = 2
hist(ceiling(rgamma(10000,shape=M*b, rate=b)))
mn_Ni = numeric()
log_Mn_Ni = numeric()
for (i in 1:1000){
  Nitm = ceiling(rgamma(Ndvs,shape=M*b, rate=b))
  # Next function calculates proportion lost (mostly 0, some 0.5 and 0.25)
  Nlst = prplst[pmin(3,ceiling(rgamma(Ndvs,shape=1, rate=1.5)))]
  # Combine to get number items consumed
  tmp = Nitm-Nlst
  # Mean and log mean number items consumed per dive
  mn_Ni[i] = mean(tmp)
  # log_Mn_Ni[i] = mean(log(tmp))
  log_Mn_Ni[i] = log(mean(tmp))
}
hist(log_Mn_Ni)
plot(fitdist(log_Mn_Ni,"norm"))

# Evaluate WA data ---------------------------------------------------------
dat = read_excel("./data/WA_SOFA_14March2017.xls")
Pdat = read_excel("./data/Prey_data_Northern.xls")

boutWA = unique(dat$BOUT)
nbouts_WA = length(boutWA)
preytp_WA = unique(dat$PREY)
ii = which(preytp_WA=="uni" | is.na(preytp_WA))
preytp_WA = preytp_WA[-ii]
NPtypes_WA = length(preytp_WA)
WA_prey = data.frame(Prey = preytp_WA[order(preytp_WA)], PreyN = seq(1,NPtypes_WA))
dat$PreyV = numeric(length=nrow(dat)); dat$PreyV = NA
ii = which(dat$PREY=="uni")
dat$PreyV[ii] = 0
for (j in 1:NPtypes_WA){
  ii = which(dat$PREY==WA_prey$Prey[j])
  dat$PreyV[ii] = WA_prey$PreyN[j]
}

UID_n = numeric(length=nbouts_WA)
UID_min = numeric(length=nbouts_WA)
ID_n = matrix(nrow = nbouts_WA, ncol=NPtypes_WA)
ID_min = matrix(nrow = nbouts_WA, ncol=NPtypes_WA)
Tot_min = numeric(length=nbouts_WA)
ID_ppnT = matrix(nrow = nbouts_WA, ncol=NPtypes_WA)
for(i in 1:nbouts_WA){
  ii = which(dat$BOUT==boutWA[i])
  FD = dat[ii,]
  ii = which(FD$PreyV==0)
  UID_n[i] = length(ii)
  UID_min[i] = sum(FD$ST[ii])
  for (j in 1:NPtypes_WA){
    ii = which(FD$PreyV==j)
    ID_n[i,j] = length(ii)
    ID_min[i,j] = sum(FD$ST[ii])
  }
  Tot_min[i] = UID_min[i] + sum(ID_min[i,])
  ID_ppnT[i,] = ID_min[i,]/Tot_min[i]
}
UID_ppnT = UID_min/Tot_min
ii = which(is.na(Tot_min) | Tot_min == 0 | UID_ppnT == 1)
hist(UID_ppnT[-ii])
hist(ID_ppnT[-ii,1])


# DELTA METHOD, for computing variance of combined variables ------------------
# (either average or weighted average)
#
#draw huge sample from 3 variables (x1, x2 and x3) 
set.seed(1)
#random variable with variance 25
x1<-rnorm(n=1000000, mean=3, sd=5)
#random variable with variance 81
x2<-rnorm(n=1000000, mean=2, sd=9)
#random variable with variance 36
x3<-rnorm(n=1000000, mean=1, sd=6)
XX = cbind(x1,x2,x3)
# let q be proportional contributiosn of each xn to combined X
q = c(.1, .4, .5)
# compute the weighted average, based on proportional contribution to combined var:

Xwt = XX[,1]*q[1] + XX[,2]*q[2] + XX[,3]*q[3]
# empiirical mean and sd:
mean(X.wt)
sd(X.wt)
# Calculated mean and sd, delta method:
mean_approx = sum(q*c(mean(x1),mean(x2),mean(x3)))
sd_approx = sqrt(var(x1)*q[1]^2 + var(x2)*q[2]^2 + var(x3)*q[3]^2)

# Un-wieghted mean:
# compute the average
x.bar<-(x1+x2+x3)/3
mean(x.bar)
mean_approx = mean(c(mean(x1),mean(x2),mean(x3)))
  
# compute the sd of the average , should be close to sqrt((25+81+36)/(3^2))
sd(x.bar)^2
sd_approx = sqrt( (var(x1)+var(x2)+var(x3))/(3^2) )


# Bhattacharyya distance and coefficient (= degree of overlap between 2 nortmal dists)
mu = c(.99,1.42)
sg = c(.16,.14)

D12 = 0.25*log(0.25*(sg[1]^2/sg[2]^2+sg[2]^2/sg[1]^2 +2))+0.25*(((mu[1]-mu[2])^2)/(sg[1]^2+sg[2]^2))
OL = exp(-D12)
OL






