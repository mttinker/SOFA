# SOFA2 Shell (Sea Otter Forage Analysis = SOFA, version 2)
# This is a shell file for running an analysis of sea otter foraging data, 
# will be replaced with Interactive Rmarkdown front end (with embedded SHiny widgets)
#
# ***NOTES***
#
#  - Save "Efficiency" as result stat?  Not sure what this would be... 
#  - Save ratio of "Unassigned" time to "Assigned" time-allocation for prey
#
# ** TO PREP FOR MODEL RUN (done by running "SOFA_prep" in main SOFA folder)
#  - allows user to select a data file to analyze, has written instructions 
#  - ***NOTE*** if data sub-folder already exists, prompt user whether to over-write files!!!
#  - creates project folder as sub-folder of "data" folder with user-specified name,
#   adds various files including:
#     - a re-named version of the user-selected forage data file called "Forage_data.xlsx",
#     (with all text fields lowercase, etc.), 
#     - the prey key spreadsheet (ready to be edited), 
#     - "reference" spreadsheets (paw size, etc.)
#     - also creates temporary prey mass and energy table spreadsheets in the data sub-folder
#     that the user can edit as desired (e.g. to use only focal region data for a species)
#    - generates rmarkdown of size-biomass data/fxns
#  *** NOTE: ALERT USER IF PROJECT DIRECTORY ALREADY EXISTS!
#    (and if user chooses to proceed, move pre-existing files to a "oldfiles" folder)
# - Run script (SOFA_run) starts by getting user to select existing project folder (created by SOFA_prep)
#     and checking to make sure that Prey_key is filled out (prey types assigned) 
#    - makes "results" sub-folder in the project folder for results rdata files (and other results obj.) 
#    - allows user to select options such as grouping levels, sub-setting, #burnin/iter, etc.
#    - initiates model fit (tracks progress?) and saves results with unique name
# - NOTE: can have multiple results rdata files for different model fits 
#     (e.g. if different group options)
# - "SOFA_sum" script to review results and output summary spreadsheets: 
#    user selects a project and results file, and then it renders an rmarkown as html
#
#
# USER INPUT-------------------------------------------------------------------
Datafile <- rstudioapi::selectFile(caption = "Select CSV File",
                               filter = "CSV Files (*.csv)",
                               existing = TRUE)

Projectname = rstudioapi::showPrompt("Project name", "Enter project name", default = "edit_me")

# Name of data file (NOTE" replace with a file selection user interface widget)
#Datafile = "../WA_SOFA_25January2018_C_NP_SP_small.xlsx"
# Projectname = "WASH_2010_2020_SA_sm"
#
MnN1 = 3 # Minimum number dives to consider for estimating prey attributes
MnN2 = 2 # Minimum number dives to consider for estimating effort allocation
# Desired number posterior samples and burn-in reps for Bayesian analysis:
nsamples <- 10000 # Increase to 10000 for best results
nburnin <- 500    # Increase to 1000 for best results
#
# END USER INPUT-----------------------------------------------------------
#
# Prep data ---------------------------------------------------------------
# Note 
require(readxl)
require(openxlsx)
dat0 = read_excel(Datafile)
colnames(dat0) = tolower(colnames(dat0))
dat = data.frame(Region=toupper(dat0$region),
                 Area=toupper(dat0$area),
                 Site=(dat0$site),
                 Period=tolower(dat0$period),
                 Date=dat0$date,
                 Sex=tolower(dat0$sex),
                 Ageclass=tolower(dat0$ageclass),
                 Pup=tolower(dat0$pup),
                 Ottername=tolower(dat0$ottername),
                 Bout=dat0$bout,Subbout=dat0$subbout,
                 Divenum=as.numeric(dat0$divenum),
                 DT = as.numeric(dat0$dt), ST = as.numeric(dat0$st),
                 Success = tolower(substr(dat0$success, 1, 1)),
                 Prey = tolower(dat0$prey),
                 N_items = as.numeric(dat0$n_items),
                 Size=as.numeric(dat0$size),
                 Qualifier=tolower(dat0$qualifier),
                 HT = as.numeric(dat0$ht),
                 Prop_lost = as.numeric(dat0$prop_lost),
                 How_lost = tolower(dat0$how_lost),
                 est_kg = as.numeric(dat0$est_kg),
                 est_cm = as.numeric(dat0$est_cm))
# Save "dat" to project directory
dir.create(paste0("./data/",Projectname),showWarnings = F)
write.xlsx(dat,file=paste0("./data/",Projectname,"/Forage_dat.xlsx"))
rm(dat0)
#
# *** Add function to create Prey_Key spreadsheet from template (so user can edit it)
#    and place that, plus all other reference spreadsheets, in project folder
#    (assumed that user then edits Prey_key and any other files)
#
# Create smaller version of data for testing (after creating Boutlist)
# ii = which(Boutlist$NdvSucc>12)
# iib = numeric()
# for (i in 1:length(ii)){
#   iib = c(iib, which(dat$Bout==Boutlist$Bout[ii[i]]))
# }
# dat = dat[iib,]
# Projectname = paste0(Projectname,"_sm")
# dir.create(paste0("./data/",Projectname),showWarnings = F)
# write.xlsx(dat,file=paste0("./data/",Projectname,"/Forage_dat.xlsx"))
#
# Process data---------------------------------------------------------------
# Load data (after selecting Projectname from available sub-folders of data folder)
#   dat = read_excel(paste0("./data/",Projectname,"/Forage_dat.xlsx"))
# Load reference spreadsheets
dfPr = read_excel(paste0("./data/",Projectname,"/Prey_Key.xlsx"))
dfPtp = read_excel(paste0("./data/",Projectname,"/Prey_Key.xlsx"),sheet = "Prey_Types")
dfPcl = read_excel(paste0("./data/",Projectname,"/Prey_Key.xlsx"),sheet = "Prey_Classes")
dfSz = read_excel(paste0("./data/",Projectname,"/Sizeclass_key.xlsx"))
dfPaw = read_excel(paste0("./data/",Projectname,"/Pawsz.xlsx"))
dfM = read_excel(paste0("./data/",Projectname,"/Mass_lng_dat.xlsx"))
dfE = read_excel(paste0("./data/",Projectname,"/Energy_dat.xlsx"))
dfElst = read_excel(paste0("./data/",Projectname,"/Prey_Spcs.xlsx"))
#
# Check Prey_key to make sure filled in properly before allowing next steps
#
source("Fdatprocess.R")
# *** Add "sub-set" option here, to filter dat if desired...
#  data filtered to one level of a categorical variable of interest (e.g. Area ==CORE),
#  so stan.data is created with data for just the subset
#  NOTE: should do this selection before Fdatprocess, to make steps easier,
#  *** Savename should include the level value ("Area_CORE" for example)
# so "Result_CORE..." instead of "Result_ALL..."
#  DO NOT ALLOW GROUPING IF SBST==1
#
rslt <- Fdatprocess(dat,dfPr,dfPtp,dfPcl,dfSz,dfPaw,dfM,dfE,dfElst)
Fdat = rslt$Fdat; Boutlist=rslt$Boutlist
Cal_dns_mn = rslt$Cal_dns_mn; Cal_dns_sg = rslt$Cal_dns_sg
logMass_sg = rslt$logMass_sg
Nbouts = rslt$Nbouts; Ndives = rslt$Ndives; Nobs = rslt$Nobs; NPtypes = rslt$NPtypes
MassLngFits = rslt$MassLngFits
rm(rslt)
#
# To plot an individual mass-length model fit, choose prey code ID (pr):
# pr = 8; ft = MassLngFits[[pr]]
# plot(exp(ft$model$x),exp(ft$model$y),main=dfPr$Description[pr],xlab="size (mm)",ylab="mass (g)")
# lines(exp(ft$prdct$x),exp(ft$prdct$ypred),col="red")
#
# Allow user to select 1-3 grouping variables (or select none for no groups), 
#   allowable group vars: Area, Site, Period, Sex, Ageclass, Ottername, Pup 
#   NOTE: filter list first such that it only presents vars with >1 levels
#  and warn user to select "bigger" containers first
Grpvar = c("Area")
Ngrpvar = length(Grpvar)
if(Ngrpvar==0){GrpOpt=0}else{GrpOpt=1}
# Determine group ID for each bout, if GrpOpt = 1
if(GrpOpt==0){
  Boutlist$Grp = rep(1,Nbouts)
  Ngrp = 1
}else{
  Grp = Boutlist %>%
    group_by_at(Grpvar) %>%
    summarize(Nbouts = length(Bout))
  Grouplist = as.data.frame(cbind(GroupID = seq(1,nrow(Grp)),Grp))
  if(Ngrpvar==1){
    Groupname = Grouplist[ , Grpvar ]
  }else{
    Groupname = apply( Grouplist[ , Grpvar ] , 1 , paste , collapse = "-" )
  }
  Grouplist = cbind(Grouplist,Groupname)
  tmp = merge(Boutlist,Grouplist[,-ncol(Grouplist)],by=Grpvar)
  tmp = tmp[order(tmp$BoutN),]; Boutlist$Grp = tmp$GroupID; rm(tmp)
  Ngrp = max(Boutlist$Grp)
}  
#
# Analyze bouts to get bout-specific stats
source("Boutprocess.R")
stan.data <- Boutprocess(Fdat,Nbouts,NPtypes,Boutlist,MnN1,MnN2,GrpOpt,Ngrp)
#
# Fit Bayesian model ------------------------------------------------------
require(parallel)
require(rstan)
require(ggplot2)
require(bayesplot)
#
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
cores = detectCores()
ncore = min(20,cores-1)
Niter = round(nsamples/ncore)+nburnin
#
if(GrpOpt==0){
  params <- c("tauB","maxPunid","CRmn","ERmn","muSZ","muSZ_u",
              "SZ","SZ_u","HT","HT_u","CR","ER","eta","PD","Pid",
              "phi1","phi2","psi1","psi2","psi1_u","psi2_u",
              "sigCR","sigSZ","sigSZ_u","sigHT","sigHT_u") #
}else{
  params <- c("tauB","tauG","maxPunid","CRmn","ERmn","muSZ","muSZ_u",
              "SZ","SZ_u","HT","HT_u","CR","ER","eta","PD","Pid",
              "phi1","phi2","psi1","psi2","psi1_u","psi2_u",
              "sigCR","sigSZ","sigSZ_u","sigHT","sigHT_u",
              "CRgmn","ERgmn","SZg","SZg_u","HTg","CRg","ERg","PDg",
              "PidG","muSZG","muSZG_u","phi1G","psi1G","psi1G_u","etaG",
              "sg1","sg2","sg3","sg4","sg5") #
}
#
# Stan model to fit:
if (GrpOpt==0){
  fitmodel = "SOFAfit.stan"
}else{
  fitmodel = "SOFAfit_Grp.stan"
}
# Compile model if necessary, or load pre-compiled version if available
testfile = paste0(substr(fitmodel ,1,nchar(fitmodel)-5),"_dso.rdata")
if(file_test("-f", testfile) ==TRUE){
  load(testfile)
}else{
  Stan_model.dso = stan_model(file = fitmodel, model_name = substr(fitmodel ,1,nchar(fitmodel)-5))
  save(Stan_model.dso,file=testfile)
}
#
#  If in markdown chunck, use results="hide" in chunk header
out <- sampling(Stan_model.dso,
                data = stan.data,    # named list of data
                pars = params,       # list of params to monitor
                init = "random",     # initial values     "random", 
                chains = ncore,         # number of Markov chains
                warmup = nburnin,       # number of warmup iterations per chain
                iter = Niter,         # total number of iterations per chain
                cores = ncore,          # number of cores (if <20, increase iter)
                refresh = 50        # show progress every 'refresh' iterations
#               control = list(adapt_delta = 0.99, max_treedepth = 15)
)
#
mcmc <- as.matrix(out)
vn = colnames(mcmc)
Nsims = nrow(mcmc)
sumstats = summary(out)$summary
vns = row.names(sumstats)
# save image file of results for post-fit processing and review:
dir.create(paste0("./data/",Projectname,"/results"),showWarnings = F)
#
if (GrpOpt==0){
  save.image(paste0("./data/",Projectname,"/results/Rslt_All_", format(Sys.time(), "%Y_%b_%d_%H"),"hr.rdata"))
}else{
  save.image(paste0("./data/",Projectname,"/results/Rslt_Grp_", format(Sys.time(), "%Y_%b_%d_%H"),"hr.rdata"))
}
#
# Summary plots ----------------------------------
#
r_mcmc = sample(Nsims,min(2500,Nsims),replace = F)
post = mcmc[r_mcmc,]
#traceplot(out,"tauG")
traceplot(out,"tauB[3]")
#
plot_title = "Posterior distributions, mean Consumption Rate (g/min) and Energy Intake Rate (kcal/min)"
suppressMessages(print(
mcmc_areas(post,
  pars = c("CRmn", "ERmn"),
  prob = 0.8) + ggtitle(plot_title,subtitle = paste0("Forage data ", Projectname))  
  # scale_x_continuous(name="Parameter value",limits=c(0,30))
))
plot_title = "Proportional contribution to diet (by biomass)"
suppressMessages(print(
mcmc_areas(post,area_method = c("equal height"),
  pars = paste0("PD[",seq(1,(NPtypes-1)),"]"),
  prob = 0.8) + 
  scale_y_discrete(name="Diet Type",labels=dfPtp$Description) +
  ggtitle(plot_title,subtitle = paste0("Forage data ", Projectname)) 
))
plot_title = "Proportional contribution to diet (by biomass)"
suppressMessages(print(
mcmc_intervals(post,
           pars = paste0("PD[",seq(1,(NPtypes-1)),"]")) +
  scale_y_discrete(name="Diet Type",labels=dfPtp$Description) +
  ggtitle(plot_title,subtitle = paste0("Forage data ", Projectname)) 
))
plot_title = "Proportional allocation of foraging effort"
suppressMessages(print(
mcmc_intervals(post,
               pars = paste0("eta[",seq(1,(NPtypes-1)),"]")) +
  scale_y_discrete(name="Diet Type",labels=dfPtp$Description) +
  ggtitle(plot_title,subtitle = paste0("Forage data ", Projectname)) 
))
plot_title = "Rate of Energy Intake (kcal/min) by prey type"
suppressMessages(print(
mcmc_intervals(mcmc,
               pars = paste0("CR[",seq(1,(NPtypes-1)),"]")) +
  scale_y_discrete(name="Diet Type",labels=dfPtp$Description) +
  ggtitle(plot_title,subtitle = paste0("Forage data ", Projectname)) 
))
plot_title = "Probability of correct ID (not recorded as 'Un-ID')"
suppressMessages(print(
mcmc_intervals(mcmc,
               pars = paste0("Pid[",seq(1,(NPtypes-1)),"]")) +
  scale_y_discrete(name="Diet Type",labels=dfPtp$Description) +
  ggtitle(plot_title,subtitle = paste0("Forage data ", Projectname)) 
))
#
ii = which(startsWith(vn,"Pid["))
jj = which(startsWith(vn,"eta["))
UNID_contrib = (mcmc[,jj] * (1-mcmc[,ii])) / rowSums((mcmc[,jj] * (1-mcmc[,ii])))

plot_title = "Contribution to Un-ID prey category"
suppressMessages(print(
mcmc_intervals(UNID_contrib) +
  scale_y_discrete(name="Diet Type",labels=dfPtp$Description) +
  ggtitle(plot_title,subtitle = paste0("Forage data ", Projectname)) 
))

# plots for by-groups:
if(GrpOpt==1){
  
  plot_title = "Posterior distributions, Energy Intake Rate (kcal/min) by Group"
  suppressMessages(print(
    mcmc_areas(mcmc,
               pars = paste0("ERgmn[",seq(1,Ngrp),"]"),
               prob = 0.8) + ggtitle(plot_title,subtitle = paste0("Forage data ", Projectname))  
    + scale_y_discrete(name="Group",labels=Grouplist$Groupname)
  ))
  
  if(length(Grpvar)==1){
    Grpvartxt = Grpvar
  }else{
    Grpvartxt = paste(Grpvar, collapse = ', ')
  }
  
  plot_title = c(paste0("Proportional contribution to diet (by biomass), by ",Grpvartxt,
                      ", Forage data ",Projectname),"","")
  plts = list()
  for(g in 1:min(3,Ngrp)){
    plts[[g]] = suppressMessages(mcmc_intervals(post,
               pars = paste0("PDg[",g,",",seq(1,(NPtypes-1)),"]"),
               prob = 0.8) + 
      scale_y_discrete(name="Diet Type",labels=dfPtp$Description) +
      ggtitle(plot_title[g],subtitle = paste0("Group level ",Grouplist[g,4])) 
    )
  }
  grid.arrange(grobs=plts,ncol=3)
  
  
}

