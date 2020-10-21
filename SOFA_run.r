# SOFA (Sea Otter Forage Analysis), Version 2.0
# Step 2, Model fitting.
# Source this script to fit the SOFA model to sea otter foraging data. 
#
require(svDialogs)
require(rstudioapi)
require(readxl)
require(openxlsx)
# NOTE: necessary to install 64 bit version of Java ('Windows Offline')
#  from https://www.java.com/en/download/manual.jsp
require(rJava) 
require(rChoiceDialogs)
require(parallel)
require(rstan)
#
# Create Generic function for stopping script in case of error:
stop_quietly <- function() {
	opt <- options(show.error.messages = FALSE)
	on.exit(options(opt))
	stop()
}
rspnse = dlg_message(c("This script is used to set up and run SOFA, fitting the model ",
											 "to a data project that has been prepared by the SOFA_prep script. ",
											 "If you have not already prepared a project, cancel and do that now. ",
											 "Otherwise, you will now be able to select a project to work with, ",
											 "set the necessary user settings and select from various options, ",
											 "and then run the model fitting (which takes a long time). ", 
											 "Once complete you can view results using the 'SOFA_sum' script. ",
											 "Continue?"), "okcancel")$res
if(rspnse == "cancel"){
	stop_quietly()
}
proj_list = dir('./projects')
Projectpath = dlg_dir('./projects/', "Select Project folder")$res
Projectname = basename(Projectpath)
Sys.sleep(.5)
# Projectname = rselect.list(proj_list, preselect = NULL, multiple = FALSE,
# 						 title = "Select Project name",
# 						 graphics = getOption("menu.graphics"))
#
MnN1 = as.numeric(dlg_input(message = "Min dives/bout for estimating prey attributes", 
														default = 3)$res)
MnN2 = as.numeric(dlg_input(message = "Min dives/bout for estimating effort allocation", 
														default = 2)$res)
nsamples = as.numeric(dlg_input(message = "Number posterior samples from Bayesian fitting", 
														default = 10000)$res)
nburnin = as.numeric(dlg_input(message = "Number of burn0in samples for Bayesian fitting", 
														default = 750)$res)
# Process data---------------------------------------------------------------
# Load data (after selecting Projectname from available sub-folders of data folder)
dat = read_excel(paste0("./projects/",Projectname,"/Forage_dat.xlsx"))
# Load reference spreadsheets
dfPr = read_excel(paste0("./projects/",Projectname,"/Prey_Key.xlsx"))
dfPtp = read_excel(paste0("./projects/",Projectname,"/Prey_Key.xlsx"),sheet = "Prey_Types")
dfPcl = read_excel(paste0("./projects/",Projectname,"/Prey_Key.xlsx"),sheet = "Prey_Classes")
dfSz = read_excel(paste0("./projects/",Projectname,"/Sizeclass_key.xlsx"))
dfPaw = read_excel(paste0("./projects/",Projectname,"/Pawsz.xlsx"))
dfM = read_excel(paste0("./projects/",Projectname,"/Mass_lng_dat.xlsx"))
dfE = read_excel(paste0("./projects/",Projectname,"/Energy_dat.xlsx"))
dfElst = read_excel(paste0("./projects/",Projectname,"/Prey_Spcs.xlsx"))
if(length(!is.na(dfPr$PreyType))<length(!is.na(dfPr$PreyCode))){
	dlg_message(c("The Prey_Key spreadsheet has not been completed properly, ",
								"please go back and edit. Stopping program."))
	stop_quietly()
}
#  Allow user to select a subset of data for analysis based on values of a group variable
rspns = dlg_message(c("You can either analyze ALL the data, or just a sub-set of records ",
											 "based on the values of one of the categorical variables. ",
											 "Note that if ALL, you will have the option of setting by-groups, ",
											 "while if you wish to analyze just a sub-set of records ", 
											 "then by-groups will not be an option"), "ok")$res
slctvars = c("Area", "Site", "Period", "Sex", "Ageclass", "Ottername", "Pup")
# Note: determine which of these potential selection variables has enough levels
ix = numeric()
for (i in 1:length(slctvars)){
	col = which(colnames(dat)==slctvars[i])
	lvls = as.matrix(table(dat[,col]))
	if(length(lvls[which(lvls[,1]>100),1])<2){
		ix = c(ix,i)
	} 
}
slctvars = slctvars[-ix]
#
rspns = dlg_message(c("Do you wish to analyze just a sub-set of the data", 
											"(if 'yes', you will next be asked to select data ", 
											"based on values of a group variable)?"), "yesno")$res
if (rspns=="no"){
	SBST = 0
}else{
	SBST = 1
}
if (SBST==1){
	slctvar = rselect.list(slctvars, preselect = NULL, multiple = FALSE,
 						 title = "Which variable?",
 						 graphics = getOption("menu.graphics"))
	col = which(colnames(dat)==slctvar)
	lvls = rownames(as.matrix(table(dat[,col])))
	lvls = rselect.list(lvls, preselect = NULL, multiple = TRUE,
										 title = paste0("Which values of ",slctvar) ,
										 graphics = getOption("menu.graphics")) 
	iis = numeric()
	slctvardef = paste0(slctvar,"-")
	for(i in 1:length(lvls)){
		slctvardef = paste0(slctvardef,lvls[i])
		ii = which(dat[,col]==lvls[i])
		iis = c(iis,ii)
	}
	dat = dat[iis,]
}
#
source("./code/Fdatprocess.R")
errtxt = c("An error occured in one or more of the data files, ",
					 "check that all data fields are formatted followig the manual,",
					 "and valid proxy species have been selected, then try again.")
wrntxt = c("Some mass-length functions are poor due to limited data")
test = tryCatch({
	rslt <- Fdatprocess(dat,dfPr,dfPtp,dfPcl,dfSz,dfPaw,dfM,dfE,dfElst)
}, warning = function(war) {
	dlg_message(wrntxt)
	print("Review diagnistic warnings and adjust options as needed.")
	print(paste("Fdatprocess warnings: ",war))
#   stop_quietly()
}, error = function(err) {
	dlg_message(errtxt)
	stop_quietly()
}, finally = {
	print("Fdatprocess succesful")
}) # END try-catch loop
rm(test)
Fdat = rslt$Fdat; Boutlist=rslt$Boutlist
Cal_dns_mn = rslt$Cal_dns_mn; Cal_dns_sg = rslt$Cal_dns_sg
logMass_sg = rslt$logMass_sg
Nbouts = rslt$Nbouts; Ndives = rslt$Ndives; Nobs = rslt$Nobs; NPtypes = rslt$NPtypes
MassLngFits = rslt$MassLngFits
rm(rslt)
# Generate mass-lenth plots by prey code
for(pr in 1:nrow(dfPr)){
	if(dfPr$PreyType[pr] != "UNID"){
		ft = MassLngFits[[pr]]
		plot(exp(ft$model$x),exp(ft$model$y),main=dfPr$Description[pr],
				 xlab="size (mm)",ylab="mass (g)")
		lines(exp(ft$prdct$x),exp(ft$prdct$ypred),col="red")
	}
}
#  Allow user to set by-groups for analysis based on grouping variables
if (SBST == 0){
	rspns = dlg_message(c("You have the option of setting 'by-groups' for analysis, ",
												"where-by one or more categorical variables are used to ", 
												"divide the data into groups, and stats will be generated ",
												"for each group level."), "ok")$res
	rspns = dlg_message(c("Do you wish to analayze data by group, with the by-group",
												"levels determined by one or more categorical variables?"), "yesno")$res
	if (rspns=="no"){
		GrpOpt = 0
	}else{
		GrpOpt = 1
	}
	if(GrpOpt==1){
		Grpvar = rselect.list(slctvars, preselect = NULL, multiple = TRUE,
												 title = "Which variables?",
												 graphics = getOption("menu.graphics"))
		Ngrpvar = length(Grpvar)
	}
}else{
	GrpOpt=0
}
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
source("./code/Boutprocess.R")
stan.data <- Boutprocess(Fdat,Nbouts,NPtypes,Boutlist,MnN1,MnN2,GrpOpt,Ngrp)
#
# Fit Bayesian model ------------------------------------------------------
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
testfile = paste0("./code/",substr(fitmodel ,1,nchar(fitmodel)-5),"_dso.rdata")
if(file_test("-f", testfile) ==TRUE){
	load(testfile)
}else{
	Stan_model.dso = stan_model(file = paste0("./code/",fitmodel),
															model_name = substr(fitmodel ,1,nchar(fitmodel)-5))
	save(Stan_model.dso,file=testfile)
}
#
rspnse = dlg_message(c("That completes set-up, the model will now be fit using ",
											 "Bayesian methods (STAN). Fitting can take several hours, ",
											 "depending on size of the data file and speed of the computer. ",
											 "Continue?"), "okcancel")$res
if(rspnse == "cancel"){
	stop_quietly()
}
#
errtxt = c("An error occured in in model fitting. Check to make sure STAN",
					 "and other libraries are properly installed, check stan_data for NA ",
					 "values or other errors, and then re-start R and try again.")
test = tryCatch({
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
}, error = function(err) {
	dlg_message(errtxt)
	stop_quietly()
}, warning = function(war) {
	print("Review diagnistic warnings from fitting and adjust options as needed.")
	print(paste("STAN WARNINGS:  ",war))
	#   stop_quietly()
}, finally = {
	print("Fitting Succesful")
}) # END try-catch loop
rm(test)
mcmc <- as.matrix(out)
vn = colnames(mcmc)
Nsims = nrow(mcmc)
sumstats = summary(out)$summary
vns = row.names(sumstats)
rm(params,Stan_model.dso)
# save image file of results for post-fit processing and review:
dir.create(paste0("./projects/",Projectname,"/results"),showWarnings = F)
#
if (GrpOpt==0){
	if (SBST==1){
		save.image(paste0("./projects/",Projectname,"/results/Rslt_",slctvardef,"_",
											format(Sys.time(), "%Y_%b_%d_%H"),"hr.rdata"))
	}else{
		save.image(paste0("./projects/",Projectname,"/results/Rslt_All_", 
											format(Sys.time(), "%Y_%b_%d_%H"),"hr.rdata"))
	}
}else{
	save.image(paste0("./projects/",Projectname,"/results/Rslt_Grp_", 
										format(Sys.time(), "%Y_%b_%d_%H"),"hr.rdata"))
}
fintxt = c("That completes model fitting, check psrf values and other diagnostics. ",
					 "The results have been saved to the 'results' sub-folder of the project. ",
					 "You can view summary plots and tables by running the 'SOFA_sum' script. ")
dlg_message(fintxt)
