# SOFA (Sea Otter Forage Analysis), Version 2.0
# Step 1, Preparation.
# Source this script to prepare for an analysis of sea otter foraging data. 
require(svDialogs)
require(rstudioapi)
require(readxl)
require(openxlsx)
# Generic function for stopping script in case of error:
stop_quietly <- function() {
	opt <- options(show.error.messages = FALSE)
	on.exit(options(opt))
	stop()
}
#
rspnse = dlg_message(c("This script is used to select a data file for running SOFA ",
											 "and set up a project folder with necessary helper files. ",
											 "You will be asked to select a properly-formatted spreadsheet ",
											 "of raw data to analyze (refer to manual) and create a Project name. ",
											 "Then you can go to the project directory (in 'projects' folder)",
											 "and edit the 'Prey_key' spreadsheet, and other files as desired, ", 
											 "before running the SOFA model fitting step (SOFA_run).",
											 "Continue?"), "okcancel")$res
if(rspnse == "cancel"){
	stop_quietly()
}
#
Sys.sleep(.5)
Datafile <-selectFile(caption = "Select a spreadsheet of raw forage data formatted correctly",
																	 filter = c("Excel files (*.xlsx)"),
																	 existing = TRUE)
#
dat0 = read_excel(Datafile)
colnames(dat0) = tolower(colnames(dat0))
# Convert data to a standardized data frame (also catch errors in formatting)
errtxt = c("The selected data file was not formatted as expected, please ",
								"refer to the manual and check that columns and fieldnames ",
								"are properly formatted, then try again.")
rslt = tryCatch({
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
}, warning = function(war) {
	dlg_message(errtxt)
	stop_quietly()
}, error = function(e) {
	dlg_message(errtxt)
	stop_quietly()
}, finally = {rm(dat0)}) # END try-catch loop
rm(rslt)
Sys.sleep(.5)
# Prompt user for a name for the project:
Projectname <- dlg_input("Enter a name for this project", "Edit_me")$res
Sys.sleep(.5)
#
# Create projects directory if it does not exist
dir.create(paste0("./projects"),showWarnings = F)
# Create a folder with the project name specified,
# and if it already exists, alert user and allow to abort or over-write
dir.create(paste0("./projects/",Projectname),showWarnings = F)
files = dir(paste0("./projects/",Projectname,"/"))
del = which(files=="oldfiles")
if(length(del)>0){
	files = files[-del]
}
if(length(files)>0){
	rspnse = dlg_message(c("There are already some files in a folder with this project name. ",
												 "Do you wish to over-write them with new data?  If you select 'OK', ",
												 "the existing files will be transferred to an 'oldfiles' folder.",
												 "Continue?"), "okcancel")$res
	if(rspnse == "cancel"){
		stop_quietly()
	}else{
		"Create backup filder and move existing files to it"
		dir.create(paste0("./projects/",Projectname,"/oldfiles"),showWarnings = F)
		file.copy(file.path(paste0("./projects/",Projectname),files),
							paste0("./projects/",Projectname,"/oldfiles"),overwrite = TRUE)
		file.remove(file.path(paste0("./projects/",Projectname),files))
	}
}
dat = dat[order(dat$Bout,dat$Subbout,dat$Divenum,dat$Success,dat$Prey,dat$Size),]
# Write data file to specified project directory:
write.xlsx(dat,file=paste0("./projects/",Projectname,"/Forage_dat.xlsx"))
# Copy other helper files to directory, 
files = c("Sizeclass_key.xlsx","Pawsz.xlsx")
file.copy(file.path(paste0("./helperfiles"),files),
					paste0("./projects/",Projectname),overwrite = TRUE)
# Create Mass_lng_dat and Energy_dat and Prey_Spcs spreadsheets
dfPL_Spc = read_excel("./preylib/Species.xlsx")
dfPL_MssLn = read_excel("./preylib/RawSizeMassData.xlsx")
dfPL_Enrgy = read_excel("./preylib/Energy.xlsx")
iis = numeric()
Ppn_ed = numeric()
iims = numeric()
iien = numeric()
PryctM = character()
PryctE = character()
SpcnmM = character()
SpcnmE = character()
for (i in 1:nrow(dfPL_Spc)){
	sp = dfPL_Spc$SppCode[i]
	iim = which(dfPL_MssLn$SppCode==sp & !is.na(dfPL_MssLn$MaxLinearDim_mm) & 
								!is.na(dfPL_MssLn$TotWetMass))
	iie = which(dfPL_Enrgy$SppCode==sp & !is.na(dfPL_Enrgy$kcal_g_edblwet))
	if(length(iim) >= 4 & length(iie)>0){
		iis = c(iis,i)
		iims = c(iims,iim)
		iien = c(iien,iie)
		PryctM = c(PryctM,rep(dfPL_Spc$PreyCat1[i],length(iim)))
		SpcnmM = c(SpcnmM,rep(dfPL_Spc$SpeciesName[i],length(iim)))
		PryctE = c(PryctE,rep(dfPL_Spc$PreyCat1[i],length(iie)))
		SpcnmE = c(SpcnmE,rep(dfPL_Spc$SpeciesName[i],length(iie)))		
		Ppn_ed = c(Ppn_ed,mean(dfPL_MssLn$Prpn_Edible[iim],na.rm = T))
	}
}
Prey_Spcs = data.frame(PreyCat1=dfPL_Spc$PreyCat1[iis],
											 SppCode = dfPL_Spc$SppCode[iis],
											 SpeciesName = dfPL_Spc$SpeciesName[iis],
											 AvgOfPrpn_Edible = Ppn_ed)
Mass_lng_dat = data.frame(PreyCat1 = PryctM,
													SpeciesName = SpcnmM,
													SppCode = dfPL_MssLn$SppCode[iims],
													Region = dfPL_MssLn$Region[iims],
													MaxLinearDim_mm = dfPL_MssLn$MaxLinearDim_mm[iims],
													TotWetMass = dfPL_MssLn$TotWetMass[iims])
Energy_dat = data.frame(PreyCat1 = PryctE,
												SpeciesName = SpcnmE,
												SppCode = dfPL_Enrgy$SppCode[iien],
												Region = dfPL_Enrgy$Region[iien],
												N = dfPL_Enrgy$N[iien],
												kcal_g_edblwet = dfPL_Enrgy$kcal_g_edblwet[iien],
												kcal_g_edblwet_SD = dfPL_Enrgy$kcal_g_edblwet_SD[iien])
write.xlsx(Prey_Spcs,file=paste0("./projects/",Projectname,"/Prey_Spcs.xlsx"))
write.xlsx(Mass_lng_dat,file=paste0("./projects/",Projectname,"/Mass_lng_dat.xlsx"))
write.xlsx(Energy_dat,file=paste0("./projects/",Projectname,"/Energy_dat.xlsx"))
#
# Set up Prey_Key spreadsheet for editing:
wb = createWorkbook()
addWorksheet(wb, "Prey_code_key") # First worksheet = prey codes
PreyCodes = unique(dat$Prey); PreyCodes = PreyCodes[order(PreyCodes)]
PreyCodes = PreyCodes[which(!is.na(PreyCodes))]
Npcode = length(PreyCodes)
N_obsprey = numeric()
for(i in 1:Npcode){
	N_obsprey[i] = length(which(dat$Prey==PreyCodes[i]))
}
prey_df = data.frame("PreyCode" = PreyCodes, 
										 "SzBrkmm" = rep("",Npcode),
										 "Description" = rep("",Npcode),
										 "N" = N_obsprey,
										 "PreyType" = rep("",Npcode))
for(i in 1:20){
	tmp = data.frame("Prox" = rep("",Npcode))
	colnames(tmp) = paste0("Prox",i)
	prey_df = cbind(prey_df,tmp)
}
writeData(wb, sheet = "Prey_code_key", x = prey_df, startCol = 1)
# Prey_Types
addWorksheet(wb, "Prey_Types")
prey_types_df = data.frame("TypeN" = 0, "PreyType" = "UNID",
													 "Description" = "UN-IDENTIFIED","Class"="")
writeData(wb, sheet = "Prey_Types", x = prey_types_df, startCol = 1)
# Add drop-down validation for prey types"
dataValidation(wb, "Prey_code_key", col = 5, rows = 2:500, type = "list", 
							 value = "'Prey_Types'!$B$2:$B$50")
# Prey_Classes
prey_class_df = read_excel("./helperfiles/Prey_Classes.xlsx")
addWorksheet(wb, "Prey_Classes")
writeData(wb, sheet = "Prey_Classes", x = prey_class_df, startCol = 1)
# Add drop-down validation for prey classes"
dataValidation(wb, "Prey_Types", col = 4, rows = 2:500, type = "list", 
							 value = "'Prey_Classes'!$B$2:$B$50")
# Prey_Library_Codes
addWorksheet(wb, "Prey_Library_Codes")
writeData(wb, sheet = "Prey_Library_Codes", x = Prey_Spcs[,1:3], startCol = 1)
# Add drop-down validation for proxy species"
dataValidation(wb, "Prey_code_key", col = 6:26, rows = 2:500, type = "list", 
							 value = "'Prey_Library_Codes'!$B$2:$B$50")
# Add comments for user:
Comments_df = read_excel("./helperfiles/Prey_Key_Comments.xlsx")
Clst = list()
for(i in 1:nrow(Comments_df)){
	nc = nchar(Comments_df$Comment[i])
	Clst[[i]] = createComment(comment = Comments_df$Comment[i],visible = F,
														width = 4,height = round(nc/25))
	if (i < 7){
		writeComment(wb, "Prey_code_key", col = i, row = 1, comment = Clst[[i]])
	}else if(i==7){
		writeComment(wb, "Prey_Types", col = 1, row = 1, comment = Clst[[i]])
	}else if(i==8){
		writeComment(wb, "Prey_Types", col = 2, row = 1, comment = Clst[[i]])
	}else if(i==9){
		writeComment(wb, "Prey_Types", col = 4, row = 1, comment = Clst[[i]])
	}else if(i==10){
		writeComment(wb, "Prey_Library_Codes", col = 2, row = 1, comment = Clst[[i]])
	}
}
saveWorkbook(wb, file=paste0("./projects/",Projectname,"/Prey_Key.xlsx"),overwrite = TRUE)
# Closing comment
rspnse = dlg_message(c("That concludes the data and project preparation step. ",
	                     "You should now go to the project folder to edit the Prey Key ",
											 "and any of the other 'helper files', as appropriate. ",
											 "Once the Prey_Key is complete, you can run the SOFA analysis ", 
											 "by Sourcing the 'SOFA_run' script"), "ok")$res
