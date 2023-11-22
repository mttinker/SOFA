# SOFA (Sea Otter Forage Analysis), Version 2.0
# Step 3, Model Summary.
# Source this script to generate an html markdown report of results. 
require(svDialogs)
require(rstudioapi)
require(readxl)
require(openxlsx)
require(tcltk)
require(parallel)
require(rstan)
require(ggplot2)
require(gridExtra)
require(bayesplot)
require(rmarkdown)
require(gtools)
require(dplyr)
require(knitr)
require(kableExtra)
rm(list = ls())
# Create Generic function for stopping script in case of error:
stop_quietly <- function() {
	opt <- options(show.error.messages = FALSE)
	on.exit(options(opt))
	stop()
}
rspnse = dlg_message(c("This script is used to generate an html report that summarizes results of  ",
											 "the analysis of sea otter foraging data using SOFA (Sea Otter Foraging Analysis). ",
											 "You begin by selecting a Project name and results file. ",
											 "Continue?"), "okcancel")$res
if(rspnse == "cancel"){
	stop_quietly()
}
Sys.sleep(.5)
proj_list = dir('./projects')
Projectpath = dlg_dir('./projects/', "Select Project folder")$res
Projectname = basename(Projectpath)
Sys.sleep(.5)
rslt_list = dir(paste0(Projectpath,"/results"))
rdata_file = tk_select.list(rslt_list, preselect = NULL, multiple = FALSE,
										title = "Select results file" ) 
if(length(rdata_file)==0){
	dlg_message(c("No data file selected"), "ok")
	stop_quietly()
}
resultsfilename = paste0("./projects/",Projectname,"/results/",rdata_file)
if(length(grep("Grp", rdata_file))>0){Grp_TF=TRUE}else{Grp_TF=FALSE}
file.copy(resultsfilename,
					# "Results.rdata",overwrite = TRUE)
					paste0("./code/","Results.rdata"),overwrite = TRUE)
#
GL1 = 1; GL2 = 0; GL3 = 0; GL4 = 1; GL5 = 0; GL6 = 0; GL7 = 1; GL8 = 0; GL9 = 0; 
GrpFgHt = 5

if(Grp_TF){
	# attach("Results.rdata"); Grouplist <- Grouplist; detach("file:Results.rdata")
	attach("./code/Results.rdata"); Grouplist <- Grouplist; detach("file:./code/Results.rdata")
	dlg_message(c("This results file has by-groups. You can select up to 9 sample group levels ",
								"for generating certain plots that show statistics by group level.",
								"(Tables of stats will still be generated for all group levels)"), "ok")
	Grplevsamp = tk_select.list(Grouplist$Groupname, preselect = NULL, multiple = TRUE,
														title = "Select group levels" ) 
	Grplevsamp = Grplevsamp[1:min(9,length(Grplevsamp))]
	GL1 = which(Grouplist$Groupname==Grplevsamp[1])
	if(length(Grplevsamp)>1){GL2 = which(Grouplist$Groupname==Grplevsamp[2])}else{GL2=0}
	if(length(Grplevsamp)>2){GL3 = which(Grouplist$Groupname==Grplevsamp[3])}else{GL3=0}
	if(length(Grplevsamp)>3){GL4 = which(Grouplist$Groupname==Grplevsamp[4])}else{GL4=0}
	if(length(Grplevsamp)>4){GL5 = which(Grouplist$Groupname==Grplevsamp[5])}else{GL5=0}
	if(length(Grplevsamp)>5){GL6 = which(Grouplist$Groupname==Grplevsamp[6])}else{GL6=0}
	if(length(Grplevsamp)>6){GL7 = which(Grouplist$Groupname==Grplevsamp[7])}else{GL7=0}
	if(length(Grplevsamp)>7){GL8 = which(Grouplist$Groupname==Grplevsamp[8])}else{GL8=0}
	if(length(Grplevsamp)>8){GL9 = which(Grouplist$Groupname==Grplevsamp[9])}else{GL9=0}
	
	if(length(Grplevsamp)<4){GrpFgHt = 5}
	if(length(Grplevsamp)>3 & length(Grplevsamp)<7){GrpFgHt = 10}
	if(length(Grplevsamp)>6){GrpFgHt = 15}

	rm(Grplevsamp)
}
rspnse = dlg_message(c("About to render an html rmarkdown report summarizing model results, ",
											 "which could take a few minutes to complete. ",
											 "Continue?"), "okcancel")$res
if(rspnse == "cancel"){
	stop_quietly()
}
vers = read.csv("./code/Version.csv"); vers = vers$Version
title = paste0("~~~ Sea otter foraging analysis (SOFA) ",vers," ~~~ ")
subtitle = paste0("Project: ", Projectname, ", Results file: ",rdata_file)
Daterun = Sys.Date()
render("./code/SOFA_summary.Rmd",
			 output_dir = paste0("./projects/",Projectname),
			 output_file = "SOFA_summary.pdf",
			 params = list(rep_title = title, rep_subtitle = subtitle, 
			 							rep_date = Daterun, show.grptxt = Grp_TF,
			 							GL1 = GL1, GL2=GL2, GL3=GL3,
			 							GL4 = GL4, GL5=GL5, GL6=GL6,
			 							GL7 = GL7, GL8=GL8, GL9=GL9,GrpFgHt=GrpFgHt)) # 
dlg_message(c("The results can be viewed by opening 'SOFA_summary' in the projecy folder"), "ok")
