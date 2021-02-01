# SOFA (Sea Otter Forage Analysis), Version 2.0
# Step 3, Model Summary.
# Source this script to generate an html markdown report of results. 
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
require(ggplot2)
require(gridExtra)
require(bayesplot)
require(rmarkdown)
require(gtools)
require(dplyr)
require(knitr)
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
rdata_file = rselect.list(rslt_list, preselect = NULL, multiple = FALSE,
										title = "Select results file" ,
										graphics = getOption("menu.graphics")) 
if(length(rdata_file)==0){
	dlg_message(c("No data file selected"), "ok")
	stop_quietly()
}
resultsfilename = paste0("./projects/",Projectname,"/results/",rdata_file)
if(length(grep("Grp", rdata_file))>0){Grp_TF=TRUE}else{Grp_TF=FALSE}
file.copy(resultsfilename,
					paste0("./code/","Results.rdata"),overwrite = TRUE)
#
GL1 = 1; GL2 = 0; GL3 = 0
if(Grp_TF){
	attach("./code/Results.rdata"); Grouplist <- Grouplist; detach("file:./code/Results.rdata")
	dlg_message(c("This results file has by-groups. You can select up to 3 sample group levels ",
								"for generating certain plots that show statistics by group level.",
								"(Tables of stats will still be generated for all group levels)"), "ok")
	Grplevsamp = rselect.list(Grouplist$Groupname, preselect = NULL, multiple = TRUE,
														title = "Select group levels" ,
														graphics = getOption("menu.graphics")) 
	Grplevsamp = Grplevsamp[1:min(3,length(Grplevsamp))]
	GL1 = which(Grouplist$Groupname==Grplevsamp[1])
	if(length(Grplevsamp)>1){GL2 = which(Grouplist$Groupname==Grplevsamp[2])}else{GL2=0}
	if(length(Grplevsamp)>2){GL3 = which(Grouplist$Groupname==Grplevsamp[3])}else{GL3=0}
	rm(Grplevsamp)
}
rspnse = dlg_message(c("About to render an html rmarkdown report summarizing model results, ",
											 "which could take a few minutes to complete. ",
											 "Continue?"), "okcancel")$res
if(rspnse == "cancel"){
	stop_quietly()
}
vers = read.csv("./code/Version.csv"); vers = vers$Version
title = paste0("~~~ Sea otter foraging analysis (SOFA) ",vers," ~~~ ", "Project: ", Projectname)
Daterun = Sys.Date()
render("./code/SOFA_summary.Rmd",
			 output_dir = paste0("./projects/",Projectname),
			 output_file = "SOFA_summary.html",
			 params = list(rep_title = title, rep_date = Daterun, show.grptxt = Grp_TF,
			 							GL1 = GL1, GL2=GL2, GL3=GL3)) # 
dlg_message(c("The results can be viewed by opening 'SOFA_summary.html' in the projecy folder"), "ok")
