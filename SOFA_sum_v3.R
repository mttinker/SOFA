# SOFA (Sea Otter Forage Analysis), Version 3.0
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

GrpFgHt = 10
GL1 = 0; GL2 = 0; GL3 = 0; GL4 = 0; 
if(Grp_TF){
	iiog = numeric()
	# attach("Results.rdata"); Grouplist <- Grouplist; detach("file:Results.rdata")
	attach("./code/Results.rdata"); Grouplist <- Grouplist; detach("file:./code/Results.rdata")
	dlg_message(c("This results file has by-groups. Next you can select 2 - 4 sample group levels ",
								"for generating plots comparing prey-level statistics by group level.",
								"(Tables of group-specific parameters will be generated for all group levels)"), "ok")
	Grplevsamp = tk_select.list(Grouplist$Groupname, preselect = NULL, multiple = TRUE,
														title = "Select group levels" ) 
	if(length(Grplevsamp)==0){ # In case user selects 0 groups, select the first 2 or 3 for them...
		for(g in 1:min(3,length(Grouplist$Groupname))){
			if(g==1){
				GL1 = 1
			}else if(g==2){
				GL2 = 2
			}else if(g==3){
				GL3 = 3
			}
		}
	}else if(length(Grplevsamp)==1){ # In case user selects only 1 group, select some others...
		for(g in 1:min(3,length(Grouplist$Groupname))){
			iiog = which(Grouplist$Groupname != Grplevsamp[1])
			if(g==1){
				GL1 = which(Grouplist$Groupname==Grplevsamp[1])
			}else if(g==2){
				GL2 = iiog[1]
			}else if(g==3){
				GL3 = iiog[2]
			}
		}
	}else{ # or, assuming the user actually follows instructions, use the first 4 they select...
		for(g in 1:min(4,length(Grplevsamp))){
			if(g==1){
				GL1 = which(Grouplist$Groupname==Grplevsamp[1])
			}else if(g==2){
				GL2 = which(Grouplist$Groupname==Grplevsamp[2])
			}else if(g==3){
				GL3 = which(Grouplist$Groupname==Grplevsamp[3])
			}else if(g==4){
				GL4 = which(Grouplist$Groupname==Grplevsamp[4])
			}
		}
	}
	rm(Grplevsamp,iiog)
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
output_dirname =  paste0("./projects/",Projectname)
output_filename = "SOFA_summary.pdf"
rmd_pathname = "./code/SOFA_summary_v3.Rmd"
tmpdir <- tempdir()
Daterun = Sys.Date()

render(rmd_pathname,
			 output_dir = tmpdir, # NOTE: need to compile to tmpdir to avoid pathname with spaces (#*^# Microsoft!)
			 output_file = output_filename,
			 params = list(rep_title = title, rep_subtitle = subtitle, 
			 							rep_date = Daterun, show.grptxt = Grp_TF,
			 							GL1 = GL1, GL2 = GL2, GL3 = GL3, GL4 = GL4, 
			 							GrpFgHt = GrpFgHt)) # 

file.copy(file.path(tmpdir,output_filename),output_dirname, overwrite = T)

on.exit(unlink(tmpdir))

dlg_message(c("The results can be viewed by opening 'SOFA_summary' in the projecy folder"), "ok")
