# SOFA (Sea Otter Forage Analysis), Version 3.1
# Step 3, Model Summary.
# Source this script to generate an html markdown report of results. 
#
rm(list = ls())
#
existing_Packages<-as.list(installed.packages()[,1])
required_Packages<-c('svDialogs', 'rstudioapi', 'readxl', 'openxlsx','tcltk', 'parallel', 'rstan', 
                     'ggplot2', 'gridExtra', 'bayesplot', 'rmarkdown', 'gtools', 'dplyr', 'shiny',
                     'knitr', 'kableExtra')
missing_Packages<- required_Packages[!required_Packages %in% existing_Packages]
if(length(missing_Packages)>0)install.packages(pkgs =  missing_Packages)
invisible(lapply(required_Packages, require, character.only=T,quietly = T))
rm(existing_Packages,missing_Packages,required_Packages)
#
# Create Generic function for stopping script in case of error:
stop_quietly <- function() {
	opt <- options(show.error.messages = FALSE)
	on.exit(options(opt))
	stop()
}
rspnse = dlg_message(c("This script is used to generate an html report that summarizes results of  ",
											 "the analysis of sea otter foraging data using SOFA (Sea Otter Foraging Analysis). ",
											 "You begin by selecting a Project name and then the results file. ",
											 "Continue?"), "okcancel")$res
if(rspnse == "cancel"){
	stop_quietly()
}
Sys.sleep(.5)
proj_list = dir('./projects')
Prjct_path = dlg_dir('./projects/', "Select Project folder")$res
Prjct_name = basename(Prjct_path)
Sys.sleep(.5)
base_dir = getwd()
setwd(paste0(Prjct_path,"/results"))
tmp = choose.files(filters = Filters[c("RData"),],multi = FALSE,caption= "Select results file")
setwd(base_dir)
rdata_file = basename(tmp)
rdata_file_noextns = substr(rdata_file, 1, nchar(rdata_file)-6) 
#
if(length(rdata_file)==0){
	dlg_message(c("No data file selected"), "ok")
	stop_quietly()
}
resultsfilename = paste0("./projects/",Prjct_name,"/results/",rdata_file)
if(length(grep("Grp", rdata_file))>0){Grp_TF=TRUE}else{Grp_TF=FALSE}
# Remove temporary "Results.rdata" file first (in case overwrite doesn't work)
if (file.exists(paste0("./code/","Results.rdata"))) {
  #Delete file if it exists
  file.remove(paste0("./code/","Results.rdata"))
}
file.copy(resultsfilename,
					paste0("./code/","Results.rdata"),overwrite = TRUE)
#
color_scheme_set("blue")
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
subtitle = paste0("Project: ", Prjct_name, ", Results file: ",rdata_file)
rslts_dirname = paste0("./projects/",Prjct_name,"/results")
rslts_filename = paste0("SOFA_sum_",rdata_file_noextns,".pdf") 
rmd_pathname = "./code/SOFA_summary_v3.Rmd"
tmp_dir <- tempdir()
Daterun = Sys.Date()
if (file.exists(file.path(tmp_dir,rslts_filename))) {
  #Delete file if it exists
  file.remove(file.path(tmp_dir,rslts_filename))
}
#
render(rmd_pathname,
			 output_dir = tmp_dir, # NOTE: need to compile to tmp_dir to avoid pathname with spaces (#*^# Microsoft!)
			 output_file = rslts_filename,
			 params = list(rep_title = title, rep_subtitle = subtitle, 
			 							rep_date = Daterun, show.grptxt = Grp_TF,
			 							GL1 = GL1, GL2 = GL2, GL3 = GL3, GL4 = GL4, 
			 							GrpFgHt = GrpFgHt)) # 
#
# Remove existing summary pdf first (in case overwrite doesn't work)
if (file.exists(file.path(rslts_dirname,rslts_filename))) {
  #Delete file if it exists
  file.remove(file.path(rslts_dirname,rslts_filename))
}

file.copy(file.path(tmp_dir,rslts_filename),rslts_dirname, overwrite = T)

on.exit(unlink(tmp_dir))

write.xlsx(sumstats, file.path(rslts_dirname,paste0(rdata_file_noextns,"_sumstats.xlsx")),
           colNames = TRUE, rowNames = TRUE)

dlg_message(paste("Results can be viewed in the SOFA_sum_(results_filename).pdf in results sub-folder of project folder. ",
                  "A table of summary stats can also be found in a spreadsheet with corresponding name"), "ok")

