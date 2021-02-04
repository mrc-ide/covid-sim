#### Author: dlaydon

## Script automates writing parameter files and pre-parameter files, and documents parameters for CovidSim. See covid-sim/docs/inputs-and-outputs.md and covid-sim/docs/intervention-description.md
## Script is a work in progress: if you see a parameter for which you can add to the documentation, please open a PR!
## Broadly, pre-parameter files contain parameters whose values are common to a series of runs. Their values reflect current knowledge of SARS-CoV-2/Covid-19. 
## Parameter files contain parameters for non-pharmaceutical interventions (NPIs), and thus will differ between a series of runs if modelling multiple scenarios.
## Both parameter and pre-parameters have the same format, which is a sequence of:

##		[Description of Parameter]
##		Parameter value

## We have endeavored to name the parameters as they appear in the Cpp code, although there remain occasional exceptions.

## == ## == ## == ## == ## == ## == ## == ## == ## == ## == ## == ## == 
## Issues / still to do: 
#		i) 		more comments for parameters; 
#		ii) 	hard-coded parameters in MakePreParamList (any that need changing between countries/populations please add to arguments accordingly)
# 		iii)	more example parameter files 	

library(here) ## assumes that getwd() will return root "covid-19-spatial-sim" folder.
OutputDir  = here("data/param_files/") ## change as appropriate
invisible(sapply(list.files(here("Rscripts/SourceFiles"), full.names = TRUE), function(x) source(x))) ## Import R source files (inc. MakePreParamList and MakeParamList)


options("scipen" = 13) #### set high penalty for scientific display. (so e.g. 10000 is not outputted as 1e+05) 


#### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== 
#### ==== Import severity progression ICDFs/quantiles (ICDF = Inverse Cumulative Distribution Function) (use script MakeSeverityProgressionICDFs.R) 

ICDFs 		= readRDS(file = file.path(here("data/param_files/"), "SeverityICDFs.rds"))
AS_Means 	= readRDS(file = file.path(here("data/param_files/"), "AS_SeverityProgMeans.rds"))

#### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== 
#### ==== Add in WAIFW matrix using POLYMOD/socialmixr
library(socialmixr)

Ages 			= seq(0,75,5)
POLYMOD_UK 		= contact_matrix(polymod, countries = "United Kingdom", age.limits = Ages)
WAIFW_Matrix 	= POLYMOD_UK$matrix
## contact_matrix only returns as high as 75+. CovidSim wants 75-80 and 80+.
## For now assume that [75,80) = 80+ = 75+. Check this
WAIFW_Matrix = rbind(WAIFW_Matrix, WAIFW_Matrix[nrow(WAIFW_Matrix),]) # add row to WAIFW_Matrix
WAIFW_Matrix = cbind(WAIFW_Matrix, WAIFW_Matrix[,ncol(WAIFW_Matrix)]) # add col to WAIFW_Matrix


#### Write non-US, non-Canada pre-parameter file
PreParamList = MakePreParamList(
		
		#WAIFW_Matrix 				= WAIFW_Matrix,
		#WAIFW_Matrix_SpatialOnly 	= WAIFW_Matrix,
		
		# Means
		LatentPeriod 					= AS_Means$LatentPeriod 		, 
		InfectiousPeriod 				= AS_Means$InfectiousPeriod 	, 
		Mean_MildToRecovery				= AS_Means$Mean_MildToRecovery	,
		Mean_ILIToRecovery 				= AS_Means$Mean_ILIToRecovery ,
		Mean_ILIToSARI 					= AS_Means$Mean_ILIToSARI ,
		Mean_ILIToDeath 				= AS_Means$Mean_ILIToDeath,
		Mean_SARIToCritical 			= AS_Means$Mean_SARIToCritical,
		Mean_SARIToDeath 				= AS_Means$Mean_SARIToDeath,
		Mean_SARIToRecovery 			= AS_Means$Mean_SARIToRecovery,
		Mean_CriticalToDeath 			= AS_Means$Mean_CriticalToDeath,
		Mean_CriticalToStepdownDeath 	= AS_Means$Mean_CriticalToStepdownDeath,
		Mean_CriticalToCritRecov 		= AS_Means$Mean_CriticalToCritRecov,
		Mean_StepdownToDeath 			= AS_Means$Mean_StepdownToDeath,
		
		# quantiles / ICDFs
		latent_icdf 					= ICDFs$latent_icdf 		, 
		infectious_icdf 				= ICDFs$infectious_icdf 	, 
		MildToRecovery_icdf				= ICDFs$MildToRecovery_icdf	, 
		ILIToRecovery_icdf 				= ICDFs$ILIToRecovery_icdf ,
		ILIToSARI_icdf 					= ICDFs$ILIToSARI_icdf ,
		ILIToDeath_icdf 				= ICDFs$ILIToDeath_icdf,
		SARIToCritical_icdf 			= ICDFs$SARIToCritical_icdf,
		SARIToDeath_icdf 				= ICDFs$SARIToDeath_icdf,
		SARIToRecovery_icdf 			= ICDFs$SARIToRecovery_icdf,
		CriticalToDeath_icdf 			= ICDFs$CriticalToDeath_icdf,
		CriticalToStepdownDeath_icdf 	= ICDFs$CriticalToStepdownDeath_icdf, 
		CriticalToCritRecov_icdf 		= ICDFs$CriticalToCritRecov_icdf,
		StepdownToDeath_icdf 			= ICDFs$StepdownToDeath_icdf,
)

WriteParamList(PreParamList, OutputDir = OutputDir, OutputFileName = "pre_R0_2.0_NewICDFs_WAIFWmat.txt")
# make default pre-param file
WriteParamList(MakePreParamList(), OutputDir = OutputDir, OutputFileName = "pre_R0_2.0_test.txt")





#### Parameter files corresponding to following set of interventions:

ArbitrarilyLargeNonStartTime = 100000 ## interventions turned off by specifying a start time greater than the simulation duration, i.e. NumSimulationDays = PreParamList[["Sampling time"]]

# i) 	All interventions
MakeAndWriteParamList(OutputDir = OutputDir, OutputFileName = "p_PC_CI_HQ_SD.txt"	)

# ii) 	No Interventions (sufficient to set start times to be larger than simulation time for interventions not implemented).
MakeAndWriteParamList(OutputDir = OutputDir, OutputFileName = "p_NoInt.txt", 
		CaseIsolationTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		HQuarantineTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		PlaceCloseTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		SocDistTimeStartBase 				= ArbitrarilyLargeNonStartTime) 

# iii) 	Case Isolation Only (sufficient to set start times to be larger than simulation time for interventions not implemented).
MakeAndWriteParamList(OutputDir = OutputDir, OutputFileName = "p_CI.txt", 
		HQuarantineTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		PlaceCloseTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		SocDistTimeStartBase 				= ArbitrarilyLargeNonStartTime)

# iv) 	Case Isolation + Home Quarantine (sufficient to set start times to be larger than simulation time for interventions not implemented).
MakeAndWriteParamList(OutputDir = OutputDir, OutputFileName = "p_CI_HQ.txt", 
		PlaceCloseTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		SocDistTimeStartBase 				= ArbitrarilyLargeNonStartTime)

# v) 	Place Closure only (sufficient to set start times to be larger than simulation time for interventions not implemented).
MakeAndWriteParamList(OutputDir = OutputDir, OutputFileName = "p_PC.txt", 
		CaseIsolationTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		HQuarantineTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		SocDistTimeStartBase 				= ArbitrarilyLargeNonStartTime) 

# vi) 	Social Distancing Only (sufficient to set start times to be larger than simulation time for interventions not implemented).
MakeAndWriteParamList(OutputDir = OutputDir, OutputFileName = "p_SD.txt", 
		CaseIsolationTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		HQuarantineTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		PlaceCloseTimeStartBase 			= ArbitrarilyLargeNonStartTime) 

#### Write parameter files with varying "efficacies" over time. 
MakeAndWriteParamList(OutputDir = OutputDir, OutputFileName = "p_PC_CI_HQ_SDEase_DCT.txt", 
		
		VaryEfficaciesOverTime = 1,
		
		#### Note: numbers here must match Num_SD_ChangeTimes; any times listed in SD_ChangeTimes that are before SocDistTimeStartBase and after SocDistDuration / SocDistDuration_byAdUnit are irrelevant.
		Num_SD_ChangeTimes = 6, SD_ChangeTimes = c(0, 60, 120, 180, 240, 300), ### by default, initialize first time to non-variable time start, then all other times to be arbitrarily large. 
		
		SD_PlaceEffects_OverTime 		= matrix(c(	1, 1, 0.2, 0.2, 
													1, 1, 0.3, 0.3, 
													1, 1, 0.4, 0.4, 
													1, 1, 0.5, 0.5, 
													1, 1, 0.6, 0.6, 
													1, 1, 0.7, 0.7), nrow = 6, byrow = TRUE),
		SD_HouseholdEffects_OverTime 	= c(1.2, 1.2, 1.15, 1.15, 1.1, 1.0), 
		SD_SpatialEffects_OverTime		= c(0.1, 0.2, 0.3, 0.4, 0.5, 0.5)	)	



