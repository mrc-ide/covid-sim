#### Author: dlaydon

## Script automates writing parameter files and pre-parameter files, and documents parameters for CovidSim. See covid-sim/docs/inputs-and-outputs.md and covid-sim/docs/intervention-description.md
## Script is a work in progress: if you see a parameter for which you can add to the documentation, please open a PR!
## Broadly, pre-parameter files contain parameters whose values are common to a series of runs. Their values reflect current knowledge of SARS-CoV-2/Covid-19. 
## Parameter files contain parameters for non-pharmaceutical interventions (NPIs), and thus will differ between a series of runs if modelling multiple scenarios.
## Both parameter and pre-parameters have the same format, which is a sequence of:

##		[Description of Parameter]
##		value

## We have endeavored to name the parameters as they appear in the Cpp code, although there remains the occasional exception.

## == ## == ## == ## == ## == ## == ## == ## == ## == ## == ## == ## == 
## Issues / still to do: 
#		i) 		more comments for parameters; 
#		ii) 	hard-coded parameters in MakePreParamList (any that need changing between countries/populations please add to arguments accordingly)
# 		iii)	more example parameter files 	

library(here) ## assumes that getwd() will return root covid-19-spatial-sim folder.

## Import R source files
invisible(sapply(list.files(here("Rscripts/SourceFiles"), full.names = TRUE), function(x) source(x)))



NUM_AGE_GROUPS 	= 17
CDF_RES			= 20
options("scipen" = 13) #### set high penalty for scientific display. (so e.g. 10000 is not outputted as 1e+05) 
OutputDir  = here("data/param_files") ## change as appropriate



PreParamList = MakePreParamList()







#### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== 
#### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== 
#### ==== Examples

#### Write non-US, non-Canada pre-parameter file
PreParamList = MakePreParamList()
WriteParamList(PreParamList, OutputDir = OutputDir, OutputFileName = "pre_R0_2.0")


#### Parameter files corresponding to following set of interventions:

ArbitrarilyLargeNonStartTime = 100000 ## interventions turned off by specifying a start time greater than the simulation duration, i.e. PreParamList[["Sampling time"]]

# i) 	All interventions
MakeAndWriteParamList(OutputDir = OutputDir, OutputFileName = "p_PC_CI_HQ_SD"	)

# ii) 	No Interventions (sufficient to set start times to be larger than simulation time for interventions not implemented).
MakeAndWriteParamList(OutputDir = OutputDir, OutputFileName = "p_NoInt", 
		CaseIsolationTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		HQuarantineTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		PlaceCloseTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		SocDistTimeStartBase 				= ArbitrarilyLargeNonStartTime) 

# iii) 	Case Isolation Only (sufficient to set start times to be larger than simulation time for interventions not implemented).
MakeAndWriteParamList(OutputDir = OutputDir, OutputFileName = "p_CI", 
		HQuarantineTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		PlaceCloseTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		SocDistTimeStartBase 				= ArbitrarilyLargeNonStartTime)

# iv) 	Case Isolation + Home Quarantine (sufficient to set start times to be larger than simulation time for interventions not implemented).
MakeAndWriteParamList(OutputDir = OutputDir, OutputFileName = "p_CI_HQ", 
		PlaceCloseTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		SocDistTimeStartBase 				= ArbitrarilyLargeNonStartTime)

# v) 	Place Closure only (sufficient to set start times to be larger than simulation time for interventions not implemented).
MakeAndWriteParamList(OutputDir = OutputDir, OutputFileName = "p_PC", 
		CaseIsolationTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		HQuarantineTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		SocDistTimeStartBase 				= ArbitrarilyLargeNonStartTime) 

# vi) 	Social Distancing Only (sufficient to set start times to be larger than simulation time for interventions not implemented).
MakeAndWriteParamList(OutputDir = OutputDir, OutputFileName = "p_SD", 
		CaseIsolationTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		HQuarantineTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		PlaceCloseTimeStartBase 			= ArbitrarilyLargeNonStartTime) 


#### Write parameter files with varying "efficacies" over time. 
MakeAndWriteParamList(OutputDir = OutputDir, OutputFileName = "p_PC_CI_HQ_SDEase_DCT", 
		
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



