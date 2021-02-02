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
#### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== 
#### ==== Make severity progression ICDFs/quantiles (ICDF = Inverse Cumulative Distribution Function)

#### Transition (state)							#### ICDF name in Parameter File	#### ICDF name in R/Cpp code	#### Mean name in Parameter File	#### Mean name in R/Cpp code	#### Notes

# 1.  Latent -> Infectious (Asymp/mild/ILI)		[Latent period inverse CDF]			latent_icdf						[Latent period]						LatentPeriod					Using State E (Latent infection) from KW.
# 2.  Asymp -> Recovery 						[Infectious period inverse CDF]		infectious_icdf					[Infectious period]					InfectiousPeriod				Using State I_A^i (Asymptomatic infection) from KW. In CovidSim, asymptomatics (or everyone if not doing severity) just use ICDF for Infectious Period, so this is Infectious Period ICDF
# 3.  Mild -> Recovery 							[MildToRecovery_icdf]				MildToRecovery_icdf				[Mean_MildToRecovery]				Mean_MildToRecovery				MildToRecovery = ILIToRecovery = ILIToSARI. Using State I_C^i (Symptomatic infection) from KW
# 4.  ILI -> Recovery							[ILIToRecovery_icdf]				ILIToRecovery_icdf				[Mean_ILIToRecovery]				Mean_ILIToRecovery				MildToRecovery = ILIToRecovery = ILIToSARI. Using State I_C^i (Symptomatic infection) from KW
# 5.  ILI -> SARI								[ILIToSARI_icdf]					ILIToSARI_icdf					[Mean_ILIToSARI]					Mean_ILIToSARI					MildToRecovery = ILIToRecovery = ILIToSARI. Using State I_C^i (Symptomatic infection) from KW
# 6.  ILI -> Death								[ILIToDeath_icdf]					ILIToDeath_icdf					[Mean_ILIToDeath]					Mean_ILIToDeath					Sum of I_C^i -> G_D -> Death
# 7.  SARI -> Critical							[SARIToCritical_icdf]				SARIToCritical_icdf				[Mean_SARIToCritical]				Mean_SARIToCritical			 	Using State ICU_pre^i (Triage to ICU) from KW, they separate general ward (what we call SARI) into pre-ICU, hospitalised leading to death, hospitalised leading to recovery
# 8.  SARI -> Death								[SARIToDeath_icdf]					SARIToDeath_icdf                [Mean_SARIToDeath]					Mean_SARIToDeath				Using State H_D (hospitalised on general ward leading to death) from KW.
# 9.  SARI -> Recovery							[SARIToRecovery_icdf] 				SARIToRecovery_icdf             [Mean_SARIToRecovery]				Mean_SARIToRecovery
# 10. Critical -> Death							[CriticalToDeath_icdf] 				CriticalToDeath_icdf            [Mean_CriticalToDeath]				Mean_CriticalToDeath			Doesn't exist in KW.
# 11. Critical -> Stepdown (to death)			[CriticalToStepdownDeath_icdf]		CriticalToStepdownDeath_icdf	[Mean_CriticalToStepdownDeath]		Mean_CriticalToStepdownDeath	Using State ICU_W_D^i (Hospitalised in ICU leading to death) from KW. Doesn't (yet) exist in CovidSim 
# 12. Critical -> Stepdown (to discharge)       [CriticalToCritRecov_icdf] 			CriticalToCritRecov_icdf    	[Mean_CriticalToCritRecov] 			Mean_CriticalToCritRecov   
# 13. Stepdown -> Death                         																									                                    Doesn't (yet) exist in CovidSim 
# 14. Stepdown -> Recovery                      [CritRecovToRecov_icdf]				CritRecovToRecov_icdf       	[Mean_CritRecovToRecov]				Mean_CritRecovToRecov     

N_Samples = 1000000
NUM_AGE_GROUPS = 17
QUANTILES = c(seq(0, 0.95, 0.05), 0.99)

# 1. Latent -> Infectious
LatentPeriod 		= 2/0.44
latent_icdf 		= qgamma(QUANTILES, shape = 2, rate = 2)

# 2. Asymp -> Recovery (Infectious Period - see Notes)
InfectiousPeriod 	= 1/0.48
infectious_icdf 	= qgamma(QUANTILES, shape = 1, rate = 1) 

# 3. Mild -> Recovery 
Mean_MildToRecovery	= rep(1/0.25, NUM_AGE_GROUPS)
MildToRecovery_icdf	= qgamma(QUANTILES, shape = 1, rate = 1) 

# 4. ILI -> Recovery
Mean_ILIToRecovery = Mean_MildToRecovery
ILIToRecovery_icdf = MildToRecovery_icdf

# 5. ILI -> SARI
Mean_ILIToSARI = Mean_MildToRecovery
ILIToSARI_icdf = MildToRecovery_icdf 

# 6. ILI -> Death. Sum of I_C^i -> G_D -> Death. Sojourn time from I_C^i -> G_D = Mean_ILIToSARI = Mean_ILIToRecovery = Mean_MildToRecovery
Mean_ILIToDeath = Mean_ILIToRecovery + (2/0.4)
# create samples
ILI_to_Death_Samples 	= rgamma(n = N_Samples, shape = 1, rate = 0.25) + rgamma(n = N_Samples, shape = 2, rate = 0.4)
Mean_ILIToDeath_Samples = mean(ILI_to_Death_Samples) ## check this equals Mean_ILIToDeath (roughly)
ILIToDeath_icdf 		= quantile(ILI_to_Death_Samples, probs = QUANTILES, names = FALSE)
Mean_ILIToDeath_Samples - Mean_ILIToDeath ## should be very small

# 7.  SARI -> Critical
Mean_SARIToCritical = rep(1/0.40, NUM_AGE_GROUPS)
SARIToCritical_icdf = qgamma(QUANTILES, shape = 1, rate = 1) 

# 8.  SARI -> Death
Mean_SARIToDeath = rep(2/0.19, NUM_AGE_GROUPS)
SARIToDeath_icdf = qgamma(QUANTILES, shape = 2, rate = 2) 

# 9.  SARI -> Recovery
Mean_SARIToRecovery = rep(1/0.09, NUM_AGE_GROUPS)
SARIToRecovery_icdf = qgamma(QUANTILES, shape = 1, rate = 1) 

# 10. Critical -> Death. Either CovidSim keeps this as one transition, in which case it is a composite of two transitions from KW. Or we model it as two transitions with a stepdown component. Then questions is whether we want the hassle of two step-down compartments. 
Mean_CriticalToDeath 			= rep(1/0.14 + 1/0.12, NUM_AGE_GROUPS)
CriticalToDeath_Samples 		= rgamma(n = N_Samples, shape = 1, rate = 0.14) + rgamma(n = N_Samples, shape = 1, rate = 0.12)
CriticalToDeath_icdf 			= quantile(CriticalToDeath_Samples, probs = QUANTILES, names = FALSE)
Mean_CriticalToDeath_Samples	= mean(CriticalToDeath_Samples)
Mean_CriticalToDeath_Samples - Mean_CriticalToDeath ## should be very small

# 11. Critical -> Stepdown (to death)
Mean_CriticalToStepdownDeath = rep(1/0.14, NUM_AGE_GROUPS)
CriticalToStepdownDeath_icdf = qgamma(QUANTILES, shape = 1, rate = 1) 

# 12. Critical -> Stepdown (to discharge) 
Mean_CriticalToCritRecov = rep(1/0.06, NUM_AGE_GROUPS)
CriticalToCritRecov_icdf = qgamma(QUANTILES, shape = 1, rate = 1) 

# 13. Stepdown -> Death - doesn't exist in CovidSim - 
Mean_StepdownToDeath = rep(1/0.12, NUM_AGE_GROUPS)
StepdownToDeath_icdf = qgamma(QUANTILES, shape = 1, rate = 1)

# 14. Stepdown -> Recovery
Mean_CritRecovToRecov = rep(2/0.16, NUM_AGE_GROUPS)
CritRecovToRecov_icdf = qgamma(QUANTILES, shape = 2, rate = 2) 

#### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== 
#### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # Add in WAIFW matrix using POLYMOD/socialmixr
library(socialmixr)

Ages 			= seq(0,75,5)
POLYMOD_UK 		= contact_matrix(polymod, countries = "United Kingdom", age.limits = Ages)
WAIFW_Matrix 	= POLYMOD_UK$matrix
## contact_matrix only returns as high as 75+. CovidSim wants 75-80 and 80+.
## For now assume that [75,80) = 80+ = 75+. Check this
WAIFW_Matrix = rbind(Matrix, Matrix[nrow(Matrix),]) # add row to matrix
WAIFW_Matrix = cbind(Matrix, Matrix[,ncol(Matrix)]) # add col to matrix



#### Write non-US, non-Canada pre-parameter file
PreParamList = MakePreParamList(
		
		WAIFW_Matrix = WAIFW_Matrix,
		
		# 1. Latent -> Infectious
		LatentPeriod 		= LatentPeriod 		, 
		latent_icdf 		= latent_icdf 		, 
		
		# 2. Asymp -> Recovery (Infectious Period - see Notes)
		InfectiousPeriod 	= InfectiousPeriod 	, 
		infectious_icdf 	= infectious_icdf 	, 

		# 3. Mild -> Recovery 
		Mean_MildToRecovery	= Mean_MildToRecovery	,
		MildToRecovery_icdf	= MildToRecovery_icdf	, 

		# 4. ILI -> Recovery                                           
		Mean_ILIToRecovery = Mean_ILIToRecovery ,
		ILIToRecovery_icdf = ILIToRecovery_icdf ,
		
		# 5. ILI -> SARI
		Mean_ILIToSARI = Mean_ILIToSARI ,
		ILIToSARI_icdf = ILIToSARI_icdf ,
		
		# 6. ILI -> Death. Sum of I_C^i -> G_D -> Death. Sojourn time from I_C^i -> G_D = Mean_ILIToSARI = Mean_ILIToRecovery = Mean_MildToRecovery
		Mean_ILIToDeath = Mean_ILIToDeath,
		ILIToDeath_icdf = ILIToDeath_icdf,
		
		# 7.  SARI -> Critical		
		Mean_SARIToCritical = Mean_SARIToCritical,
		SARIToCritical_icdf = SARIToCritical_icdf,
		
		# 8.  SARI -> Death
		Mean_SARIToDeath = Mean_SARIToDeath,
		SARIToDeath_icdf = SARIToDeath_icdf,

		# 9.  SARI -> Recovery
		Mean_SARIToRecovery = Mean_SARIToRecovery,
		SARIToRecovery_icdf = SARIToRecovery_icdf,

		# 10. Critical -> Death. Either CovidSim keeps this as one transition, in which case it is a composite of two transitions from KW. Or we model it as two transitions with a stepdown component. Then questions is whether we want the hassle of two step-down compartments. 
		Mean_CriticalToDeath = Mean_CriticalToDeath,
		CriticalToDeath_icdf = CriticalToDeath_icdf,
		
		# 11. Critical -> Stepdown (to death)
		Mean_CriticalToStepdownDeath = Mean_CriticalToStepdownDeath,
		CriticalToStepdownDeath_icdf = CriticalToStepdownDeath_icdf, 

		# 12. Critical -> Stepdown (to discharge) 
		Mean_CriticalToCritRecov = Mean_CriticalToCritRecov,
		CriticalToCritRecov_icdf = CriticalToCritRecov_icdf,
		
		# 13. Stepdown -> Death - doesn't exist in CovidSim - 
		Mean_StepdownToDeath = Mean_StepdownToDeath,
		StepdownToDeath_icdf = StepdownToDeath_icdf,
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



