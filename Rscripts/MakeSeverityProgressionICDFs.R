#### ==== #### ==== #### ==== #### ==== 
#### ==== Author: dlaydon
#### Script calculates delay distributions and means for severity progressions for CovidSim.
#### Script uses fitted parameters values from Knock, Whittles et al (abbreviated to KW henceforth, URL = https://www.medrxiv.org/content/10.1101/2021.01.11.21249564v1.full.pdf)
#### Fitted Gamma/Erlang shape and rate parameters found in KW Table S2, and labelled in Fig S2 of KW. 


library(here) ## assumes that getwd() will return root "covid-19-spatial-sim" folder.
OutputDir  = here("data/param_files/") ## change as appropriate

## Import Erlang/Gamma rate and shape parameters for progression / sojourn times
ErlangInfo = read.csv(here("data/param_files/support_progression.csv"))
ProgParams = ErlangInfo$value
names(ProgParams) = as.character(ErlangInfo$parameter)

# Import transition probabilities
TransitionProbabilites = read.csv(here("data/param_files/support_severity.csv"))

# Extract transition probabilities to match naming in CovidSim.


#### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== 
#### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== 
#### ==== Make severity progression ICDFs/quantiles (ICDF = Inverse Cumulative Distribution Function)

#		TRANSITION									ICDF name in pre-parameter file		ICDF name in Cpp code			Mean name in pre-parameter file		Mean name in Cpp code			NOTES

# 1.  	Latent -> Infectious (Asymp/mild/ILI)		[Latent period inverse CDF]			latent_icdf						[Latent period]						LatentPeriod					Using State E (Latent infection) from KW.
# 2.  	Asymp -> Recovery 							[Infectious period inverse CDF]		infectious_icdf					[Infectious period]					InfectiousPeriod				Using State I_A^i (Asymptomatic infection) from KW. In CovidSim, asymptomatics (or everyone if not doing severity) just use ICDF for Infectious Period, so this is Infectious Period ICDF
# 3.  	Mild -> Recovery 							[MildToRecovery_icdf]				MildToRecovery_icdf				[Mean_MildToRecovery]				Mean_MildToRecovery				MildToRecovery = ILIToRecovery = ILIToSARI. Using State I_C^i (Symptomatic infection) from KW
# 4.  	ILI -> Recovery								[ILIToRecovery_icdf]				ILIToRecovery_icdf				[Mean_ILIToRecovery]				Mean_ILIToRecovery				MildToRecovery = ILIToRecovery = ILIToSARI. Using State I_C^i (Symptomatic infection) from KW
# 5.  	ILI -> SARI									[ILIToSARI_icdf]					ILIToSARI_icdf					[Mean_ILIToSARI]					Mean_ILIToSARI					MildToRecovery = ILIToRecovery = ILIToSARI. Using State I_C^i (Symptomatic infection) from KW
# 6.  	ILI -> Death								[ILIToDeath_icdf]					ILIToDeath_icdf					[Mean_ILIToDeath]					Mean_ILIToDeath					Sum of I_C^i -> G_D -> Death
# 7.  	SARI -> Critical							[SARIToCritical_icdf]				SARIToCritical_icdf				[Mean_SARIToCritical]				Mean_SARIToCritical			 	Using State ICU_pre^i (Triage to ICU) from KW, they separate general ward (what we call SARI) into pre-ICU, hospitalised leading to death, hospitalised leading to recovery
# 8.  	SARI -> Death								[SARIToDeath_icdf]					SARIToDeath_icdf                [Mean_SARIToDeath]					Mean_SARIToDeath				Using State H_D^i (hospitalised on general ward leading to death) from KW.
# 9.  	SARI -> Recovery							[SARIToRecovery_icdf] 				SARIToRecovery_icdf             [Mean_SARIToRecovery]				Mean_SARIToRecovery				Using State H_R^i
# 10. 	Critical -> Death							[CriticalToDeath_icdf] 				CriticalToDeath_icdf            [Mean_CriticalToDeath]				Mean_CriticalToDeath			Using State ICU_D^i.
# 11. 	Critical -> Stepdown (to death)				[CriticalToStepdownDeath_icdf]		CriticalToStepdownDeath_icdf	[Mean_CriticalToStepdownDeath]		Mean_CriticalToStepdownDeath	Using State ICU_W_D^i (Hospitalised in ICU leading to death) from KW. Doesn't (yet) exist in CovidSim 
# 12. 	Critical -> Stepdown (to discharge)  	    [CriticalToCritRecov_icdf] 			CriticalToCritRecov_icdf    	[Mean_CriticalToCritRecov] 			Mean_CriticalToCritRecov   
# 13. 	Stepdown -> Death                         																									                                    	Doesn't (yet) exist in CovidSim 
# 14. 	Stepdown -> Recovery                    	[CritRecovToRecov_icdf]				CritRecovToRecov_icdf       	[Mean_CritRecovToRecov]				Mean_CritRecovToRecov     

### Note all Erlang/Gamma distributions calculated below chosen to have shape = rate so they have a mean value of 1. This enables age-specific means to be inputted if need be without fitting separate shape parameters.


N_Samples 		= 1000000
NUM_AGE_GROUPS 	= 17
QUANTILES 		= c(seq(0, 0.95, 0.05), 0.99)

# 1. Latent -> Infectious. Using State E (Latent infection) from KW.
LatentPeriod 		= 2/0.44
latent_icdf 		= qgamma(QUANTILES, shape = 2, rate = 2)

# 2. Asymp -> Recovery. Using State I_A^i (Asymptomatic infection) from KW. (See notes above)
InfectiousPeriod 	= 1/0.48
infectious_icdf 	= qgamma(QUANTILES, shape = 1, rate = 1) 

# 3. Mild -> Recovery. Using State I_C^i (Symptomatic infection) from KW
Mean_MildToRecovery	= rep(1/0.25, NUM_AGE_GROUPS)
MildToRecovery_icdf	= qgamma(QUANTILES, shape = 1, rate = 1) 

# 4. ILI -> Recovery. Using State I_C^i (Symptomatic infection) from KW
Mean_ILIToRecovery = rep(1/0.25, NUM_AGE_GROUPS) 	# = Mean_MildToRecovery
ILIToRecovery_icdf = MildToRecovery_icdf

# 5. ILI -> SARI. Using State I_C^i (Symptomatic infection) from KW
Mean_ILIToSARI = rep(1/0.25, NUM_AGE_GROUPS) 		# = Mean_MildToRecovery
ILIToSARI_icdf = MildToRecovery_icdf 

# 6. ILI -> Death. Sum of I_C^i -> G_D -> Death. Sojourn time from I_C^i -> G_D = Mean_ILIToSARI = Mean_ILIToRecovery = Mean_MildToRecovery
Mean_ILIToDeath = Mean_ILIToRecovery + (2/0.4)
# create samples
ILI_to_Death_Samples 	= rgamma(n = N_Samples, shape = 1, rate = 0.25) + rgamma(n = N_Samples, shape = 2, rate = 0.4)
Mean_ILIToDeath_Samples = mean(ILI_to_Death_Samples) ## check this equals Mean_ILIToDeath (roughly)
ILIToDeath_icdf 		= quantile(ILI_to_Death_Samples, probs = QUANTILES, names = FALSE)
Mean_ILIToDeath_Samples - Mean_ILIToDeath ## should be very small

# 7.  SARI -> Critical. Using State ICU_pre^i (Triage to ICU)
#Mean_SARIToCritical = rep(1/0.40, NUM_AGE_GROUPS)
Mean_SARIToCritical = rep(ProgParams["k_ICU_pre"] / ProgParams["gamma_ICU_pre"], NUM_AGE_GROUPS)
SARIToCritical_icdf = qgamma(QUANTILES, shape = ProgParams["k_ICU_pre"], rate = ProgParams["k_ICU_pre"]) 

# 8.  SARI -> Death. Using State H_D^i
#Mean_SARIToDeath = rep(2/0.19, NUM_AGE_GROUPS)
Mean_SARIToDeath = rep(ProgParams["k_H_D"] / ProgParams["gamma_H_D"], NUM_AGE_GROUPS)
SARIToDeath_icdf = qgamma(QUANTILES, shape = ProgParams["k_H_D"], rate = ProgParams["k_H_D"]) 

# 9.  SARI -> Recovery. Using State H_R^i
#Mean_SARIToRecovery = rep(1/0.09, NUM_AGE_GROUPS)
Mean_SARIToRecovery = rep(ProgParams["k_H_R"] / ProgParams["gamma_H_R"], NUM_AGE_GROUPS)
SARIToRecovery_icdf = qgamma(QUANTILES, shape = ProgParams["k_H_R"], rate = ProgParams["k_H_R"]) 

# 10. Critical -> Death. Using State ICU_D^i.
Mean_CriticalToDeath 		= rep(ProgParams["k_ICU_D"] / ProgParams["gamma_ICU_D"], NUM_AGE_GROUPS)
CriticalToDeath_icdf 		= qgamma(QUANTILES, shape = ProgParams["k_ICU_D"], rate = ProgParams["k_ICU_D"]) 

# 11. Critical -> Stepdown (to death). Using State ICU_W_D^i
Mean_CriticalToStepdownDeath = rep(ProgParams["k_ICU_W_D"] / ProgParams["gamma_ICU_W_D"], NUM_AGE_GROUPS)         
CriticalToStepdownDeath_icdf = qgamma(QUANTILES, shape = ProgParams["k_ICU_W_D"], rate = ProgParams["k_ICU_W_D"]) 

# 12. Critical -> Stepdown (to discharge). Using State ICU_W_R^i
Mean_CriticalToCritRecov = rep(ProgParams["k_ICU_W_R"] / ProgParams["gamma_ICU_W_R"], NUM_AGE_GROUPS)            
CriticalToCritRecov_icdf = qgamma(QUANTILES, shape = ProgParams["k_ICU_W_R"], rate = ProgParams["k_ICU_W_R"])    

# 13. Stepdown -> Death. Using State W_D^i. Need to put into CovidSim.
Mean_StepdownToDeath = rep(ProgParams["k_W_D"] / ProgParams["gamma_W_D"], NUM_AGE_GROUPS)        
StepdownToDeath_icdf = qgamma(QUANTILES, shape = ProgParams["k_W_D"], rate = ProgParams["k_W_D"])

# 14. Stepdown -> Recovery. Using State W_R^i.
Mean_CritRecovToRecov = rep(ProgParams["k_W_R"] / ProgParams["gamma_W_R"], NUM_AGE_GROUPS)        
CritRecovToRecov_icdf = qgamma(QUANTILES, shape = ProgParams["k_W_R"], rate = ProgParams["k_W_R"])




## Save output
AgeSpecificMeans = data.frame(LatentPeriod, InfectiousPeriod, Mean_MildToRecovery, 
		Mean_ILIToRecovery, Mean_ILIToSARI, Mean_ILIToSARI, Mean_ILIToDeath, 
		Mean_SARIToCritical, Mean_SARIToDeath, Mean_SARIToRecovery, 
		Mean_CriticalToDeath, Mean_CriticalToStepdownDeath, Mean_CriticalToCritRecov, 
		Mean_StepdownToDeath, Mean_CritRecovToRecov) # note not actually making means age-specific anymore, but option there if need be. 
rownames(AgeSpecificMeans) = AgeGroupNames
saveRDS(AgeSpecificMeans, file = file.path(OutputDir, "AS_SeverityProgMeans.rds"))

## Save output
ICDFs = data.frame(Quantile = seq(0,1, 0.05), latent_icdf, infectious_icdf, MildToRecovery_icdf, 
		ILIToRecovery_icdf, ILIToSARI_icdf, ILIToSARI_icdf, ILIToDeath_icdf, 
		SARIToCritical_icdf, SARIToDeath_icdf, SARIToRecovery_icdf, 
		CriticalToDeath_icdf, CriticalToStepdownDeath_icdf, CriticalToCritRecov_icdf, 
		StepdownToDeath_icdf, CritRecovToRecov_icdf) # note not actually making means age-specific anymore, but option there if need be. 
saveRDS(ICDFs, file = file.path(OutputDir, "SeverityICDFs.rds"))



#AgeSpecificSeverityScalings = matrix(c(
#				
#				0.039, 0.243, 0.039, 0.282, 0.091, 0, 
#				0.001, 0.289, 0.037, 0.286, 0.083, 0, 
#				0.006, 0.338, 0.035, 0.291, 0.077, 0, 
#				0.009, 0.389, 0.035, 0.299, 0.074, 0, 
#				0.026, 0.443, 0.036, 0.310, 0.074, 0, 
#				0.040, 0.503, 0.039, 0.328, 0.076, 0, 
#				0.042, 0.570, 0.045, 0.353, 0.080, 0, 
#				0.045, 0.653, 0.055, 0.390, 0.086, 0, 
#				0.050, 0.756, 0.074, 0.446, 0.093, 0, 
#				0.074, 0.866, 0.107, 0.520, 0.102, 0, 
#				0.138, 0.954, 0.157, 0.604, 0.117, 0, 
#				0.198, 1.000, 0.238, 0.705, 0.148, 0, 
#				0.247, 0.972, 0.353, 0.806, 0.211, 0, 
#				0.414, 0.854, 0.502, 0.899, 0.332, 0, 
#				0.638, 0.645, 0.675, 0.969, 0.526, 0, 
#				1.000, 0.402, 0.832, 1.000, 0.753, 0, 
#				0.873, 0.107, 1.000, 0.918, 1.000, 0), ncol = 6, byrow = TRUE)

AgeGroupNames = c("[0,5)", "[5,10)", "[10,15)", "[15,20)", "[20,25)", "[25,30)", "[30,35)", "[35,40)", "[40,45)", 
		"[45,50)", "[50,55)", "[55,60)", "[60,65)", "[65,70)", "[70,75)", "[75,80)", "80+") 



