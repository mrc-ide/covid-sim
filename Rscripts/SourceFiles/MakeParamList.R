MakeParamList = function(
		
		#### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== 
		#### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== 
		#### ==== 			NOTE: Defaults have all interventions enabled. 
		
		#### Interventions are turned on when BOTH start time has passed (and start time + duration has not), AND when trigger threshold exceeded. See https://github.com/mrc-ide/covid-sim/blob/master/docs/intervention-description.md for further discussion.
		#### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== 
		
		NumAdUnits = 3, DoInterventionDelaysByAdUnit = 0, VaryEfficaciesOverTime = 0, 
		
		# = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = 
		# = # = # = # = # = # = # = # = # = # = 	CASE ISOLATION PARAMETERS
		
		##### NOTE: Default is to have case isolation ENABLED
		
		CaseIsolationTimeStartBase = 6,		# a value larger than PreParamList[["Sampling time"]] indicates that case isolation will not start, regardless of trigger threshold(s) CaseIsolation_CellIncThresh / CaseIsolation_CellIncThresh. Otherwise intervention will start once threshold exceeded. 
		CaseIsolationDelay = 1, 			# delay between case detection and case isolation
		CaseIsolationDuration = 7, 			# CaseIsolationDuration = number of days a case will self-isolate. Different from CaseIsolationPolicyDuration = number of days policy is in effect.
		CaseIsolationPolicyDuration = 364, 	# CaseIsolationDuration = number of days a case will self-isolate. Different from CaseIsolationPolicyDuration = number of days policy is in effect.
		CI_Delay_byAdUnit 			= rep(CaseIsolationDelay			, NumAdUnits),  # delay between case detection and case isolation by admin unit
		CI_PolicyDuration_byAdUnit 	= rep(CaseIsolationPolicyDuration	, NumAdUnits), 	# duration of case isolation policy by admin unit.
		
		CaseIsolation_CellIncThresh = 0, CaseIsolationProp = 0.9, CaseIsolationEffectiveness = 0.25, CaseIsolationHouseEffectiveness = 0.5, 
		
		#### Note: numbers here must match Num_CI_ChangeTimes; any times listed in CI_ChangeTimes that are before CaseIsolationTimeStartBase and after CaseIsolationPolicyDuration / CI_PolicyDuration_byAdUnit are irrelevant.
		Num_CI_ChangeTimes = 1, CI_ChangeTimes = c(CaseIsolationTimeStartBase, rep(1000000, Num_CI_ChangeTimes - 1)), ### by default, initialize first time to non-variable time start, then all other times to be arbitrarily large. 
		CI_SpatialAndPlaceEffects_OverTime 		= rep(CaseIsolationEffectiveness			, Num_CI_ChangeTimes), 
		CI_HouseholdEffects_OverTime 			= rep(CaseIsolationHouseEffectiveness		, Num_CI_ChangeTimes),
		CI_Prop_OverTime						= rep(CaseIsolationProp						, Num_CI_ChangeTimes),
		CI_CellIncThresh_OverTime				= rep(CaseIsolation_CellIncThresh			, Num_CI_ChangeTimes),
		
		
		# = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = 
		# = # = # = # = # = # = # = # = # = # = 	HOUSEHOLD QUARANTINE PARAMETERS
		
		##### NOTE: Default is to have household quarantine ENABLED
		
		HQuarantineTimeStartBase = 6, ### a value larger than PreParamList[["Sampling time"]] indicates that case isolation will not start, regardless of trigger threshold(s) CaseIsolation_CellIncThresh / CaseIsolation_CellIncThresh. Otherwise intervention will start once threshold exceeded.  
		HQuarantineDelay = 0, HQuarantinePolicyDuration = 364,
		HQuarantineDelay_byAdUnit 			= rep(HQuarantineDelay			, NumAdUnits), 
		HQuarantinePolicyDuration_byAdUnit 	= rep(HQuarantinePolicyDuration	, NumAdUnits),
		
		HHQuar_CellIncThresh = 0, HQuarantineHouseDuration = 14, 
		HQuarantineHouseEffect = 1.5, HQuarantinePlaceEffect = c(0.25, 0.25, 0.25, 0.25), HQuarantineSpatialEffect = 0.25,
		HQuarantinePropHouseCompliant = 0.75, HQuarantinePropIndivCompliant = 1.0,
		
		#### Note: numbers here must match Num_HQ_ChangeTimes; any times listed in HQ_ChangeTimes that are before HQuarantineTimeStartBase and after HQuarantinePolicyDuration / HQuarantinePolicyDuration_byAdUnit are irrelevant.
		Num_HQ_ChangeTimes = 1, HQ_ChangeTimes = c(HQuarantineTimeStartBase, rep(1000000, Num_HQ_ChangeTimes - 1)) , # by default, initialize first time to non-variable time start, then all other times to be arbitrarily large. 
		HQ_HouseholdEffects_OverTime 	= rep(HQuarantineHouseEffect		, Num_HQ_ChangeTimes),
		HQ_SpatialEffects_OverTime 		= rep(HQuarantineSpatialEffect		, Num_HQ_ChangeTimes),
		HQ_PlaceEffects_OverTime		= matrix(rep(HQuarantinePlaceEffect	, Num_HQ_ChangeTimes), nrow = Num_HQ_ChangeTimes, byrow = TRUE),
		
		HQ_Household_PropComply_OverTime 	= rep(HQuarantinePropHouseCompliant	, Num_HQ_ChangeTimes), # household compliance
		HQ_Individual_PropComply_OverTime 	= rep(HQuarantinePropIndivCompliant	, Num_HQ_ChangeTimes), # individual compliance
		
		HQ_CellIncThresh_OverTime 			= rep(HHQuar_CellIncThresh			, Num_HQ_ChangeTimes), # thresholds
		
		
		# = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = 
		# = # = # = # = # = # = # = # = # = # = 	PLACE CLOSURE PARAMETERS
		
		##### NOTE: Default is to have place closure ENABLED
		
		PlaceCloseTimeStartBase = 7, PlaceCloseTimeStartBase2 = 100000, PlaceCloseDelayMean = 1, PlaceCloseDurationBase = 364, PlaceCloseRadius = 1,  
		PlaceCloseByAdminUnit = 0, PlaceCloseAdminUnitDivisor = 1, PlaceCloseAdunitPlaceTypes = c(1, 1, 1, 0), PlaceCloseCasePropThresh = 1, PlaceCloseAdunitPropThresh = 1,
		PlaceCloseDelayMean_byAdUnit 	= rep(PlaceCloseDelayMean		, NumAdUnits),
		PlaceCloseDuration_byAdUnit 	= rep(PlaceCloseDurationBase	, NumAdUnits),
		PlaceCloseHouseholdRelContact 		= 1.5, 
		PlaceCloseSpatialRelContact 		= 1.25, 
		PropPlacesRemainingOpen_byPlaceType	= c(0, 0, 0.25, 1), 
		PlaceClosePropAttending = c(0.05, 0.05, 0.05, 0.05), ### Partial closure: for places that are closed, what proportion of their members do attend? Set to zero if no partiall admittance. 
		PlaceCloseIncTrig = 0, PlaceCloseFracIncTrig = 0, PlaceCloseCellIncThresh = 0, PlaceCloseCellIncStopThresh = 0,
		
		#### Note: numbers here must match Num_PC_ChangeTimes; any times listed in PC_ChangeTimes that are before PlaceCloseTimeStartBase and after PlaceCloseDurationBase / PlaceCloseDuration_byAdUnit are irrelevant.
		Num_PC_ChangeTimes = 1,	PC_ChangeTimes = c(PlaceCloseTimeStartBase, rep(1000000, Num_PC_ChangeTimes - 1))		, ### by default, initialize first time to non-variable time start, then all other times to be arbitrarily large. 
		
		PC_IncThresh_OverTime 			= rep(PlaceCloseIncTrig				, Num_PC_ChangeTimes)	,
		PC_FracIncThresh_OverTime 		= rep(PlaceCloseFracIncTrig			, Num_PC_ChangeTimes) 	, 
		PC_CellIncThresh_OverTime 		= rep(PlaceCloseCellIncThresh		, Num_PC_ChangeTimes) 	, 
		PC_HouseholdEffects_OverTime 	= rep(PlaceCloseHouseholdRelContact	, Num_PC_ChangeTimes)	,
		PC_SpatialEffects_OverTime 		= rep(PlaceCloseSpatialRelContact	, Num_PC_ChangeTimes)	,
		
		PC_PlaceEffects_OverTime 		= matrix(rep(PropPlacesRemainingOpen_byPlaceType, Num_PC_ChangeTimes), nrow = Num_PC_ChangeTimes, byrow = TRUE),
		PC_Durs_OverTime				= rep(PlaceCloseDurationBase		, Num_PC_ChangeTimes)	,
		PC_PropAttending_OverTime		= matrix(rep(PlaceClosePropAttending			, Num_PC_ChangeTimes), nrow = Num_PC_ChangeTimes, byrow = TRUE),
		
		# = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = 
		# = # = # = # = # = # = # = # = # = # = 	SOCIAL DISTANCING PARAMETERS
		
		##### NOTE: Default is to have social distancing ENABLED
		
		SocDistTimeStartBase = 6, SocDistDuration = 364, 
		SocDistDelay_byAdUnit 		= rep(0, NumAdUnits),
		SocDistDuration_byAdUnit 	= rep(SocDistDuration, NumAdUnits),
		
		SocDistCellIncThresh = 0, SocDistCellIncStopThresh = 0,
		
		SocDistPlaceEffect 		= c(1, 1, 0.5, 0.5)	, EnhancedSocDistPlaceEffect 		= c(0.25, 0.25, 0.25, 0.25), 
		SocDistHouseholdEffect 	= 1.25				, EnhancedSocDistHouseholdEffect 	= 1, 
		SocDistSpatialEffect 	= 0.1				, EnhancedSocDistSpatialEffect 		= 0.25, 
		
		EnhancedSocDistProportionCompliantNonAge = 0, EnhancedSocDistProportionCompliantByAge = rep(0, NUM_AGE_GROUPS),
		
		#### Note: numbers here must match Num_SD_ChangeTimes; any times listed in SD_ChangeTimes that are before SocDistTimeStartBase and after SocDistDuration / SocDistDuration_byAdUnit are irrelevant.
		Num_SD_ChangeTimes = 1, SD_ChangeTimes = c(SocDistTimeStartBase, rep(1000000, Num_SD_ChangeTimes - 1)), ### by default, initialize first time to non-variable time start, then all other times to be arbitrarily large. 
		
		SD_PlaceEffects_OverTime 				= matrix(rep(SocDistPlaceEffect			, Num_SD_ChangeTimes), nrow = Num_SD_ChangeTimes, byrow = TRUE),
		Enhanced_SD_PlaceEffects_OverTime 		= matrix(rep(EnhancedSocDistPlaceEffect	, Num_SD_ChangeTimes), nrow = Num_SD_ChangeTimes, byrow = TRUE),
		SD_HouseholdEffects_OverTime 			= rep(SocDistHouseholdEffect			, Num_SD_ChangeTimes), 
		Enhanced_SD_HouseholdEffects_OverTime	= rep(EnhancedSocDistHouseholdEffect	, Num_SD_ChangeTimes),
		SD_SpatialEffects_OverTime				= rep(SocDistSpatialEffect				, Num_SD_ChangeTimes), 
		Enhanced_SD_SpatialEffects_OverTime		= rep(EnhancedSocDistSpatialEffect		, Num_SD_ChangeTimes), 
		
		SD_CellIncThresh_OverTime 				= rep(SocDistCellIncThresh				, Num_SD_ChangeTimes), 
		
		# = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = 
		# = # = # = # = # = # = # = # = # = # = 	DIGITAL CONTACT TRACING PARAMETERS
		
		##### NOTE: Default is to have digital contact tracing ENABLED
		
		DoDigitalContactTracing = 1, DigitalContactTracingTimeStartBase = 0, # i.e. will start as soon as trigger threshold reached.
		OutputDigitalContactTracing = 1, # Output the a file showing the number of people under isolation due to contact tracing in each admin unit at any given time. Only relevant if DoDigitalContactTracing == 1.
		ClusterDigitalContactUsers = 0, #by default, don't cluster by location. #If we cluster by household, then we select a proportion of households to be potential app users, if not then we select people over the whole population
		PropPopUsingDigitalContactTracing = 0.6, 
		ProportionSmartphoneUsersByAge = c(0, 0, 0, 0.96, 0.96, 0.96, 0.96, 0.91, 0.91, 0.91, 0.91, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55),  # (Taken from p. 189 here: https://www.ofcom.org.uk/__data/assets/pdf_file/0026/143981/technology-tracker-2019-uk-data-tables.pdf)
		# This mobile phone data would be used specifically for the digital app. If we include actual phone usage data and don't cluster by household, then we adjust the probability of any individual to be a user in order to obtain the specified proportion - however this is limited by a theoretical maximum and if our desired probability exceeds this, we'll get an error.  We don't currently adjust the probabilities if clustering by household - the proportion of households containing users will match the desired proportion but due to age dependent usage, the actual proportion of the population who are users will be lower.
		# NOTE: if we set the proportion of users by age to 1 for all age groups, this means everyone can be tracked and is more similar to standard contact tracing.
		
		ProportionDigitalContactsIsolate = 0.4, DelayFromIndexCaseDetectionToDCTIsolation = 0, LengthDigitalContactIsolation = 14, DigitalContactTracingDelay = 0,
		ScalingFactorSpatialDigitalContacts = 10, ScalingFactorPlaceDigitalContacts = 3, 
		DigitalContactTracing_CellIncThresh = 100, DigitalContactTracingPolicyDuration = 121,
		
		DCT_Duration_byAdUnit 	= rep(DigitalContactTracingPolicyDuration, NumAdUnits),
		DCT_Delay_byAdUnit 		= rep(0, NumAdUnits), ### note there isn't a non-admin unit version of this variable in cpp code yet. 
		DCTIsolateIndexCases 	= 1, 
		
		DelayToTestIndexCase = 0, DelayToTestDCTContacts = 0, 
		
		SensitivityDCT = 1, SpecificityDCT = 1,	FindContactsOfDCTContacts = 0, DoDCTTest = 0, RemoveContactsOfNegativeIndexCase = 0, 
		DCTCaseIsolationEffectiveness = 0.25, DCTCaseIsolationHouseEffectiveness = 0.25, 
		
		Num_DCT_ChangeTimes = 1, DCT_ChangeTimes = c(DigitalContactTracingTimeStartBase, rep(1000000, Num_DCT_ChangeTimes - 1)), ### by default, initialize first time to non-variable time start, then all other times to be arbitrarily large. 
		DCT_SpatialAndPlaceEffects_OverTime = rep(DCTCaseIsolationEffectiveness		, Num_DCT_ChangeTimes), 
		DCT_HouseholdEffects_OverTime 		= rep(DCTCaseIsolationHouseEffectiveness, Num_DCT_ChangeTimes), 
		DCT_Prop_OverTime 					= rep(ProportionDigitalContactsIsolate	, Num_DCT_ChangeTimes)

)
{
	#### Add checks - e.g. number of change times must match number of levels. Add warnings for change times set before and after start and start + duration of policy
	
	#### if change times not provided, default is initialize time to first start time and then arbitrarily large change times (i.e. they won't change). If change times provided, default is to initialize Num_XX_ChangeTimes to equal length XX_ChangeTimes
	if (is.null(CI_ChangeTimes	)) CI_ChangeTimes 	= c(CaseIsolationTimeStartBase			, rep(1000000, Num_CI_ChangeTimes 	- 1)) else Num_CI_ChangeTimes 	= length(CI_ChangeTimes)
	if (is.null(HQ_ChangeTimes	)) HQ_ChangeTimes 	= c(HQuarantineTimeStartBase			, rep(1000000, Num_HQ_ChangeTimes 	- 1)) else Num_HQ_ChangeTimes 	= length(HQ_ChangeTimes)
	if (is.null(PC_ChangeTimes	)) PC_ChangeTimes 	= c(PlaceCloseTimeStartBase				, rep(1000000, Num_PC_ChangeTimes 	- 1)) else Num_PC_ChangeTimes 	= length(PC_ChangeTimes)
	if (is.null(SD_ChangeTimes	)) SD_ChangeTimes 	= c(SocDistTimeStartBase				, rep(1000000, Num_SD_ChangeTimes 	- 1)) else Num_SD_ChangeTimes 	= length(SD_ChangeTimes)
	
	
	ParamList = list()
	
	ParamList[["Include intervention delays by admin unit"]] = DoInterventionDelaysByAdUnit
	ParamList[["Vary efficacies over time"]] = VaryEfficaciesOverTime
	
	# = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = 
	# = # = # = # = # = # = # = # = # = # = 	CASE ISOLATION PARAMETERS
	
	ParamList[["Case isolation start time"]] 					= CaseIsolationTimeStartBase	# Start of case isolation policy.
	ParamList[["Duration of case isolation"]] 					= CaseIsolationDuration  		# CaseIsolationDuration = number of days a case will self-isolate. Different from CaseIsolationPolicyDuration = number of days policy is in effect.
	ParamList[["Duration of case isolation policy"]] 			= CaseIsolationPolicyDuration	# CaseIsolationDuration = number of days a case will self-isolate. Different from CaseIsolationPolicyDuration = number of days policy is in effect.
	ParamList[["Delay to start case isolation"]] 				= CaseIsolationDelay 		    # delay between case detection and case isolation
	
	ParamList[["Delay to case isolation by admin unit"]] 		= CI_Delay_byAdUnit
	ParamList[["Duration of case isolation by admin unit"]] 	= CI_PolicyDuration_byAdUnit   ### delay between case detection and case isolation by admin unit
	
	ParamList[["Case isolation trigger incidence per cell"]] 	= CaseIsolation_CellIncThresh
	ParamList[["Proportion of detected cases isolated"]] 		= CaseIsolationProp
	
	ParamList[["Residual contacts after case isolation"]] 				= CaseIsolationEffectiveness
	ParamList[["Residual household contacts after case isolation"]] 	= CaseIsolationHouseEffectiveness
	ParamList[["Number of change times for levels of case isolation"]] 	= Num_CI_ChangeTimes
	
	ParamList[["Change times for levels of case isolation"]] 					= CI_ChangeTimes #### Note: numbers here must match "Number of change times for levels of case isolation"; that any times listed here that are before "Case isolation start time" and after "Duration of case isolation policy" are irrelevant.
	ParamList[["Residual contacts after case isolation over time"]] 			= CI_SpatialAndPlaceEffects_OverTime
	ParamList[["Residual household contacts after case isolation over time"]] 	= CI_HouseholdEffects_OverTime
	ParamList[["Proportion of detected cases isolated over time"]] 				= CI_Prop_OverTime
	ParamList[["Case isolation trigger incidence per cell over time"]] 			= CI_CellIncThresh_OverTime
	
	
	# = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = 
	# = # = # = # = # = # = # = # = # = # = 	HOUSEHOLD QUARANTINE PARAMETERS
	
	ParamList[["Household quarantine start time"]] 			= HQuarantineTimeStartBase
	ParamList[["Delay to start household quarantine"]] 		= HQuarantineDelay
	ParamList[["Duration of household quarantine policy"]] 	= HQuarantinePolicyDuration
	
	ParamList[["Delay to household quarantine by admin unit"]] 						= HQuarantineDelay_byAdUnit
	ParamList[["Duration of household quarantine by admin unit"]] 					= HQuarantinePolicyDuration_byAdUnit
	
	ParamList[["Household quarantine trigger incidence per cell"]] 					= HHQuar_CellIncThresh
	ParamList[["Length of time households are quarantined"]] 						= HQuarantineHouseDuration
	ParamList[["Relative household contact rate after quarantine"]] 				= HQuarantineHouseEffect
	ParamList[["Residual place contacts after household quarantine by place type"]] = HQuarantinePlaceEffect
	ParamList[["Residual spatial contacts after household quarantine"]] 			= HQuarantineSpatialEffect
	ParamList[["Household level compliance with quarantine"]] 						= HQuarantinePropHouseCompliant
	ParamList[["Individual level compliance with quarantine"]] 						= HQuarantinePropIndivCompliant
	
	ParamList[["Number of change times for levels of household quarantine"]] 					= Num_HQ_ChangeTimes
	ParamList[["Change times for levels of household quarantine"]] 								= HQ_ChangeTimes
	ParamList[["Relative household contact rates over time after quarantine"]] 					= HQ_HouseholdEffects_OverTime
	ParamList[["Residual place contacts over time after household quarantine by place type"]] 	= HQ_PlaceEffects_OverTime
	ParamList[["Residual spatial contacts over time after household quarantine"]] 				= HQ_SpatialEffects_OverTime
	
	ParamList[["Household level compliance with quarantine over time"]] 		= HQ_Household_PropComply_OverTime
	ParamList[["Individual level compliance with quarantine over time"]] 		= HQ_Individual_PropComply_OverTime
	ParamList[["Household quarantine trigger incidence per cell over time"]] 	= HQ_CellIncThresh_OverTime
	
	
	
	
	# = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = 
	# = # = # = # = # = # = # = # = # = # = 	PLACE CLOSURE PARAMETERS
	
	ParamList[["Place closure start time"]] = PlaceCloseTimeStartBase
	ParamList[["Place closure second start time"]] = PlaceCloseTimeStartBase2
	ParamList[["Delay to start place closure"]] = PlaceCloseDelayMean
	ParamList[["Duration of place closure"]] = PlaceCloseDurationBase
	
	ParamList[["Delay to place closure by admin unit"]] = PlaceCloseDelayMean_byAdUnit
	ParamList[["Duration of place closure by admin unit"]] = PlaceCloseDuration_byAdUnit
	
	ParamList[["Place closure in administrative units rather than rings"]] 	= PlaceCloseByAdminUnit
	ParamList[["Administrative unit divisor for place closure"]] 			= PlaceCloseAdminUnitDivisor
	ParamList[["Place types to close for admin unit closure (0/1 array)"]] 	= PlaceCloseAdunitPlaceTypes
	ParamList[["Minimum radius for place closure"]] 						= PlaceCloseRadius
	
	ParamList[["Cumulative proportion of place members needing to become sick for admin unit closure"]] = PlaceCloseCasePropThresh
	ParamList[["Proportion of places in admin unit needing to pass threshold for place closure"]] 		= PlaceCloseAdunitPropThresh
	ParamList[["Proportion of places remaining open after closure by place type"]] 						= PropPlacesRemainingOpen_byPlaceType
	ParamList[["Proportional attendance after closure by place type"]] 									= PlaceClosePropAttending
	
	ParamList[["Relative household contact rate after closure"]] 										= PlaceCloseHouseholdRelContact
	ParamList[["Relative spatial contact rate after closure"]] 											= PlaceCloseSpatialRelContact
	
	ParamList[["Place closure incidence threshold"]] 			= PlaceCloseIncTrig 			## needs to be 0 for global triggers
	ParamList[["Place closure fractional incidence threshold"]] = PlaceCloseFracIncTrig			## needs to be 0 for global triggers or if abs incidence threshold used
	ParamList[["Trigger incidence per cell for place closure"]] = PlaceCloseCellIncThresh 		## change this for global too ###
	ParamList[["Trigger incidence per cell for end of place closure"]] = PlaceCloseCellIncStopThresh 		## change this for global too ###
	
	ParamList[["Number of change times for levels of place closure"]] 							= Num_PC_ChangeTimes
	ParamList[["Change times for levels of place closure"]] 									= PC_ChangeTimes #### Note: numbers here must match "Number of change times for levels of place closure"; that any times listed here that are before "Place closure start time" and after "Duration of place closure" are irrelevant.
	ParamList[["Proportion of places remaining open after closure by place type over time"]] 	= PC_PlaceEffects_OverTime 
	ParamList[["Relative household contact rates over time after place closure"]] 				= PC_HouseholdEffects_OverTime
	ParamList[["Relative spatial contact rates over time after place closure"]] 				= PC_SpatialEffects_OverTime
	ParamList[["Place closure incidence threshold over time"]] 									= PC_IncThresh_OverTime 
	ParamList[["Place closure fractional incidence threshold over time"]] 						= PC_FracIncThresh_OverTime 		
	ParamList[["Trigger incidence per cell for place closure over time"]] 						= PC_CellIncThresh_OverTime 		
	ParamList[["Duration of place closure over time"]] 											= PC_Durs_OverTime #### Note: closure durations longer than interval between change times will be truncated
	ParamList[["Proportional attendance after closure by place type over time" ]] 				= PC_PropAttending_OverTime
	
	# = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = 
	# = # = # = # = # = # = # = # = # = # = 	SOCIAL DISTANCING PARAMETERS
	
	ParamList[["Social distancing start time"]] 				= SocDistTimeStartBase
	ParamList[["Duration of social distancing"]] 				= SocDistDuration
	ParamList[["Delay to social distancing by admin unit"]] 	= SocDistDelay_byAdUnit
	ParamList[["Duration of social distancing by admin unit"]] 	= SocDistDuration_byAdUnit
	
	ParamList[["Trigger incidence per cell for social distancing"]] 					= SocDistCellIncThresh
	ParamList[["Trigger incidence per cell for end of social distancing"]] 			= SocDistCellIncStopThresh
	
	ParamList[["Relative place contact rate given social distancing by place type"]] 			= SocDistPlaceEffect
	ParamList[["Relative household contact rate given social distancing"]] 						= SocDistHouseholdEffect
	ParamList[["Relative spatial contact rate given social distancing"]] 						= SocDistSpatialEffect
	ParamList[["Relative place contact rate given enhanced social distancing by place type"]] 	= EnhancedSocDistPlaceEffect
	ParamList[["Relative household contact rate given enhanced social distancing"]] 			= EnhancedSocDistHouseholdEffect
	ParamList[["Relative spatial contact rate given enhanced social distancing"]] 				= EnhancedSocDistPlaceEffect
	
	ParamList[["Minimum radius for social distancing"]] = 1
	ParamList[["Proportion compliant with enhanced social distancing"]] 				= EnhancedSocDistProportionCompliantNonAge
	ParamList[["Proportion compliant with enhanced social distancing by age group"]] 	= EnhancedSocDistProportionCompliantByAge
	
	### "after change" parameters - now superseded by "over time" parameters. 
	#ParamList[["Delay for change in effectiveness of social distancing"]] 							= 1000
	#ParamList[["Relative place contact rate given social distancing by place type after change"]] 	= c(1, 1, 0.75, 0.75)
	#ParamList[["Relative household contact rate given social distancing after change"]] 			= 1.25
	#ParamList[["Relative spatial contact rate given social distancing after change"]] 				= 0.25
	
	ParamList[["Number of change times for levels of social distancing"]] 	= Num_SD_ChangeTimes #### Must match "Change times for levels of social distancing"
	ParamList[["Change times for levels of social distancing"]] 			= SD_ChangeTimes #### Note: numbers here must match "Number of change times for levels of social distancing"; that any times listed here that are before "Social distancing start time" and after "Duration of social distancing" are irrelevant.
	
	ParamList[["Relative place contact rates over time given social distancing by place type"]] 			= SD_PlaceEffects_OverTime #### Want this to supercede "Relative place contact rate given social distancing by place type". Should be matrix of dimension "Number of change times for levels of social distancing" by Number of place types.
	ParamList[["Relative place contact rates over time given enhanced social distancing by place type"]] 	= Enhanced_SD_PlaceEffects_OverTime
	ParamList[["Relative household contact rates over time given social distancing"]] 						= SD_HouseholdEffects_OverTime ####  Ideally want this to supercede "Relative household contact rate given social distancing" but need to preserved backwards compatibility for now.
	ParamList[["Relative household contact rates over time given enhanced social distancing"]] 				= Enhanced_SD_HouseholdEffects_OverTime
	ParamList[["Relative spatial contact rates over time given social distancing"]] 						= SD_SpatialEffects_OverTime ####  Ideally want this to supercede "Relative spatial contact rate given social distancing" but need to preserved backwards compatibility for now.
	ParamList[["Relative spatial contact rates over time given enhanced social distancing"]] 				= Enhanced_SD_SpatialEffects_OverTime
	ParamList[["Trigger incidence per cell for social distancing over time"]] 								= SD_CellIncThresh_OverTime
	
	
	return(ParamList)
}
