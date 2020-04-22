#pragma once

#ifndef COVIDSIM_PARAM_H_INCLUDED_
#define COVIDSIM_PARAM_H_INCLUDED_

#include "Country.h"
#include "Constants.h"

/**
 * @brief Stores the parameters for the simulation.
 *
 */
typedef struct PARAM {


	int N; /**< Population size */
	int NH; // Number of households
	int NR; /**< Number of Realisations */
	int NRN; /**< Number of non-extinct realisations */
	int NRactual;
	int NRactE;
	int NRactNE;
	int UpdatesPerSample; // Number of time steps between samples
	int NumSamples; // Total number of samples that will be made
	int KernelType;
	int MoveKernelType;
	int AirportKernelType;
	unsigned int BinFileLen;
	int DoBin, DoSaveSnapshot, DoLoadSnapshot;
	double SnapshotSaveTime, SnapshotLoadTime, clP1, clP2, clP3, clP4, clP5, clP6;
	int NC; // Number of cells
	int NMC; // Number of microcells
	int NMCL; // Number of microcells wide/high a cell is; i.e. NMC = NC * NMCL * NMCL
	int NCP; /**< Number of populated cells  */
	int NMCP, ncw, nch, nmcw, nmch, DoUTM_coords, nsp, DoSeasonality, DoCorrectAgeDist, DoPartialImmunity;
	int DoAdUnits, NumAdunits, DoAdunitBoundaries, AdunitLevel1Divisor, AdunitLevel1Mask, AdunitBitmapDivisor, CountryDivisor;
	int DoAdunitOutput, DoAdunitBoundaryOutput, DoAdunitDemog, DoCorrectAdunitPop, DoSpecifyPop, AdunitLevel1Lookup[ADUNIT_LOOKUP_SIZE];
	int DoOutputPlaceDistForOneAdunit, OutputPlaceDistAdunit, OutputDensFile;
	int DoOneGen, OutputEveryRealisation, BitmapMovieFrame, MaxCorrSample, DoLatent, InfQueuePeakLength, NumThreads, MaxNumThreads;
	int bwidth, bheight; // Size in pixels of the map area in the bitmap output
	int bheight2; // Height in pixels of the entire bitmap output, including both the spectrum at the top and the map area
	int bminx, bminy;
	int OutputBitmap; // Whether to output a bitmap
	int DoSI, DoHeteroDensity, DoPeriodicBoundaries, DoImmuneBitmap, OutputBitmapDetected; //added OutputBitmapDetected - ggilani 04/08/15
	int DoHouseholds, DoPlaces, PlaceTypeNum, Nplace[NUM_PLACE_TYPES], SmallEpidemicCases, DoPlaceGroupTreat;
	int NumInitialInfections[MAX_NUM_SEED_LOCATIONS], DoRandomInitialInfectionLoc, DoAllInitialInfectioninSameLoc;
	int MinPopDensForInitialInfection, NumSeedLocations, MaxPopDensForInitialInfection, InitialInfectionsAdminUnitId[MAX_NUM_SEED_LOCATIONS],InitialInfectionsAdminUnit[MAX_NUM_SEED_LOCATIONS];
	int DoAge, DoSymptoms, LoadSaveNetwork, IncThreshPop, GlobalIncThreshPop;
	int OutputOnlyNonExtinct, DoInfectiousnessProfile, DoInfectionTree, DoWholeHouseholdImmunity, DoSpatial, DoDeath;
	int DoAirports, Nairports, Air_popscale, DoSchoolFile, DoRealSymptWithdrawal, CaseAbsentChildAgeCutoff, DoEarlyCaseDiagnosis, DoInterventionFile;
	int PlaceTypeNoAirNum; // If DoAirports then this is the number of non-airport place types (< PlaceTypeNum), else == PlaceTypeNum (~ no airport places).
	int HotelPlaceType; // If DoAirports then this is place type for hotel (>= PlaceTypeNoAirNum, < PlaceTypeNum), else == PlaceTypeNum (~ unused).
	long seed1, seed2, seed3, seed4;
	long newseed1, newseed2; //added these to allow for seeds to be reset - ggilani 09/03/17
	int ResetSeeds,KeepSameSeeds, ResetSeedsPostIntervention, ResetSeedsFlag, TimeToResetSeeds;
	double SpatialBoundingBox[4], LocationInitialInfection[MAX_NUM_SEED_LOCATIONS][2], InitialInfectionsAdminUnitWeight[MAX_NUM_SEED_LOCATIONS], TimeStepsPerDay;
	double FalsePositiveRate, FalsePositivePerCapitaIncidence, FalsePositiveAgeRate[NUM_AGE_GROUPS];
	double latent_icdf[CDF_RES + 1], infectious_icdf[CDF_RES + 1], infectious_prof[INFPROF_RES + 1], infectiousness[MAX_INFECTIOUS_STEPS];

	double MildToRecovery_icdf[CDF_RES + 1], ILIToRecovery_icdf[CDF_RES + 1], SARIToRecovery_icdf[CDF_RES + 1], CriticalToCritRecov_icdf[CDF_RES + 1], CritRecovToRecov_icdf[CDF_RES + 1];
	double ILIToSARI_icdf[CDF_RES + 1], SARIToCritical_icdf[CDF_RES + 1], ILIToDeath_icdf[CDF_RES + 1], SARIToDeath_icdf[CDF_RES + 1], CriticalToDeath_icdf[CDF_RES + 1];
	/// means for above icdf's.
	double Mean_MildToRecovery, Mean_ILIToRecovery, Mean_SARIToRecovery, Mean_CriticalToCritRecov, Mean_CritRecovToRecov, Mean_TimeToTest, Mean_TimeToTestOffset, Mean_TimeToTestCriticalOffset, Mean_TimeToTestCritRecovOffset;
	double Mean_ILIToSARI, Mean_SARIToCritical, Mean_CriticalToDeath, Mean_SARIToDeath, Mean_ILIToDeath;
	double Prop_Mild_ByAge[NUM_AGE_GROUPS], Prop_ILI_ByAge[NUM_AGE_GROUPS], Prop_SARI_ByAge[NUM_AGE_GROUPS], Prop_Critical_ByAge[NUM_AGE_GROUPS];
	double CFR_SARI_ByAge[NUM_AGE_GROUPS], CFR_Critical_ByAge[NUM_AGE_GROUPS], CFR_ILI_ByAge[NUM_AGE_GROUPS];

	double TimeStep; // The length of a time step, in days
	double SampleTime; // The number of days to run for
	double SampleStep; // The length of a sampling step, in days
	double BitmapAspectScale; // Height of bitmap / Width of bitmap
	int ts_age;
	int DoSeverity; // Non-zero (true) if severity analysis should be done
	double scalex, scaley; // Number of pixels per degree in bitmap output
	double width, height; // Size of spatial domain in degrees
	double cwidth, cheight; // Size of spatial domain in cells
	double mcwidth, mcheight; // Size of spatial domain in microcells
	double KernelShape, KernelScale, KernelP3, KernelP4, KernelDelta, MoveKernelShape, MoveKernelScale, MoveKernelP3, MoveKernelP4;
	double AirportKernelShape, AirportKernelScale, AirportKernelP3, AirportKernelP4, AirportTrafficScale;
	double R0, R0scale, LocalBeta;
	double LatentPeriod; // In days. Mean of icdf (inverse cumulative distribution function).
	double InfectiousPeriod; // In days. Mean of icdf (inverse cumulative distribution function).
	double R0household, R0places, R0spatial;
	double Seasonality[DAYS_PER_YEAR];
	double InfectiousnessSD, R0DensityScalePower, InfectiousnessGamA, InfectiousnessGamR, SuscReductionFactorPerInfection;
	double ProportionSymptomatic[NUM_AGE_GROUPS], LatentToSymptDelay, SymptInfectiousness;
	double SymptSpatialContactRate, SymptPlaceTypeContactRate[NUM_PLACE_TYPES], InhibitInterAdunitPlaceAssignment[NUM_PLACE_TYPES];
	double SymptPlaceTypeWithdrawalProp[NUM_PLACE_TYPES], CaseAbsenteeismDuration, CaseAbsenteeismDelay;
	double CaseAbsentChildPropAdultCarers;
	double RelativeTravelRate[NUM_AGE_GROUPS], RelativeSpatialContact[NUM_AGE_GROUPS];
	double AgeSusceptibility[NUM_AGE_GROUPS], AgeInfectiousness[NUM_AGE_GROUPS], InitialImmunity[NUM_AGE_GROUPS];
	double WAIFW_Matrix[NUM_AGE_GROUPS][NUM_AGE_GROUPS];
	double HotelPropLocal, JourneyDurationDistrib[MAX_TRAVEL_TIME], LocalJourneyDurationDistrib[MAX_TRAVEL_TIME];
	double MeanJourneyTime, MeanLocalJourneyTime;
	int InvJourneyDurationDistrib[1025], InvLocalJourneyDurationDistrib[1025];
	double HouseholdTrans, HouseholdSizeDistrib[MAX_ADUNITS][MAX_HOUSEHOLD_SIZE], HouseholdTransPow;
	double HouseholdDenomLookup[MAX_HOUSEHOLD_SIZE];
	int PlaceTypeAgeMin[NUM_PLACE_TYPES], PlaceTypeAgeMax[NUM_PLACE_TYPES], PlaceTypeMaxAgeRead[NUM_PLACE_TYPES];
	int PlaceTypeAgeMin2[NUM_PLACE_TYPES], PlaceTypeAgeMax2[NUM_PLACE_TYPES];
	int PlaceTypeAgeMin3[NUM_PLACE_TYPES], PlaceTypeAgeMax3[NUM_PLACE_TYPES];
	int PlaceTypeNearestNeighb[NUM_PLACE_TYPES], PlaceTypeKernelType[NUM_PLACE_TYPES];
	double PlaceTypePropAgeGroup[NUM_PLACE_TYPES], PlaceTypePropAgeGroup2[NUM_PLACE_TYPES];
	double PlaceTypePropAgeGroup3[NUM_PLACE_TYPES], PlaceTypeKernelShape[NUM_PLACE_TYPES], PlaceTypeKernelScale[NUM_PLACE_TYPES];
	double PlaceTypeKernelP3[NUM_PLACE_TYPES], PlaceTypeKernelP4[NUM_PLACE_TYPES], PlaceTypeTrans[NUM_PLACE_TYPES];
	double PlaceTypeMeanSize[NUM_PLACE_TYPES], PlaceTypePropBetweenGroupLinks[NUM_PLACE_TYPES], PlaceTypeSizeSD[NUM_PLACE_TYPES]; //added PlaceTypeSizeSD for lognormal distribution - ggilani 09/02/17
	double PlaceTypeSizePower[NUM_PLACE_TYPES], PlaceTypeSizeOffset[NUM_PLACE_TYPES], PlaceTypeSizeMax[NUM_PLACE_TYPES];
	double PlaceTypeGroupSizeParam1[NUM_PLACE_TYPES], PlaceExclusivityMatrix[NUM_PLACE_TYPES * NUM_PLACE_TYPES]; //changed PlaceExclusivityMatrix from [NUM_PLACE_TYPES][NUM_PLACE_TYPES]
	double PropAgeGroup[MAX_ADUNITS][NUM_AGE_GROUPS], PopByAdunit[MAX_ADUNITS][2], MeanAnnualDeathRate;
	double InvLifeExpecDist[MAX_ADUNITS][1001];

	double PlaceCloseTimeStart, PlaceCloseTimeStart2, PlaceCloseDurationBase, PlaceCloseDuration, PlaceCloseDuration2, PlaceCloseDelayMean, PlaceCloseRadius, PlaceCloseRadius2;
	double PlaceCloseEffect[NUM_PLACE_TYPES], PlaceCloseSpatialRelContact, PlaceCloseHouseholdRelContact;
	double PlaceCloseCasePropThresh, PlaceCloseAdunitPropThresh, PlaceCloseFracIncTrig;
	int DoHolidays, NumHolidays;
	double HolidayEffect[NUM_PLACE_TYPES], HolidayStartTime[DAYS_PER_YEAR], HolidayDuration[DAYS_PER_YEAR];
	double ColourPeriod, BoundingBox[4], BitmapScale;
	double TreatSuscDrop, TreatInfDrop, TreatDeathDrop, TreatSympDrop, TreatDelayMean, TreatTimeStart, TreatPlaceGeogDuration;
	double TreatProphCourseLength, TreatCaseCourseLength, TreatPropRadial, TreatRadius, TreatRadius2, TreatCellIncThresh;
	double CaseIsolation_CellIncThresh, HHQuar_CellIncThresh, DigitalContactTracing_CellIncThresh;
	double TreatPropCases, TreatPropCaseHouseholds, TreatHouseholdsDuration;
	double TreatPlaceProbCaseId[NUM_PLACE_TYPES], TreatPlaceTotalProp[NUM_PLACE_TYPES];
	double TreatMaxCoursesBase, TreatNewCoursesRate, TreatNewCoursesStartTime, TreatMaxCourses;
	double VaccSuscDrop, VaccSuscDrop2, VaccInfDrop, VaccMortDrop, VaccSympDrop, VaccDelayMean, VaccTimeStart, VaccTimeEfficacySwitch, VaccTimeStartGeo;
	double VaccTimeToEfficacy, VaccProp, VaccRadius, VaccRadius2, VaccMinRadius, VaccMinRadius2, VaccPropCaseHouseholds, VaccHouseholdsDuration, VaccMaxCoursesBase;
	double VaccNewCoursesRate, VaccNewCoursesStartTime, VaccMaxCourses, VaccNewCoursesEndTime, VaccEfficacyDecay, VaccCellIncThresh, VaccCampaignInterval, VaccCoverageIncreasePeriod;
	int VaccDosePerDay;
	double PreAlertControlPropCasesId, PostAlertControlPropCasesId, ControlPropCasesId;
	double MoveRestrRadius, MoveRestrRadius2;
	double MoveDelayMean, MoveRestrEffect, MoveRestrDuration, MoveRestrTimeStart;
	double AirportCloseTimeStart, AirportCloseDuration, AirportCloseEffectiveness;

	double CaseIsolationTimeStart, CaseIsolationDuration, CaseIsolationEffectiveness, CaseIsolationHouseEffectiveness;
	double CaseIsolationDelay, CaseIsolationPolicyDuration, CaseIsolationProp;

	double HQuarantineTimeStart, HQuarantineHouseDelay, HQuarantineHouseDuration, HQuarantinePolicyDuration, HQuarantinePropIndivCompliant;
	double HQuarantinePropHouseCompliant, HQuarantinePlaceEffect[NUM_PLACE_TYPES], HQuarantineSpatialEffect, HQuarantineHouseEffect;

	double SocDistTimeStart, SocDistDuration, SocDistHouseholdEffect, SocDistPlaceEffect[NUM_PLACE_TYPES], SocDistSpatialEffect;
	double EnhancedSocDistHouseholdEffect, EnhancedSocDistPlaceEffect[NUM_PLACE_TYPES], EnhancedSocDistSpatialEffect, EnhancedSocDistProportionCompliant[NUM_AGE_GROUPS];

	double SocDistChangeDelay, SocDistDuration2, SocDistHouseholdEffect2, SocDistPlaceEffect2[NUM_PLACE_TYPES], SocDistSpatialEffect2;
	double EnhancedSocDistHouseholdEffect2, EnhancedSocDistPlaceEffect2[NUM_PLACE_TYPES], EnhancedSocDistSpatialEffect2;

	double SocDistDurationCurrent, SocDistHouseholdEffectCurrent, SocDistPlaceEffectCurrent[NUM_PLACE_TYPES], SocDistSpatialEffectCurrent;
	double EnhancedSocDistHouseholdEffectCurrent, EnhancedSocDistPlaceEffectCurrent[NUM_PLACE_TYPES], EnhancedSocDistSpatialEffectCurrent;

	double SocDistRadius, SocDistRadius2;

	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** 
	///// **** VARIABLE EFFICACIES OVER TIME
	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** 

	int VaryEfficaciesOverTime;

	/**< SOCIAL DISTANCING	*/
	/**< non-enhanced	*/
	int Num_SD_ChangeTimes; //// must be at most MAX_NUM_INTERVENTION_CHANGE_TIMES
	int SD_ChangeTimes						[MAX_NUM_INTERVENTION_CHANGE_TIMES]; /**< change times for intensity of (enhanced) social distancing */
	double SD_SpatialEffects_OverTime		[MAX_NUM_INTERVENTION_CHANGE_TIMES]; //// time-varying equivalent of SocDistSpatialEffectCurrent
	double SD_HouseholdEffects_OverTime		[MAX_NUM_INTERVENTION_CHANGE_TIMES]; //// time-varying equivalent of SocDistHouseholdEffectCurrent
	double SD_PlaceEffects_OverTime			[MAX_NUM_INTERVENTION_CHANGE_TIMES][NUM_PLACE_TYPES];	//// indexed by i) change time; ii) place type;  //// time-varying equivalent of SocDistPlaceEffectCurrent
	double SD_PlaceEffects_OverTime_dummy	[MAX_NUM_INTERVENTION_CHANGE_TIMES * NUM_PLACE_TYPES];	//// Hack as can't read in matrix above. Will read this in then populate the above.
	int SD_CellIncThresh_OverTime			[MAX_NUM_INTERVENTION_CHANGE_TIMES]; //// time-varying equivalent of SocDistCellIncThresh

	/**< enhanced	*/
	double Enhanced_SD_SpatialEffects_OverTime		[MAX_NUM_INTERVENTION_CHANGE_TIMES]; //// time-varying equivalent of EnhancedSocDistSpatialEffectCurrent
	double Enhanced_SD_HouseholdEffects_OverTime	[MAX_NUM_INTERVENTION_CHANGE_TIMES]; //// time-varying equivalent of EnhancedSocDistHouseholdEffectCurrent
	double Enhanced_SD_PlaceEffects_OverTime		[MAX_NUM_INTERVENTION_CHANGE_TIMES][NUM_PLACE_TYPES];	//// indexed by i) change time; ii) place type;  time-varying equivalent of EnhancedSocDistPlaceEffectCurrent
	double Enhanced_SD_PlaceEffects_OverTime_dummy	[MAX_NUM_INTERVENTION_CHANGE_TIMES * NUM_PLACE_TYPES];	//// Hack as can't read in matrix above. Will read this in then populate the above.  

	/**< CASE ISOLATION	*/
	int Num_CI_ChangeTimes; //// must be at most MAX_NUM_INTERVENTION_CHANGE_TIMES
	int CI_ChangeTimes							[MAX_NUM_INTERVENTION_CHANGE_TIMES]; /**< change times for intensity of case isolation */
	double CI_SpatialAndPlaceEffects_OverTime	[MAX_NUM_INTERVENTION_CHANGE_TIMES]; //// time-varying equivalent of CaseIsolationEffectiveness
	double CI_HouseholdEffects_OverTime			[MAX_NUM_INTERVENTION_CHANGE_TIMES]; //// time-varying equivalent of CaseIsolationHouseEffectiveness
	double CI_Prop_OverTime						[MAX_NUM_INTERVENTION_CHANGE_TIMES]; //// time-varying equivalent of CaseIsolationProp
	double CI_CellIncThresh_OverTime				[MAX_NUM_INTERVENTION_CHANGE_TIMES]; //// time-varying equivalent of CaseIsolation_CellIncThresh

	/**< HOUSEHOLD QUARANTINE	*/
	int Num_HQ_ChangeTimes; //// must be at most MAX_NUM_INTERVENTION_CHANGE_TIMES
	int HQ_ChangeTimes							[MAX_NUM_INTERVENTION_CHANGE_TIMES]; /**< change times for intensity of household quarantine */
	double HQ_SpatialEffects_OverTime			[MAX_NUM_INTERVENTION_CHANGE_TIMES]; //// time-varying equivalent of HQuarantineSpatialEffect
	double HQ_HouseholdEffects_OverTime			[MAX_NUM_INTERVENTION_CHANGE_TIMES]; //// time-varying equivalent of HQuarantineHouseEffect
	double HQ_PlaceEffects_OverTime				[MAX_NUM_INTERVENTION_CHANGE_TIMES][NUM_PLACE_TYPES];	//// indexed by i) change time; ii) place type; time-varying equivalent of HQuarantinePlaceEffect
	double HQ_PlaceEffects_OverTime_dummy		[MAX_NUM_INTERVENTION_CHANGE_TIMES * NUM_PLACE_TYPES];	//// Hack as can't read in matrix above. Will read this in then populate the above.
	double HQ_Individual_PropComply_OverTime	[MAX_NUM_INTERVENTION_CHANGE_TIMES]; //// time-varying equivalent of HQuarantinePropIndivCompliant
	double HQ_Household_PropComply_OverTime		[MAX_NUM_INTERVENTION_CHANGE_TIMES]; //// time-varying equivalent of HQuarantinePropHouseCompliant
  double HQ_CellIncThresh_OverTime				[MAX_NUM_INTERVENTION_CHANGE_TIMES]; //// time-varying equivalent of HHQuar_CellIncThresh

	/**< PLACE CLOSURE	*/
	int Num_PC_ChangeTimes; //// must be at most MAX_NUM_INTERVENTION_CHANGE_TIMES
	int PC_ChangeTimes						[MAX_NUM_INTERVENTION_CHANGE_TIMES]; /**< change times for intensity of place closure */
	double PC_SpatialEffects_OverTime		[MAX_NUM_INTERVENTION_CHANGE_TIMES]; //// time-varying equivalent of PlaceCloseSpatialRelContact
	double PC_HouseholdEffects_OverTime		[MAX_NUM_INTERVENTION_CHANGE_TIMES]; //// time-varying equivalent of PlaceCloseHouseholdRelContact
	double PC_PlaceEffects_OverTime			[MAX_NUM_INTERVENTION_CHANGE_TIMES][NUM_PLACE_TYPES];	//// indexed by i) change time; ii) place type; //// time-varying equivalent of PlaceCloseEffect
	double PC_PlaceEffects_OverTime_dummy	[MAX_NUM_INTERVENTION_CHANGE_TIMES * NUM_PLACE_TYPES];	//// Hack as can't read in matrix above. Will read this in then populate the above.
	int PC_IncThresh_OverTime				[MAX_NUM_INTERVENTION_CHANGE_TIMES]; //// time-varying equivalent of PlaceCloseIncTrig / PlaceCloseIncTrig1
	double PC_FracIncThresh_OverTime		[MAX_NUM_INTERVENTION_CHANGE_TIMES]; //// time-varying equivalent of PlaceCloseFracIncTrig
	int PC_CellIncThresh_OverTime			[MAX_NUM_INTERVENTION_CHANGE_TIMES]; //// time-varying equivalent of PlaceCloseCellIncThresh


	/**< DIGITAL CONTACT TRACING	*/
	int Num_DCT_ChangeTimes; //// must be at most MAX_NUM_INTERVENTION_CHANGE_TIMES
	int DCT_ChangeTimes							[MAX_NUM_INTERVENTION_CHANGE_TIMES]; /**< change times for intensity of digital contact tracing */
	double DCT_SpatialAndPlaceEffects_OverTime	[MAX_NUM_INTERVENTION_CHANGE_TIMES]; //// time-varying equivalent of DCTCaseIsolationEffectiveness
	double DCT_HouseholdEffects_OverTime		[MAX_NUM_INTERVENTION_CHANGE_TIMES]; //// time-varying equivalent of DCTCaseIsolationHouseEffectiveness
	double DCT_Prop_OverTime					[MAX_NUM_INTERVENTION_CHANGE_TIMES]; //// time-varying equivalent of ProportionDigitalContactsIsolate

	double KeyWorkerProphTimeStart, KeyWorkerProphDuration, KeyWorkerPropInKeyPlaces[NUM_PLACE_TYPES], KeyWorkerHouseProp;
	double KeyWorkerProphRenewalDuration, KeyWorkerProphRadius, KeyWorkerProphRadius2;

	double TreatTimeStartBase, VaccTimeStartBase, MoveRestrTimeStartBase, PlaceCloseTimeStartBase, PlaceCloseTimeStartBase2,PlaceCloseTimeStartPrevious;
	double AirportCloseTimeStartBase, HQuarantineTimeStartBase, CaseIsolationTimeStartBase, SocDistTimeStartBase, KeyWorkerProphTimeStartBase, DigitalContactTracingTimeStartBase;
	double InfectionImportRate1, InfectionImportRate2, InfectionImportChangeTime, ImportInfectionTimeProfile[MAX_DUR_IMPORT_PROFILE];
	double PreControlClusterIdTime, PreControlClusterIdCalTime, PreControlClusterIdHolOffset, PreIntervIdCalTime,PreIntervTime,SeedingScaling;
	int PreControlClusterIdCaseThreshold, PreControlClusterIdUseDeaths, PreControlClusterIdDuration, DoAlertTriggerAfterInterv, AlertTriggerAfterIntervThreshold,StopCalibration,ModelCalibIteration;
	int DoPerCapitaTriggers, DoGlobalTriggers, DoAdminTriggers, DoICUTriggers, MoveRestrCellIncThresh, DoHQretrigger;

	int PlaceCloseCellIncThresh, PlaceCloseCellIncThresh1, PlaceCloseCellIncThresh2, TriggersSamplingInterval, PlaceCloseIndepThresh, SocDistCellIncThresh, VaccPriorityGroupAge[2];
	int PlaceCloseCellIncStopThresh, SocDistCellIncStopThresh;
	int PlaceCloseAdunitPlaceTypes[NUM_PLACE_TYPES];

	int DoPlaceCloseOnceOnly, DoSocDistOnceOnly, DoMoveRestrOnceOnly, DoKeyWorkerProphOnceOnly;

	int VaccMaxRounds, VaccByAdminUnit, VaccAdminUnitDivisor, TreatByAdminUnit, TreatAdminUnitDivisor, MoveRestrByAdminUnit, MoveRestrAdminUnitDivisor, PlaceCloseByAdminUnit, PlaceCloseAdminUnitDivisor;
	int KeyWorkerProphCellIncThresh, KeyWorkerPopNum, KeyWorkerPlaceNum[NUM_PLACE_TYPES], KeyWorkerNum, KeyWorkerIncHouseNum;
	int DoBlanketMoveRestr, PlaceCloseIncTrig, PlaceCloseIncTrig1, PlaceCloseIncTrig2, TreatMaxCoursesPerCase, DoImportsViaAirports, DoMassVacc, DurImportTimeProfile;
	unsigned short int usHQuarantineHouseDuration, usVaccTimeToEfficacy, usVaccTimeEfficacySwitch; //// us = unsigned short versions of their namesakes, multiplied by P.TimeStepsPerDay
	unsigned short int usCaseIsolationDuration, usCaseIsolationDelay, usCaseAbsenteeismDuration, usCaseAbsenteeismDelay;

	int DoDistributionVaccination;
	//Added DoRecordInfEvents and MaxInfEvents in order to give the user a choice as to whether to output infection events as a line list: ggilani - 10/10/14
	int DoRecordInfEvents, MaxInfEvents, RecordInfEventsPerRun;
	double KernelPowerScale, KernelOffsetScale;
	int LimitNumInfections, MaxNumInfections;

	//Added parameters to deal with digital contact tracing - ggilani 09/03/2020
	int DoDigitalContactTracing, ClusterDigitalContactUsers, NDigitalContactUsers, NDigitalHouseholdUsers, FindContactsOfDCTContacts, DoDCTTest;
	double PropPopUsingDigitalContactTracing, ScalingFactorSpatialDigitalContacts, ScalingFactorPlaceDigitalContacts, DigitalContactTracingDelay, LengthDigitalContactIsolation, ProportionDigitalContactsIsolate, ProportionSmartphoneUsersByAge[NUM_AGE_GROUPS];
	double DelayFromIndexCaseDetectionToDCTIsolation, DelayToTestIndexCase, DelayToTestDCTContacts, SpecificityDCT, SensitivityDCT;
	double DigitalContactTracingPolicyDuration, DCTCaseIsolationHouseEffectiveness, DCTCaseIsolationEffectiveness;
	int OutputDigitalContactTracing, OutputDigitalContactDist, IncludeHouseholdDigitalContactTracing, IncludePlaceGroupDigitalContactTracing, DCTIsolateIndexCases, RemoveContactsOfNegativeIndexCase;

	int DoOriginDestinationMatrix; //added: ggilani 28/01/15
	int DoInterventionDelaysByAdUnit;


	int OutputAge, OutputR0, OutputControls, OutputCountry, OutputAdUnitVar, OutputHousehold, OutputInfType, OutputNonSeverity, OutputSeverityAdminUnit, OutputNonSummaryResults;

	int MeanChildAgeGap; // Average gap between ages of children in a household, in years
	int MinAdultAge; // The youngest age, in years, at which someone is considered to be an adult
	int MaxMFPartnerAgeGap; // The largest number of years older than a female partner that a male partner can be
	int MaxFMPartnerAgeGap; // The largest number of years older than a male partner that a female partner can be
	int MinParentAgeGap; // The minimum number of years older than a child that a parent must be
	int MaxParentAgeGap; // The maximum number of years older than a child that a parent can be
	int MaxChildAge; // The maximum age, in years, of a child
	double OneChildTwoPersProb;
	double TwoChildThreePersProb;
	double OnePersHouseProbOld;
	double TwoPersHouseProbOld;
	double OnePersHouseProbYoung;
	double TwoPersHouseProbYoung;
	double OneChildProbYoungestChildUnderFive;
	double TwoChildrenProbYoungestUnderFive;
	double ProbYoungestChildUnderFive;
	double ZeroChildThreePersProb;
	double OneChildFourPersProb;
	double YoungAndSingleSlope;
	int YoungAndSingle;
	int NoChildPersAge;
	int OldPersAge;
	double ThreeChildFivePersProb;
	int OlderGenGap;
} param;

extern param P;

#endif // COVIDSIM_PARAM_H_INCLUDED_
