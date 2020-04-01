#pragma once

#ifndef SPATIALSIM_PARAM_H_INCLUDED_
#define SPATIALSIM_PARAM_H_INCLUDED_

#include "Country.h"
#include "Constants.h"

/**
 * @brief Stores the parameters for the simulation.
 * 
 */
typedef struct PARAM {


	int N; /**< Population size */
	int NH;
	int NR; /**< Number of Realisations */
	int NRN; /**< Number of non-extinct realisations */
	int NRactual;
	int NRactE;
	int NRactNE;
	int UpdatesPerSample;
	int NumSamples;
	int KernelType;
	int MoveKernelType;
	int AirportKernelType;
	unsigned int BinFileLen;
	int DoBin, DoSaveSnapshot, DoLoadSnapshot, DoTimeSeries;
	double SnapshotSaveTime, SnapshotLoadTime, clP1, clP2, clP3, clP4, clP5, clP6;
	int NC, NMC, NMCL; 
	int NCP; /**< Number of populated cells  */
	int NMCP, ncw, nch, nmcw, nmch, DoUTM_coords, nsp, DoSIS, DoInitEquilib, DoSeasonality,DoCorrectAgeDist;
	int DoAdUnits, NumAdunits, DoAdunitBoundaries, AdunitLevel1Divisor, AdunitLevel1Mask, AdunitBitmapDivisor, CountryDivisor;
	int DoAdunitOutput, DoAdunitBoundaryOutput, DoAdunitDemog, DoCorrectAdunitPop, DoSpecifyPop, AdunitLevel1Lookup[ADUNIT_LOOKUP_SIZE];
	int DoOutputPlaceDistForOneAdunit, OutputPlaceDistAdunit, OutputDensFile;
	int DoOneGen, OutputEveryRealisation, BitmapMovieFrame, MaxCorrSample, DoLatent, InfQueuePeakLength, NumThreads, MaxNumThreads;
	int bwidth, bheight, bheight2, bminx, bminy, OutputBitmap, DoSI, DoHeteroDensity, DoPeriodicBoundaries, DoImmuneBitmap, OutputBitmapDetected; //added OutputBitmapDetected - ggilani 04/08/15
	int DoHouseholds, DoPlaces, PlaceTypeNum, Nplace[NUM_PLACE_TYPES], SmallEpidemicCases, DoPlaceGroupTreat;
	int NumInitialInfections[MAX_NUM_SEED_LOCATIONS], DoRandomInitialInfectionLoc, DoAllInitialInfectioninSameLoc;
	int MinPopDensForInitialInfection, NumSeedLocations, MaxPopDensForInitialInfection, InitialInfectionsAdminUnitId[MAX_NUM_SEED_LOCATIONS],InitialInfectionsAdminUnit[MAX_NUM_SEED_LOCATIONS];
	int DoAge, DoSymptoms, LoadSaveNetwork, IncThreshPop, GlobalIncThreshPop;
	int OutputOnlyNonExtinct, DoInfectiousnessProfile, DoInfectionTree, DoWholeHouseholdImmunity, DoSpatial, DoDeath, UpdatesPerDemogUpdate;
	int DoAirports, Nairports, Air_popscale, DoSchoolFile, DoRealSymptWithdrawal, CaseAbsentChildAgeCutoff, DoEarlyCaseDiagnosis, DoInterventionFile;
	long seed1, seed2, seed3, seed4;
	long newseed1, newseed2, newseed3, newseed4; //added these to allow for seeds to be reset - ggilani 09/03/17
	int ResetSeeds,KeepSameSeeds, ResetSeedsPostIntervention, ResetSeedsFlag, TimeToResetSeeds;
	double SpatialBoundingBox[4], LocationInitialInfection[MAX_NUM_SEED_LOCATIONS][2], InitialInfectionsAdminUnitWeight[MAX_NUM_SEED_LOCATIONS], TimeStepsPerDay, TimeStepsPerYear;
	double FalsePositiveRate, FalsePositivePerCapitaIncidence, FalsePositiveAgeRate[NUM_AGE_GROUPS];
	double latent_icdf[CDF_RES + 1], infectious_icdf[CDF_RES + 1], infectious_prof[INFPROF_RES + 1], infectiousness[MAX_INFECTIOUS_STEPS];

	double MildToRecovery_icdf[CDF_RES + 1], ILIToRecovery_icdf[CDF_RES + 1], SARIToRecovery_icdf[CDF_RES + 1], CriticalToCritRecov_icdf[CDF_RES + 1], CritRecovToRecov_icdf[CDF_RES + 1];
	double ILIToSARI_icdf[CDF_RES + 1], SARIToCritical_icdf[CDF_RES + 1], CriticalToDeath_icdf[CDF_RES + 1];
	/// means for above icdf's. 
	double Mean_MildToRecovery, Mean_ILIToRecovery, Mean_SARIToRecovery, Mean_CriticalToCritRecov, Mean_CritRecovToRecov;
	double Mean_ILIToSARI, Mean_SARIToCritical, Mean_CriticalToDeath;
	double Prop_Mild_ByAge[NUM_AGE_GROUPS], Prop_ILI_ByAge[NUM_AGE_GROUPS], Prop_SARI_ByAge[NUM_AGE_GROUPS], Prop_Critical_ByAge[NUM_AGE_GROUPS];
	double CFR_SARI_ByAge[NUM_AGE_GROUPS], CFR_Critical_ByAge[NUM_AGE_GROUPS];

	double T, TimeStep, SampleTime, SampleStep, BitmapAspectScale;
	int ts_age, TimeStepsPerDayInt, DoSeverity;
	double scalex, scaley, width, height, cwidth, cheight, mcwidth, mcheight;
	double KernelShape, KernelScale, KernelP3, KernelP4, KernelDelta, MoveKernelShape, MoveKernelScale, MoveKernelP3, MoveKernelP4;
	double AirportKernelShape, AirportKernelScale, AirportKernelP3, AirportKernelP4, AirportTrafficScale, CellPop2;
	double R0, R0scale, ContactsPerDay, LocalBeta, LatentPeriod, InfectiousPeriod, R0household, R0places, R0spatial;	//// LatentPeriod and InfectiousPeriod are means of icdf's (inverse cumulative distribution functions). 
	double Seasonality[DAYS_PER_YEAR];
	double InfectiousnessSD, R0DensityScalePower, InfectiousnessGamA, InfectiousnessGamR, SuscReductionFactorPerInfection, InfectiousnessBetaA, InfectiousnessBetaB;
	double LethalInfectiousPeriod, ProportionSymptomatic[NUM_AGE_GROUPS], LatentToSymptDelay, SymptInfectiousness;
	double SymptSpatialContactRate, SymptPlaceTypeContactRate[NUM_PLACE_TYPES], InhibitInterAdunitPlaceAssignment[NUM_PLACE_TYPES];
	double SymptPlaceTypeWithdrawalProp[NUM_PLACE_TYPES], CaseAbsenteeismDuration, CaseAbsenteeismDelay;
	double CaseAbsentChildPropAdultCarers;
	double DiseaseMortality;
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
	double MortalityByAge[MAX_ADUNITS][NUM_AGE_GROUPS * AGE_GROUP_WIDTH], CumulPropDead[MAX_ADUNITS][NUM_AGE_GROUPS * AGE_GROUP_WIDTH + 1], InvLifeExpecDist[MAX_ADUNITS][1001];
	double PlaceCloseTimeStart, PlaceCloseTimeStart2, PlaceCloseDurationBase, PlaceCloseDuration, PlaceCloseDuration2, PlaceCloseDelayMean, PlaceCloseRadius, PlaceCloseRadius2;
	double PlaceCloseEffect[NUM_PLACE_TYPES], PlaceCloseSpatialRelContact, PlaceCloseHouseholdRelContact;
	double PlaceCloseCasePropThresh, PlaceCloseAdunitPropThresh, PlaceCloseFracIncTrig;
	int DoHolidays, NumHolidays;
	double HolidayEffect[NUM_PLACE_TYPES], HolidayStartTime[DAYS_PER_YEAR], HolidayDuration[DAYS_PER_YEAR];
	double ColourPeriod, BoundingBox[4], BitmapScale, BitmapStartTime, BitmapStopTime;
	double TreatSuscDrop, TreatInfDrop, TreatDeathDrop, TreatSympDrop, TreatDelayMean, TreatTimeStart, TreatPlaceGeogDuration;
	double TreatProphCourseLength, TreatCaseCourseLength, TreatPropRadial, TreatRadius, TreatRadius2, TreatCellIncThresh; 
	double CaseIsolation_CellIncThresh, HHQuar_CellIncThresh, DigitalContactTracing_CellIncThresh;
	double TreatPropCases, TreatPropCaseHouseholds, TreatHouseholdsDuration;
	double TreatPlaceProbCaseId[NUM_PLACE_TYPES], TreatPlaceTotalProp[NUM_PLACE_TYPES];
	double TreatMaxCoursesBase, TreatNewCoursesRate, TreatNewCoursesStartTime, TreatMaxCourses;
	double VaccSuscDrop, VaccSuscDrop2, VaccInfDrop, VaccMortDrop, VaccSympDrop, VaccDelayMean, VaccTimeStart, VaccTimeEfficacySwitch, VaccTimeStartGeo;
	double VaccTimeToEfficacy, VaccProp, VaccRadius, VaccRadius2, VaccMinRadius, VaccMinRadius2, VaccPropCaseHouseholds, VaccHouseholdsDuration, VaccMaxCoursesBase;
	double VaccNewCoursesRate, VaccNewCoursesStartTime, VaccMaxCourses, VaccNewCoursesEndTime, VaccEfficacyDecay, VaccCellIncThresh, VaccCampaignInterval, VaccCoverageIncreasePeriod;
	int MinVaccAge, VaccDosePerDay; 
	double PreAlertControlPropCasesId, PostAlertControlPropCasesId, ControlPropCasesId;
	double MoveRestrRadius, MoveRestrRadius2;
	double MoveDelayMean, MoveRestrEffect, MoveRestrDuration, MoveRestrTimeStart;
	double AirportCloseTimeStart, AirportCloseDuration, AirportCloseEffectiveness;
	double CaseIsolationTimeStart, CaseIsolationDuration, CaseIsolationEffectiveness, CaseIsolationHouseEffectiveness;
	double CaseIsolationDelay, CaseIsolationPolicyDuration, CaseIsolationProp;
	double HQuarantineTimeStart, HQuarantineHouseDelay, HQuarantineHouseDuration, HQuarantinePolicyDuration, HQuarantinePropIndivCompliant;
	double HQuarantinePropHouseCompliant, HQuarantinePlaceEffect[NUM_PLACE_TYPES], HQuarantineSpatialEffect, HQuarantineHouseEffect;

	double SocDistTimeStart, SocDistDuration, SocDistHouseholdEffect, SocDistPlaceEffect[NUM_PLACE_TYPES], SocDistSpatialEffect;
	double ESocDistHouseholdEffect, ESocDistPlaceEffect[NUM_PLACE_TYPES], ESocDistSpatialEffect, ESocProportionCompliant[NUM_AGE_GROUPS];

	double SocDistChangeDelay, SocDistDuration2, SocDistHouseholdEffect2, SocDistPlaceEffect2[NUM_PLACE_TYPES], SocDistSpatialEffect2;
	double ESocDistHouseholdEffect2, ESocDistPlaceEffect2[NUM_PLACE_TYPES], ESocDistSpatialEffect2;

	double SocDistDurationC, SocDistHouseholdEffectC, SocDistPlaceEffectC[NUM_PLACE_TYPES], SocDistSpatialEffectC;
	double ESocDistHouseholdEffectC, ESocDistPlaceEffectC[NUM_PLACE_TYPES], ESocDistSpatialEffectC;

	double SocDistRadius, SocDistRadius2;

	double KeyWorkerProphTimeStart, KeyWorkerProphDuration, KeyWorkerPropInKeyPlaces[NUM_PLACE_TYPES], KeyWorkerHouseProp;
	double KeyWorkerProphRenewalDuration, KeyWorkerProphRadius, KeyWorkerProphRadius2;
	double TreatTimeStartBase, VaccTimeStartBase, MoveRestrTimeStartBase, PlaceCloseTimeStartBase, PlaceCloseTimeStartBase2;
	double AirportCloseTimeStartBase, HQuarantineTimeStartBase, CaseIsolationTimeStartBase, SocDistTimeStartBase, KeyWorkerProphTimeStartBase, DigitalContactTracingTimeStartBase;;
	double InfectionImportRate1, InfectionImportRate2, InfectionImportChangeTime, ImportInfectionTimeProfile[MAX_DUR_IMPORT_PROFILE];
	double PreControlClusterIdTime, PreControlClusterIdCalTime, PreControlClusterIdHolOffset;
	int PreControlClusterIdCaseThreshold, PreControlClusterIdUseDeaths, PreControlClusterIdDuration;
	int DoPerCapitaTriggers, DoGlobalTriggers, DoAdminTriggers, DoICUTriggers, MoveRestrCellIncThresh, DoHQretrigger;
	int PlaceCloseCellIncThresh, TriggersSamplingInterval, PlaceCloseIndepThresh, SocDistCellIncThresh, VaccPriorityGroupAge[2];
	int PlaceCloseCellIncStopThresh, SocDistCellIncStopThresh;
	int PlaceCloseAdunitPlaceTypes[NUM_PLACE_TYPES], DoPlaceCloseOnceOnly, DoSocDistOnceOnly, DoMoveRestrOnceOnly, DoKeyWorkerProphOnceOnly;
	int VaccMaxRounds, VaccByAdminUnit, VaccAdminUnitDivisor, TreatByAdminUnit, TreatAdminUnitDivisor, MoveRestrByAdminUnit, MoveRestrAdminUnitDivisor, PlaceCloseByAdminUnit, PlaceCloseAdminUnitDivisor;
	int KeyWorkerProphCellIncThresh, KeyWorkerPopNum, KeyWorkerPlaceNum[NUM_PLACE_TYPES], KeyWorkerNum, KeyWorkerIncHouseNum;
	int DoBlanketMoveRestr, PlaceCloseIncTrig, TreatMaxCoursesPerCase, DoImportsViaAirports, DoMassVacc, DurImportTimeProfile;
	unsigned short int usHQuarantineHouseDuration, usVaccTimeToEfficacy, usVaccTimeEfficacySwitch; //// us = unsigned short versions of their namesakes, multiplied by P.TimeStepsPerDay
	unsigned short int usCaseIsolationDuration, usCaseIsolationDelay, usCaseAbsenteeismDuration, usCaseAbsenteeismDelay;

	double RoutineImmunisationStartTime, RoutineImmunisationEffectiveCoverage, RoutineImmunisationMinAge, RoutineImmunisationMaxAge, RoutineCoverageNonDistrCountry;
	double SIAMinAge, SIAMaxAge, SIAStartTime, SIARepeatInterval, SIAEffectiveCoverage, SIADuration, VaccPhialLifetime;
	unsigned short int usRoutineImmunisationStartTime, usSIAStartTime, usRoutineImmunisationMaxAge, usRoutineImmunisationMinAge;
	unsigned short int usSIAMaxAge, usSIAMinAge, usSIARepeatInterval, usSIADuration;
	int DoDistributionVaccination, SIADoAllCountries, VaccDosesPerPhial;
	//Added DoRecordInfEvents and MaxInfEvents in order to give the user a choice as to whether to output infection events as a line list: ggilani - 10/10/14
	int DoRecordInfEvents, MaxInfEvents, RecordInfEventsPerRun;
	double KernelPowerScale, KernelOffsetScale;
	int LimitNumInfections, MaxNumInfections;

	//Added parameters to deal with digital contact tracing - ggilani 09/03/2020
	int DoDigitalContactTracing, ClusterDigitalContactUsers, NDigitalContactUsers, NDigitalHouseholdUsers;
	double PropPopUsingDigitalContactTracing, ScalingFactorPlaceDigitalContacts, ScalingFactorSpatialDigitalContacts, DigitalContactTracingDelay, LengthDigitalContactIsolation, ProportionDigitalContactsIsolate, ProportionSmartphoneUsersByAge[NUM_AGE_GROUPS];
	double DigitalContactTracingPolicyDuration;
	int OutputDigitalContactTracing;
	int DoOriginDestinationMatrix; //added: ggilani 28/01/15
	double KernelRadiusMax, KernelRadiusMax2, PropWithinCellTransmission;
	int DoUpdateCaseDetection;
	double TimeToUpdateCaseDetection, UpdatedCaseDetectionRate;
	int DoInterventionDelaysByAdUnit;
	

	int OutputAge, OutputR0, OutputControls, OutputCountry, OutputAdUnitVar, OutputHousehold, OutputInfType, OutputNonSeverity, OutputSeverityAdminUnit, OuptutNonSummaryResults; 
} param;

extern param P;

#endif // SPATIALSIM_PARAM_H_INCLUDED_
