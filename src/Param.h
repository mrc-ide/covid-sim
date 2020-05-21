#ifndef COVIDSIM_PARAM_H_INCLUDED_
#define COVIDSIM_PARAM_H_INCLUDED_

#include "Country.h"
#include "Constants.h"
#include "MicroCellPosition.hpp"

/** @brief Enumeration of bitmap formats. */
enum BitmapFormats
{
  BF_PNG = 0,  // PNG - default if IMAGE_MAGICK or _WIN32 defined
  BF_BMP = 1   // BMP - fall-back
};

/// Size of spatial domain in various units
struct DomainSize
{
	/// The width
	double width_;

	/// The height
	double height_;
};

struct SocialDistancing
{
  double spatial_effect_;

  double household_effect_;

  /// indexed by place type;
  double place_effect_[NUM_PLACE_TYPES];

  int cell_inc_thresh_;

  // enhanced
  double enhanced_spatial_effect_;

  double enhanced_household_effect_;

  /// indexed by place type;
  double enhanced_place_effect_[NUM_PLACE_TYPES];
};

struct SocialDistancings
{
  // non-enhanced
  /// must be at most MAX_NUM_INTERVENTION_CHANGE_TIMES
  int num_change_times_;

  /// change times for intensity of (enhanced) social distancing
  double change_times_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  /// time-varying equivalent of current_social_distancing_.spatial_effect_
  double spatial_effects_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  /// time-varying equivalent of current_social_distancing_.household_effect_
  double household_effects_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  /// indexed by i) change time; ii) place type;
  /// time-varying equivalent of current_social_distancing_.place_effect_
  double place_effects_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES][NUM_PLACE_TYPES];

  /// time-varying equivalent of current_social_distancing_.cell_inc_thresh_
  int cell_inc_thresh_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  // enhanced
  /// time-varying equivalent of current_social_distancing_.enhanced_spatial_effect_
  double enhanced_spatial_effects_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  /// time-varying equivalent of current_social_distancing_.enhanced_household_effect_
  double enhanced_household_effects_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  /// indexed by i) change time; ii) place type;
  /// time-varying equivalent of current_social_distancing_.enhanced_place_effect_
  double enhanced_place_effects_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES][NUM_PLACE_TYPES];
};

struct CaseIsolation
{
  double spatial_and_place_effect_;

  double household_effect_;

  double prop_;

  double cell_inc_thresh_;
};

struct CaseIsolations
{
  /// must be at most MAX_NUM_INTERVENTION_CHANGE_TIMES
  int num_change_times_;

  /// change times for intensity of case isolation
  double change_times_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  /// time-varying equivalent of current_case_isolation.spatial_and_place_effect_
  double spatial_and_place_effects_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  /// time-varying equivalent of current_case_isolation.household_effect_
  double household_effects_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  /// time-varying equivalent of current_case_isolation.prop_
  double prop_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  /// time-varying equivalent of current_case_isolation.cell_inc_thresh_
  double cell_inc_thresh_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];
};

struct HouseholdQuarantine
{
  double spatial_effect_;

  double household_effect_;

  /// indexed by place type;
  double place_effect_[NUM_PLACE_TYPES];

  double individual_prop_comply_;

  double household_prop_comply_;

  double cell_inc_thresh_;
};

struct HouseholdQuarantines
{
  /// must be at most MAX_NUM_INTERVENTION_CHANGE_TIMES
  int num_change_times_;

  /// change times for intensity of household quarantine
  double change_times_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  /// time-varying equivalent of current_household_quarantine_.spatial_effect_
  double spatial_effects_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  /// time-varying equivalent of current_household_quarantine_.household_effect_
  double household_effects_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  /// indexed by i) change time; ii) place type;
  /// time-varying equivalent of current_household_quarantine_.place_effect_
  double place_effects_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES][NUM_PLACE_TYPES];

  /// time-varying equivalent of current_household_quarantine_.individual_prop_comply_
  double individual_prop_comply_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  /// time-varying equivalent of current_household_quarantine_.household_prop_comply_
  double household_prop_comply_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  /// time-varying equivalent of current_household_quarantine_.cell_inc_thresh_
  double cell_inc_thresh_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];
};

struct PlaceClosure
{
  double spatial_effect_;

  double household_effect_;

  /// indexed by place type;
  double place_effect_[NUM_PLACE_TYPES];

  double prop_attending_[NUM_PLACE_TYPES];

  int inc_thresh_;

  double frac_inc_thresh_;

  int cell_inc_thresh_;

  double duration_;
};

struct PlaceClosures
{
  /// must be at most MAX_NUM_INTERVENTION_CHANGE_TIMES
  int num_change_times_;

  /// change times for intensity of place closure
  double change_times_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  /// time-varying equivalent of current_place_closure_.spatial_effect_
  double spatial_effects_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  /// time-varying equivalent of current_place_closure_.household_effect_
  double household_effects_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  /// indexed by i) change time; ii) place type;
  /// time-varying equivalent of current_place_closure_.place_effect_
  double place_effects_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES][NUM_PLACE_TYPES];

  double prop_attending_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES][NUM_PLACE_TYPES];

  /// time-varying equivalent of PlaceCloseIncTrig / current_place_closure_.inc_thresh_
  int inc_thresh_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  /// time-varying equivalent of current_place_closure_.frac_inc_thresh_
  double frac_inc_thresh_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  /// time-varying equivalent of current_place_closure_.cell_inc_thresh_
  int cell_inc_thresh_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  /// time-varying equivalent of current_place_closure_.duration_
  double durs_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];
};

struct DigitalContactTracing
{
  double spatial_and_place_effect_;

  double household_effect_;

  double prop_;

  int max_to_trace_;
};

struct DigitalContactTracings
{
  /// must be at most MAX_NUM_INTERVENTION_CHANGE_TIMES
  int num_change_times_;

  /// change times for intensity of digital contact tracing
  double change_times_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  /// time-varying equivalent of current_digital_contact_tracing_.spatial_and_place_effect_
  double spatial_and_place_effects_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  /// time-varying equivalent of current_digital_contact_tracing_.household_effect_
  double household_effects_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  /// time-varying equivalent of current_digital_contact_tracing_.prop_
  double prop_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];

  int max_to_trace_over_time_[MAX_NUM_INTERVENTION_CHANGE_TIMES];
};

/**
 * @brief Stores the parameters for the simulation.
 *
 */
struct Param
{
	int PopSize; /**< Population size */
	int NH; // Number of households
	int NumRealisations; /**< Number of Realisations */
	int NumNonExtinctRealisations; /**< Number of non-extinct realisations */
	int NRactual;
	int NRactE;
	int NRactNE;
	int UpdatesPerSample; // Number of time steps between samples
	int NumSamples; // Total number of samples that will be made
	int KernelType;
	int NKR; // Size of kernel lookup table
	int NK_HR; // Factor to expand hi-res kernel lookup table by
	int MoveKernelType;
	int AirportKernelType;
	unsigned int BinFileLen;
	int DoBin, DoSaveSnapshot, DoLoadSnapshot;
	double SnapshotSaveTime, SnapshotLoadTime, clP1, clP2, clP3, clP4, clP5, clP6;
	int NC; // Number of cells
	int NMC; // Number of microcells
	int NMCL; // Number of microcells wide/high a cell is; i.e. NMC = NC * NMCL * NMCL
	int NCP; /**< Number of populated cells  */
	int NMCP, ncw, nch, DoUTM_coords, nsp, DoSeasonality, DoCorrectAgeDist, DoPartialImmunity;

	int get_number_of_micro_cells_wide() const;
	int get_number_of_micro_cells_high() const;
	MicroCellPosition get_micro_cell_position_from_cell_index(int cell_index) const;
	int get_micro_cell_index_from_position(MicroCellPosition position) const;
	bool is_in_bounds(MicroCellPosition position) const;

	int DoAdUnits, NumAdunits, DoAdunitBoundaries, AdunitLevel1Divisor, AdunitLevel1Mask, AdunitBitmapDivisor, CountryDivisor;
	int DoAdunitOutput, DoAdunitBoundaryOutput, DoAdunitDemog, DoCorrectAdunitPop, DoSpecifyPop, AdunitLevel1Lookup[ADUNIT_LOOKUP_SIZE];
	int DoOutputPlaceDistForOneAdunit, OutputPlaceDistAdunit, OutputDensFile;
	int DoOneGen, OutputEveryRealisation, BitmapMovieFrame, MaxCorrSample, DoLatent, InfQueuePeakLength, NumThreads, MaxNumThreads;
	int bwidth, bheight; // Size in pixels of the map area in the bitmap output
	int bheight2; // Height in pixels of the entire bitmap output, including both the spectrum at the top and the map area
	int bminx, bminy;
	int OutputBitmap; // Whether to output a bitmap
	BitmapFormats BitmapFormat; // Format of bitmap (platform dependent and command-line /BM: specified).
	int DoSI, DoHeteroDensity, DoPeriodicBoundaries, DoImmuneBitmap, OutputBitmapDetected; //added OutputBitmapDetected - ggilani 04/08/15
	int DoHouseholds, DoPlaces, PlaceTypeNum, Nplace[NUM_PLACE_TYPES], SmallEpidemicCases, DoPlaceGroupTreat;
	int NumInitialInfections[MAX_NUM_SEED_LOCATIONS], DoRandomInitialInfectionLoc, DoAllInitialInfectioninSameLoc;
	int MinPopDensForInitialInfection, NumSeedLocations, MaxPopDensForInitialInfection, InitialInfectionsAdminUnitId[MAX_NUM_SEED_LOCATIONS],InitialInfectionsAdminUnit[MAX_NUM_SEED_LOCATIONS];
	int DoAge, DoSymptoms, LoadSaveNetwork, IncThreshPop, GlobalIncThreshPop;
	int OutputOnlyNonExtinct, DoInfectiousnessProfile, DoInfectionTree, DoWholeHouseholdImmunity, DoSpatial, DoDeath;
	int DoAirports, Nairports, Air_popscale, DoSchoolFile, DoRealSymptWithdrawal, CaseAbsentChildAgeCutoff, DoEarlyCaseDiagnosis, DoInterventionFile;
	int PlaceTypeNoAirNum; // If DoAirports then this is the number of non-airport place types (< PlaceTypeNum), else == PlaceTypeNum (~ no airport places).
	int HotelPlaceType; // If DoAirports then this is place type for hotel (>= PlaceTypeNoAirNum, < PlaceTypeNum), else == PlaceTypeNum (~ unused).
	long setupSeed1, setupSeed2; // RNG seeds from the command line, used to initialise the RNG for setup
	long runSeed1, runSeed2; // RNG seeds from the command line, used to initialise the RNG for running the model
	long nextSetupSeed1, nextSetupSeed2; // The next RNG seeds to use when we need to reinitialise the RNG for setup
	long nextRunSeed1, nextRunSeed2; // The next RNG seeds to use when we need to reinitialise the RNG for the model
	int ResetSeeds,KeepSameSeeds, ResetSeedsPostIntervention, ResetSeedsFlag, TimeToResetSeeds;
	double LongitudeCutLine; // Longitude to image earth is cut at to produce a flat map.  Default -360 degrees (effectively -180).  Use to ensure countries have a contiguous boundary
	double SpatialBoundingBox[4], LocationInitialInfection[MAX_NUM_SEED_LOCATIONS][2], InitialInfectionsAdminUnitWeight[MAX_NUM_SEED_LOCATIONS], TimeStepsPerDay;
	double FalsePositiveRate, FalsePositivePerCapitaIncidence, FalsePositiveAgeRate[NUM_AGE_GROUPS];
	double latent_icdf[CDF_RES + 1], infectious_icdf[CDF_RES + 1], infectious_prof[INFPROF_RES + 1], infectiousness[MAX_INFECTIOUS_STEPS];

	double MildToRecovery_icdf[CDF_RES + 1], ILIToRecovery_icdf[CDF_RES + 1], SARIToRecovery_icdf[CDF_RES + 1], CriticalToCritRecov_icdf[CDF_RES + 1], CritRecovToRecov_icdf[CDF_RES + 1];
	double ILIToSARI_icdf[CDF_RES + 1], SARIToCritical_icdf[CDF_RES + 1], ILIToDeath_icdf[CDF_RES + 1], SARIToDeath_icdf[CDF_RES + 1], CriticalToDeath_icdf[CDF_RES + 1];
	/// means for above icdf's.
	double Mean_MildToRecovery[NUM_AGE_GROUPS], Mean_ILIToRecovery[NUM_AGE_GROUPS], Mean_SARIToRecovery[NUM_AGE_GROUPS], Mean_CriticalToCritRecov[NUM_AGE_GROUPS], Mean_CritRecovToRecov[NUM_AGE_GROUPS];
	double Mean_TimeToTest, Mean_TimeToTestOffset, Mean_TimeToTestCriticalOffset, Mean_TimeToTestCritRecovOffset;
	double Mean_ILIToSARI[NUM_AGE_GROUPS], Mean_SARIToCritical[NUM_AGE_GROUPS], Mean_CriticalToDeath[NUM_AGE_GROUPS], Mean_SARIToDeath[NUM_AGE_GROUPS], Mean_ILIToDeath[NUM_AGE_GROUPS];
	double Prop_Mild_ByAge[NUM_AGE_GROUPS], Prop_ILI_ByAge[NUM_AGE_GROUPS], Prop_SARI_ByAge[NUM_AGE_GROUPS], Prop_Critical_ByAge[NUM_AGE_GROUPS];
	double CFR_SARI_ByAge[NUM_AGE_GROUPS], CFR_Critical_ByAge[NUM_AGE_GROUPS], CFR_ILI_ByAge[NUM_AGE_GROUPS];

	double TimeStep; // The length of a time step, in days
	double SampleTime; // The number of days to run for
	double SampleStep; // The length of a sampling step, in days
	double BitmapAspectScale; // Height of bitmap / Width of bitmap
	int ts_age;
	int DoSeverity; // Non-zero (true) if severity analysis should be done
	double scalex, scaley; // Number of pixels per degree in bitmap output
	DomainSize in_degrees_; ///< Size of spatial domain in degrees
	DomainSize in_cells_; ///< Size of spatial domain in cells
	DomainSize in_microcells_; ///< Size of spatial domain in microcells
	double KernelShape, KernelScale, KernelP3, KernelP4, KernelDelta, MoveKernelShape, MoveKernelScale, MoveKernelP3, MoveKernelP4;
	double AirportKernelShape, AirportKernelScale, AirportKernelP3, AirportKernelP4, AirportTrafficScale;
	double R0, R0scale, LocalBeta;
	double LatentPeriod; // In days. Mean of icdf (inverse cumulative distribution function).
	double InfectiousPeriod; // In days. Mean of icdf (inverse cumulative distribution function).
	double R0household, R0places, R0spatial;
	double Seasonality[DAYS_PER_YEAR];
	double SusceptibilitySD,InfectiousnessSD, R0DensityScalePower;
	double ProportionSymptomatic[NUM_AGE_GROUPS], LatentToSymptDelay, SymptInfectiousness;
	double SymptSpatialContactRate, SymptPlaceTypeContactRate[NUM_PLACE_TYPES], InhibitInterAdunitPlaceAssignment[NUM_PLACE_TYPES];
	double SymptPlaceTypeWithdrawalProp[NUM_PLACE_TYPES], CaseAbsenteeismDuration, CaseAbsenteeismDelay;
	int PlaceCloseRoundHousehold; // Default 1 (close places around a household), 0 (off)
	int AbsenteeismPlaceClosure; // Default 0 (off), 1 (on) track place closures in more detail
	int MaxAbsentTime; // In days.  Max number of days absent, range [0, MAX_ABSENT_TIME].  Default 0 if !P.AbsenteeismPlaceClosure, otherwise MAX_ABSENT_TIME
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
	double PropAgeGroup[MAX_ADUNITS][NUM_AGE_GROUPS], PopByAdunit[MAX_ADUNITS][2];
	double InvLifeExpecDist[MAX_ADUNITS][1001];

	double PlaceCloseTimeStart, PlaceCloseTimeStart2, PlaceCloseDuration, PlaceCloseDuration2, PlaceCloseDelayMean, PlaceCloseRadius, PlaceCloseRadius2;
	double PlaceCloseCasePropThresh, PlaceCloseAdunitPropThresh;
	int DoHolidays, NumHolidays;
	double HolidayEffect[NUM_PLACE_TYPES], HolidayStartTime[DAYS_PER_YEAR], HolidayDuration[DAYS_PER_YEAR];
	double ColourPeriod, BoundingBox[4], BitmapScale;
	double TreatSuscDrop, TreatInfDrop, TreatDeathDrop, TreatSympDrop, TreatDelayMean, TreatTimeStart, TreatPlaceGeogDuration;
	double TreatProphCourseLength, TreatCaseCourseLength, TreatPropRadial, TreatRadius, TreatRadius2, TreatCellIncThresh;
	double DigitalContactTracing_CellIncThresh;
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

	double CaseIsolationDuration;
	double CaseIsolationDelay, CaseIsolationPolicyDuration;

	double HQuarantineTimeStart, HQuarantineDelay, HQuarantineHouseDuration, HQuarantinePolicyDuration;

	int EnhancedSocDistClusterByHousehold;
	double SocDistTimeStart, SocDistDuration, SocDistHouseholdEffect, SocDistPlaceEffect[NUM_PLACE_TYPES], SocDistSpatialEffect;
	double EnhancedSocDistHouseholdEffect, EnhancedSocDistPlaceEffect[NUM_PLACE_TYPES], EnhancedSocDistSpatialEffect, EnhancedSocDistProportionCompliant[NUM_AGE_GROUPS];

	double SocDistChangeDelay, SocDistDuration2;

	double SocDistDurationCurrent;

  SocialDistancing current_social_distancing_;
  SocialDistancing changed_social_distancing_;
  CaseIsolation current_case_isolation_;
  HouseholdQuarantine current_household_quarantine_;
  PlaceClosure current_place_closure_;
  DigitalContactTracing current_digital_contact_tracing_;

	double SocDistRadius, SocDistRadius2;

	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****
	///// **** VARIABLE EFFICACIES OVER TIME
	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****

	int VaryEfficaciesOverTime;

  /// SOCIAL DISTANCING
  SocialDistancings social_distancing_;

  /// CASE ISOLATION
  CaseIsolations case_isolation_;

  /// HOUSEHOLD QUARANTINE
  HouseholdQuarantines household_quarantine_;

  /// PLACE CLOSURE
  PlaceClosures place_closure_;

  /// DIGITAL CONTACT TRACING
  DigitalContactTracings digital_contact_tracing_;

	double KeyWorkerProphTimeStart, KeyWorkerProphDuration, KeyWorkerPropInKeyPlaces[NUM_PLACE_TYPES], KeyWorkerHouseProp;
	double KeyWorkerProphRenewalDuration, KeyWorkerProphRadius, KeyWorkerProphRadius2;

	double TreatTimeStartBase, VaccTimeStartBase, MoveRestrTimeStartBase, PlaceCloseTimeStartBase, PlaceCloseTimeStartBase2,PlaceCloseTimeStartPrevious;
	double AirportCloseTimeStartBase, HQuarantineTimeStartBase, CaseIsolationTimeStartBase, SocDistTimeStartBase, KeyWorkerProphTimeStartBase, DigitalContactTracingTimeStartBase;
	double InfectionImportRate1, InfectionImportRate2, InfectionImportChangeTime, ImportInfectionTimeProfile[MAX_DUR_IMPORT_PROFILE];
	double PreControlClusterIdTime, PreControlClusterIdCalTime, PreControlClusterIdHolOffset, PreIntervIdCalTime,PreIntervTime,SeedingScaling;
	int PreControlClusterIdCaseThreshold, PreControlClusterIdCaseThreshold2, PreControlClusterIdUseDeaths, PreControlClusterIdDuration, DoAlertTriggerAfterInterv, AlertTriggerAfterIntervThreshold,StopCalibration,ModelCalibIteration;
	int DoPerCapitaTriggers, DoGlobalTriggers, DoAdminTriggers, DoICUTriggers, MoveRestrCellIncThresh, DoHQretrigger;

	int PlaceCloseCellIncThresh, PlaceCloseCellIncThresh2, TriggersSamplingInterval, PlaceCloseIndepThresh, VaccPriorityGroupAge[2];
	int PlaceCloseCellIncStopThresh, SocDistCellIncStopThresh;
	int PlaceCloseAdunitPlaceTypes[NUM_PLACE_TYPES];

	int DoPlaceCloseOnceOnly, DoSocDistOnceOnly, DoMoveRestrOnceOnly, DoKeyWorkerProphOnceOnly;

	int VaccMaxRounds, VaccByAdminUnit, VaccAdminUnitDivisor, TreatByAdminUnit, TreatAdminUnitDivisor, MoveRestrByAdminUnit, MoveRestrAdminUnitDivisor, PlaceCloseByAdminUnit, PlaceCloseAdminUnitDivisor;
	int KeyWorkerProphCellIncThresh, KeyWorkerPopNum, KeyWorkerPlaceNum[NUM_PLACE_TYPES], KeyWorkerNum, KeyWorkerIncHouseNum;
	int DoBlanketMoveRestr, PlaceCloseIncTrig, PlaceCloseIncTrig2, TreatMaxCoursesPerCase, DoImportsViaAirports, DoMassVacc, DurImportTimeProfile;
	unsigned short int usHQuarantineHouseDuration, usVaccTimeToEfficacy, usVaccTimeEfficacySwitch; //// us = unsigned short versions of their namesakes, multiplied by P.TimeStepsPerDay
	unsigned short int usCaseIsolationDuration, usCaseIsolationDelay, usCaseAbsenteeismDuration, usCaseAbsenteeismDelay;

	//Added DoRecordInfEvents and MaxInfEvents in order to give the user a choice as to whether to output infection events as a line list: ggilani - 10/10/14
	int DoRecordInfEvents, MaxInfEvents, RecordInfEventsPerRun;
	double KernelPowerScale, KernelOffsetScale;
	int LimitNumInfections, MaxNumInfections;

	//Added parameters to deal with digital contact tracing - ggilani 09/03/2020
	int DoDigitalContactTracing, ClusterDigitalContactUsers, NDigitalContactUsers, NDigitalHouseholdUsers, FindContactsOfDCTContacts, DoDCTTest;
	double PropPopUsingDigitalContactTracing, ScalingFactorSpatialDigitalContacts, ScalingFactorPlaceDigitalContacts, DigitalContactTracingDelay, LengthDigitalContactIsolation, ProportionSmartphoneUsersByAge[NUM_AGE_GROUPS];
	double DelayFromIndexCaseDetectionToDCTIsolation, DelayToTestIndexCase, DelayToTestDCTContacts, SpecificityDCT, SensitivityDCT;
	double DigitalContactTracingPolicyDuration;
	int OutputDigitalContactTracing, OutputDigitalContactDist, DCTIsolateIndexCases, RemoveContactsOfNegativeIndexCase;

	int DoOriginDestinationMatrix; //added: ggilani 28/01/15
	int DoInterventionDelaysByAdUnit;


	int OutputAge, OutputR0, OutputControls, OutputCountry, OutputAdUnitVar, OutputHousehold, OutputInfType, OutputNonSeverity, OutputSeverityAdminUnit, OutputSeverityAge, OutputNonSummaryResults;

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
};

extern Param P;

#endif // COVIDSIM_PARAM_H_INCLUDED_
