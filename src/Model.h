#pragma once

#ifndef SPATIALSIM_MODEL_H_INCLUDED_
#define SPATIALSIM_MODEL_H_INCLUDED_

#include "Country.h"
#include "MachineDefines.h"
#include "Constants.h"


//// need to test that inequalities in IncubRecoverySweep can be replaced if you initialize to USHRT_MAX, rather than zero. 
//// need to output quantities by admin unit

#pragma pack(push, 2)

typedef struct PERSON {

	int pcell; //// place cell, 
	int mcell; //// microcell
	int listpos, infector, hh /*= household*/;
	unsigned short int age;
	unsigned short int to_die, Travelling, detected; //added hospitalisation flag: ggilani 28/10/2014, added flag to determined whether this person's infection is detected or not
	int digitalContactTraced;
	unsigned short int esocdist_comply, quar_start_time, isolation_start_time, isolation_stop_time, keyworker; //added isolation stop time to help with contact tracing
	unsigned char quar_comply;
	unsigned short int absent_start_time, absent_stop_time;
	unsigned short int PlaceGroupLinks[NUM_PLACE_TYPES];
	int PlaceLinks[NUM_PLACE_TYPES]; //// indexed by i) place type. Value is the number of that place type (e.g. school no. 17; office no. 310 etc.)
	float infectiousness, susc, vacc_eff; 
	short int inf, infect_type;
	unsigned short int infection_time, latent_time;		// Set in DoInfect function. infection time is time of infection; latent_time is time at which you become infectious (i.e. infection time + latent period for this person). latent_time will also refer to time of onset with ILI or Mild symptomatic disease. 
	unsigned short int recovery_time;	// set in DoIncub function (note recovery_time can be death_time also) 
	unsigned short int treat_start_time, treat_stop_time, vacc_start_time;  //// set in TreatSweep function. 
	unsigned short int dct_start_time, dct_end_time; //digital contact tracing start and end time: ggilani 10/03/20
	short int num_treats;
	char Severity_Current, Severity_Final; //// Note we allow Severity_Final to take values: Severity_Mild, Severity_ILI, Severity_SARI, Severity_Critical (not e.g. Severity_Dead or Severity_RecoveringFromCritical)
	unsigned short int SARI_time, Critical_time, RecoveringFromCritical_time; //// /*mild_time, ILI_time,*/ Time of infectiousness onset same for asymptomatic, Mild, and ILI infection so don't need mild_time etc. 

	//added variables for digital contact tracing
	int digitalContactTracingUser, *contacts;
	short ncontacts;

} person;

typedef struct HOUSEHOLD {
	int FirstPerson;
	unsigned short int nh; // number people in household
	float loc_x, loc_y;
	unsigned short int nhr;
} household;

/**
 * @brief The global state of the model.
 *
 * TODO: Detailed explanation.
 */
typedef struct POPVAR {

	int NL, S, L, I, R, D, cumI, cumR, cumD, cumC, cumTC, cumFC, cumInf_h, cumInf_n, cumInf_s, cumDC, trigDC;
	int H, cumH; //Added total and cumulative hospitalisation: ggilani 28/10/14
	int CT, cumCT, CC, cumCC, DCT, cumDCT; //Added total and cumulative contact tracing: ggilani 15/06/17, and equivalents for digital contact tracing: ggilani 11/03/20
	int cumC_country[MAX_COUNTRIES]; //added cumulative cases by country: ggilani 12/11/14
	int cumHQ, cumAC, cumAA, cumAH, cumACS, cumAPC, cumAPA, cumAPCS;
	//// age specific versions of above variables. e.g. cumI is cumulative infections. cumIa is cumulative infections by age group. 
	int Na[NUM_AGE_GROUPS], cumIa[NUM_AGE_GROUPS], cumCa[NUM_AGE_GROUPS], cumDa[NUM_AGE_GROUPS];
	int cumI_adunit[MAX_ADUNITS], cumC_adunit[MAX_ADUNITS], cumD_adunit[MAX_ADUNITS], cumT_adunit[MAX_ADUNITS], cumH_adunit[MAX_ADUNITS], H_adunit[MAX_ADUNITS], cumDC_adunit[MAX_ADUNITS]; //added cumulative hospitalisation per admin unit: ggilani 28/10/14, cumulative detected cases per adunit: ggilani 03/02/15
	int cumCT_adunit[MAX_ADUNITS], CT_adunit[MAX_ADUNITS], cumCC_adunit[MAX_ADUNITS], CC_adunit[MAX_ADUNITS], trigDC_adunit[MAX_ADUNITS]; //added cumulative and CT per admin unit: ggilani 15/06/17
	int cumDCT_adunit[MAX_ADUNITS], DCT_adunit[MAX_ADUNITS]; //added cumulative and overall digital contact tracing per adunit: ggilani 11/03/20
	int cumItype[INFECT_TYPE_MASK], cumI_keyworker[2], cumC_keyworker[2], cumT_keyworker[2];
	int* inv_cell_inf; //// think indexed by i) person. 
	int *inf_queue[MAX_NUM_THREADS], n_queue[MAX_NUM_THREADS * CACHE_LINE_SIZE]; 	// n_queue is number of people in the queue, inf_queue is the actual queue (i.e. list) of people. 1st index is thread, 2nd is person.
	int* p_queue[NUM_PLACE_TYPES], *pg_queue[NUM_PLACE_TYPES], np_queue[NUM_PLACE_TYPES];		// np_queue is number of places in place queue (by place type), p_queue, and pg_queue is the actual place and place-group queue (i.e. list) of places. 1st index is place type, 2nd is place.
	int NumPlacesClosed[NUM_PLACE_TYPES], n_mvacc, mvacc_cum;
	float* cell_inf;  //// List of spatial infectiousnesses by person within cell. 
	double sumRad2, maxRad2, cumT, cumV, cumVG, cumUT, cumTP, cumV_daily, cumVG_daily; //added cumVG, cumVG_daily
	int* CellMemberArray, *CellSuscMemberArray;
	int** InvAgeDist;
	int* mvacc_queue;
	int dum[CACHE_LINE_SIZE];
	int* h_queue[MAX_ADUNITS], nh_queue[MAX_ADUNITS], *hd_queue[MAX_ADUNITS], nhd_queue[MAX_ADUNITS]; //queues for hospitalisation by admin unit: ggilani 30/10/14. h_queue and hd_queue actually 2D. Weirdly one dimension allocated on stack and other is pointer. d refers to discharge. n to length of queue. 
	int* ct_queue[MAX_ADUNITS], nct_queue[MAX_ADUNITS]; // queues for contact tracing: ggilani 12/06/17
	int* dct_queue[MAX_ADUNITS], ndct_queue[MAX_ADUNITS]; //queues for digital contact tracing: ggilani 10/03/20
	double* origin_dest[MAX_ADUNITS]; //added intermediate storage for calculation of origin-destination matrix: ggilani 02/02/15

									  ///// Prevalence quantities (+ by admin unit)
	int Mild, ILI, SARI, Critical, CritRecov, /*cumulative incidence*/ cumMild, cumILI, cumSARI, cumCritical, cumCritRecov;
	int Mild_adunit[MAX_ADUNITS], ILI_adunit[MAX_ADUNITS], SARI_adunit[MAX_ADUNITS], Critical_adunit[MAX_ADUNITS], CritRecov_adunit[MAX_ADUNITS];
	/// cum incidence quantities. (+ by admin unit)
	int cumMild_adunit[MAX_ADUNITS], cumILI_adunit[MAX_ADUNITS], cumSARI_adunit[MAX_ADUNITS], cumCritical_adunit[MAX_ADUNITS], cumCritRecov_adunit[MAX_ADUNITS];

	//// above quantities need to be amended in following parts of code: 
	//// i) InitModel (set to zero); Done
	//// ii) RecordSample: (collate from threads); 
	//// iii) RecordSample: add to incidence / Timeseries). 
	//// iv) Print out statement in RunModel but you don't need to add to that yet. 
	///// need to update these quantities in InitModel (DONE), Record Sample (DONE) (and of course places where you need to increment, decrement). 



} popvar;

/**
 * @brief Recorded time-series variables (typically populated from the `POPVAR` state) 
 * 
 * In the function `RecordSample` we transform (copy parts, calculate summary statistics) 
 * of the `POPVAR` state into a time-stamped `RESULTS` structure.
 */
typedef struct RESULTS {

	double t, S, L, I, R, D, incC, incTC, incFC, incL, incI, incR, incD, incDC; 
	double H, incH; //added total hospitalisation and incidence of hospitalisation: ggilani 28/10/14
	double CT, incCT, CC, incCC, DCT, incDCT; //added total numbers being contact traced and incidence of contact tracing: ggilani 15/06/17, and for digital contact tracing: ggilani 11/03/20
	double incC_country[MAX_COUNTRIES]; //added incidence of cases
	double cumT, cumUT, cumTP, cumV, cumTmax, cumVmax, cumDC, extinct, cumVG; //added cumVG
	double incHQ, incAC, incAH, incAA, incACS, incAPC, incAPA, incAPCS;
	double incIa[NUM_AGE_GROUPS], incCa[NUM_AGE_GROUPS], incDa[NUM_AGE_GROUPS];
	double incItype[INFECT_TYPE_MASK], Rtype[INFECT_TYPE_MASK], Rage[NUM_AGE_GROUPS], Rdenom;
	double rmsRad, maxRad, PropPlacesClosed[NUM_PLACE_TYPES], PropSocDist;
	double incI_adunit[MAX_ADUNITS], incC_adunit[MAX_ADUNITS], cumT_adunit[MAX_ADUNITS], incD_adunit[MAX_ADUNITS], cumD_adunit[MAX_ADUNITS], incH_adunit[MAX_ADUNITS], H_adunit[MAX_ADUNITS], incDC_adunit[MAX_ADUNITS]; //added incidence of hospitalisation per day: ggilani 28/10/14, incidence of detected cases per adunit,: ggilani 03/02/15
	double incCT_adunit[MAX_ADUNITS], CT_adunit[MAX_ADUNITS], incCC_adunit[MAX_ADUNITS], CC_adunit[MAX_ADUNITS], incDCT_adunit[MAX_ADUNITS], DCT_adunit[MAX_ADUNITS]; //added incidence of contact tracing and number of people being contact traced per admin unit: ggilani 15/06/17
	double incI_keyworker[2], incC_keyworker[2], cumT_keyworker[2];
	float* bmi2, *bmi3, *bmi4;

	///@{
	/** Severity States track the COVID-19 states (e.g., mild, critical, etc.) */
	double Mild, ILI, SARI, Critical, CritRecov;				// Prevalence				//// Must be: i) initialised to zero in SetUpModel. ii) outputted in SaveResults iii) outputted in SaveSummaryResults
	double incMild, incILI, incSARI, incCritical, incCritRecov;	// Incidence				//// Must be: i) initialised to zero in SetUpModel. ii) calculated in RecordSample iii) outputted in SaveResults. 
	double cumMild, cumILI, cumSARI, cumCritical, cumCritRecov;	// cumulative incidence		//// Must be: i) initialised to zero in SetUpModel. ii) outputted in SaveResults
	///@}

	/////// Severity States by admin unit
	double Mild_adunit[MAX_ADUNITS], ILI_adunit[MAX_ADUNITS], SARI_adunit[MAX_ADUNITS], Critical_adunit[MAX_ADUNITS], CritRecov_adunit[MAX_ADUNITS];				// Prevalence by admin unit
	double incMild_adunit[MAX_ADUNITS], incILI_adunit[MAX_ADUNITS], incSARI_adunit[MAX_ADUNITS], incCritical_adunit[MAX_ADUNITS], incCritRecov_adunit[MAX_ADUNITS];	// incidence by admin unit
	double cumMild_adunit[MAX_ADUNITS], cumILI_adunit[MAX_ADUNITS], cumSARI_adunit[MAX_ADUNITS], cumCritical_adunit[MAX_ADUNITS], cumCritRecov_adunit[MAX_ADUNITS]; // cumulative incidence by admin unit

	/////// possibly need quantities by age (later)
	//// state varialbes (S, L, I, R) and therefore (Mild, ILI) etc. changed in i) SetUpModel (initialised to zero); ii) 

	//// above quantities need to be amended in following parts of code: 
	//// i) InitModel (set to zero); Done
	//// ii) RecordSample: (collate from threads); 
	//// iii) RecordSample: add to incidence / Timeseries). 
	//// iv) Print out statement in RunModel but you don't need to add to that yet. 
	///// SaveResults and SaveSummary results. 
	///// need to update these quantities in InitModel (DONE), Record Sample (DONE) (and of course places where you need to increment, decrement). 



} results;

/**
 * Supports producing individual infection events from the simulation (and is not used that 
 * much because it was developed for Ebola, and slows the simulation).
 * 
 * Added Events struct to allow us to log and write out infection events: ggilani 10/10/14
 */ 
typedef struct EVENTS {
	double infectee_x, infectee_y, t, t_infector;
	int run, infectee_ind, infector_ind, type, infectee_adunit, listpos, infectee_cell, infector_cell, thread;
} events;

/*
  HQ - quarantined households
  AH - Quarantined (and perhaps sick) working adults
  AC - Non-quarantined working adult cases absent thru sickness.
  AA - Absent working adults who are caring for sick children (only assigned if no non-working, quarantine or sick adults available).
  ACS - Children below care requirement cut-off age who are absent (sick or quarantined)

  APC - Non-quarantined working adult cases absent due to closure of their workplace (excl teachers).
  APA - Absent working adults who are caring for children at home due to school closure (only assigned if no non-working, quarantine or sick adults available).
  ACS - Children below care requirement cut-off age who are absent due to school closure


  AH x rq + AC + AA +(APC+APA) x rc = total adult absence
  rq=ratio of quarantine time to duration of absence due to illness
  rc=ratio of school/workplace closure duration of absence due to illness
*/

/**
 * @brief Used for computing spatial interactions more efficiently.
 */ 
typedef struct INDEXLIST {
	int id;
	float prob;
} indexlist;

/**
 * @brief Airport state.
 * 
 * Not used for COVID-19 right now. Might be more relevant for USA and 
 * other countries that have lots of internal flights. Slows the simulation.
 */ 
typedef struct AIRPORT {
	int num_mcell, num_place, Inv_prop_traffic[129], Inv_DestMcells[1025], Inv_DestPlaces[1025];
	unsigned short int country, adunit, num_connected, control, *conn_airports;
	float int_traffic, total_traffic, loc_x, loc_y, control_start_time, control_end_time;
	float* prop_traffic;
	indexlist* DestMcells, *DestPlaces;
} airport;

/**
 * @brief The basic unit of the simulation and is associated to a geographical location. 
 * 
 * Interventions (e.g., school closures) are tracked at this level. It contains a list of its 
 * members (people), places (schools, universities, workplaces etc.), road networks, links to 
 * airports etc.
 */ 
typedef struct MICROCELL {
	/* Note use of short int here limits max run time to USHRT_MAX*TimeStep - e.g. 65536*0.25=16384 days=44 yrs.
	   Global search and replace of 'unsigned short int' with 'int' would remove this limit, but use more memory.
	*/
	int n /*Number of people in microcell*/, adunit;
	int* members;
	unsigned short int country;

	int* places[NUM_PLACE_TYPES];
	unsigned short int np[NUM_PLACE_TYPES];
	unsigned short int moverest, placeclose, socdist, keyworkerproph, move_trig, place_trig, socdist_trig, keyworkerproph_trig;
	unsigned short int move_start_time, move_end_time;
	unsigned short int place_end_time, socdist_end_time, keyworkerproph_end_time;
	unsigned short int treat, vacc, treat_trig, vacc_trig;
	unsigned short int treat_start_time, treat_end_time;
	unsigned short int vacc_start_time;
	indexlist* AirportList;
} microcell;

/**
 * @brief Holds microcells.
 * 
 * Keeps track of susceptible, latent and infected people (in addition to details like who 
 * is vaccinated, treated etc.) Also contains data for the spatial gravity model for social 
 * interactions (probability distributions).
*/
typedef struct CELL {
	int n, S, L, I, R, D, cumTC, S0, tot_treat, tot_vacc;
	int* members, *susceptible, *latent, *infected; //// pointers to people in cell. e.g. *susceptible identifies where the final susceptible member of cel is. 
	int* InvCDF;
	float tot_prob, *cum_trans, *max_trans;
	short int CurInterv[MAX_INTERVENTION_TYPES];
} cell;

/**
 * @brief Represents an institution that people may belong to.
 * 
 * PLACE be an elementary school, high schools, universities, workplaces etc. Places 
 * belong to a microcell (and therefore have a spatial location). Places may have state 
 * (i.e., closed or open). Mechanisms exist for absenteeism tracking (but are not currently used). 
 * The `members` array lists all individuals who belong to a place.
 * Places can have different groups (to model differential interaction strengths between groups 
 * in the same place).
 */
typedef struct PLACE {
	int n, mcell;
	unsigned short int ng, treat, control_trig, country;
	unsigned short int close_start_time, close_end_time, treat_end_time;
	unsigned short int* AvailByAge;
#ifdef ABSENTEEISM_PLACE_CLOSURE
	unsigned short int Absent[MAX_ABSENT_TIME], AbsentLastUpdateTime;
#endif
	float loc_x, loc_y;
	int* group_start, *group_size, *members;
} place;

/**
 * @brief Deprecated intervention mechanism.
 * 
 * Not currently being used, but may be reinstated.
 */ 
typedef struct INTERVENTION {
	int InterventionType, DoAUThresh, NoStartAfterMin;
	double StartTime, StopTime, MinDuration, RepeatInterval, TimeOffset;
	double StartThresholdHigh, StartThresholdLow, StopThreshold, Level, LevelCellVar, LevelAUVar, LevelCountryVar, ControlParam, LevelClustering;
	unsigned int MaxRounds, MaxResource;
} intervention;

/**
 * @brief A political entity that administers a geographical area.
 */
typedef struct ADMINUNIT {
	int id, cnt_id, NI, NIactive[MAX_INTERVENTION_TYPES], n; //added n - number of people in admin unit: ggilani 05/01/15
	int ActiveInterv[MAX_INTERVENTION_TYPES][MAX_INTERVENTIONS_PER_ADUNIT];
	intervention InterventionList[MAX_INTERVENTIONS_PER_ADUNIT];
	char cnt_name[100], ad_name[200];
	int NP, place_close_trig;
	double CaseIsolationTimeStart, HQuarantineTimeStart, DigitalContactTracingTimeStart;
	double SocialDistanceTimeStart, PlaceCloseTimeStart; //added these to admin unit in the hope of getting specific start times for Italy: ggilani 16/03/20
	//adding in admin level delays and durations for admin units: ggilani 17/03/20
	double SocialDistanceDelay, HQuarantineDelay, CaseIsolationDelay, PlaceCloseDelay;
	double SocialDistanceDuration, HQuarantineDuration, CaseIsolationDuration, PlaceCloseDuration;
	int* dct_queue, ndct_queue, *dct, ndct; //arrays for admin unit based digital contact tracing: ggilani 10/03/20
	double* origin_dest; //storage for origin-destination matrix between admin units: ggilani 28/01/15
	double caseDetectRate; //case detection rate: ggilani 03/02/15
	double caseDetectInit;
} adminunit;

#pragma pack(pop)

extern person* Hosts;
extern household* Households;
extern popvar State, StateT[MAX_NUM_THREADS];
extern cell* Cells, ** CellLookup;
extern microcell* Mcells, ** McellLookup;
extern place** Places;
extern adminunit AdUnits[MAX_ADUNITS];

//// Time Series defs: 
//// TimeSeries is an array of type results, used to store (unsurprisingly) a time series of every quantity in results. Mostly used in RecordSample.
//// TSMeanNE and TSVarNE are the mean and variance of non-extinct time series. TSMeanE and TSVarE are the mean and variance of extinct time series. TSMean and TSVar are pointers that point to either extinct or non-extinct. 
extern results* TimeSeries, *TSMean, *TSVar, *TSMeanNE, *TSVarNE, *TSMeanE, *TSVarE; //// TimeSeries used in RecordSample, RecordInfTypes, SaveResults. TSMean and TSVar

extern airport* Airports;
extern events* InfEventLog;
extern int* nEvents;
extern int* RevCellLookup;


extern double ** PopDensity, *mcell_dens;
extern int* mcell_adunits, *mcell_num, *mcell_country;
extern double inftype[INFECT_TYPE_MASK], inftype_av[INFECT_TYPE_MASK], infcountry[MAX_COUNTRIES], infcountry_av[MAX_COUNTRIES], infcountry_num[MAX_COUNTRIES];
extern double indivR0[MAX_SEC_REC][MAX_GEN_REC], indivR0_av[MAX_SEC_REC][MAX_GEN_REC];
extern double inf_household[MAX_HOUSEHOLD_SIZE + 1][MAX_HOUSEHOLD_SIZE + 1], denom_household[MAX_HOUSEHOLD_SIZE + 1];
extern double inf_household_av[MAX_HOUSEHOLD_SIZE + 1][MAX_HOUSEHOLD_SIZE + 1], AgeDist[NUM_AGE_GROUPS], AgeDist2[NUM_AGE_GROUPS];
extern double case_household[MAX_HOUSEHOLD_SIZE + 1][MAX_HOUSEHOLD_SIZE + 1], case_household_av[MAX_HOUSEHOLD_SIZE + 1][MAX_HOUSEHOLD_SIZE + 1];
extern double PropPlaces[NUM_AGE_GROUPS * AGE_GROUP_WIDTH][NUM_PLACE_TYPES];
extern double PropPlacesC[NUM_AGE_GROUPS * AGE_GROUP_WIDTH][NUM_PLACE_TYPES], AirTravelDist[MAX_DIST];
extern double PeakHeightSum, PeakHeightSS, PeakTimeSum, PeakTimeSS;

extern int DoInitUpdateProbs;

#endif // SPATIALSIM_MODEL_H_INCLUDED_
