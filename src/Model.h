#ifndef COVIDSIM_MODEL_H_INCLUDED_
#define COVIDSIM_MODEL_H_INCLUDED_

#include "Country.h"
#include "MachineDefines.h"
#include "Constants.h"
#include "InfStat.h"


//// need to test that inequalities in IncubRecoverySweep can be replaced if you initialize to USHRT_MAX, rather than zero.
//// need to output quantities by admin unit

#pragma pack(push, 2)

struct Person
{
	int pcell;			// place cell, Cells[person->pcell] holds this person
	int mcell;			// microcell, Mcells[person->mcell] holds this person
	int hh;				// Household[person->hh] holds this person
	int infector;		// If >=0, Hosts[person->infector] was who infected this person
	int listpos;		// Goes up to at least MAX_SEC_REC, also used as a temp variable?

	int PlaceLinks[NUM_PLACE_TYPES]; //// indexed by i) place type. Value is the number of that place type (e.g. school no. 17; office no. 310 etc.) Place[i][person->PlaceLinks[i]], can be up to P.Nplace[i]
	float infectiousness, susc,ProbAbsent,ProbCare;

	unsigned int esocdist_comply : 1;
	unsigned int keyworker : 1;				// also used to binary index cumI_keyworker[] and related arrays
	unsigned int to_die : 1;
	unsigned int detected : 1; //added hospitalisation flag: ggilani 28/10/2014, added flag to determined whether this person's infection is detected or not

	unsigned char Travelling;	// Range up to MAX_TRAVEL_TIME
	unsigned char age;
	unsigned char quar_comply;		// can be 0, 1, or 2
	unsigned char num_treats;		// set to 0 and tested < 2. but never modified?
	Severity Severity_Current, Severity_Final; //// Note we allow Severity_Final to take values: Severity_Mild, Severity_ILI, Severity_SARI, Severity_Critical (not e.g. Severity_Dead or Severity_RecoveringFromCritical)

	unsigned short int PlaceGroupLinks[NUM_PLACE_TYPES];	// These can definitely get > 255
	short int infect_type;		// INFECT_TYPE_MASK
	InfStat inf;

	unsigned short int detected_time; //added hospitalisation flag: ggilani 28/10/2014, added flag to determined whether this person's infection is detected or not
	unsigned short int absent_start_time, absent_stop_time;
	unsigned short int quar_start_time, isolation_start_time;
	unsigned short int infection_time, latent_time;		// Set in DoInfect function. infection time is time of infection; latent_time is a misnomer - it is the time at which person become infectious (i.e. infection time + latent period for this person). latent_time will also refer to time of onset with ILI or Mild symptomatic disease.
	unsigned short int recovery_or_death_time;	// set in DoIncub function
	unsigned short int SARI_time, Critical_time, RecoveringFromCritical_time; //// /*mild_time, ILI_time,*/ Time of infectiousness onset same for asymptomatic, Mild, and ILI infection so don't need mild_time etc.
	unsigned short int treat_start_time, treat_stop_time, vacc_start_time;  //// set in TreatSweep function.
	unsigned int digitalContactTraced : 1;
	unsigned int index_case_dct : 2;
	unsigned int digitalContactTracingUser : 1;
	unsigned short int dct_start_time, dct_end_time, dct_trigger_time, dct_test_time; //digital contact tracing start and end time: ggilani 10/03/20
	int ncontacts; //added this in to record total number of contacts each index case records: ggilani 13/04/20

};

struct Household
{
	int FirstPerson;
	unsigned short int nh; // number people in household
	float loc_x, loc_y;
	unsigned short int nhr;
};

/*
In the main InfectSweep loop, we cannot safely set
Hosts[infectee].infector and Hosts[infectee].infect_type, as concurrent
threads might be trying to set the values differently. We therefore
make a queue of `infection`s in `inf_queue` containing the information
we need, so that we can set the values after the main loop has finished.
*/
struct Infection
{
	int infector;
	int infectee;
	short int infect_type;
};

/**
 * @brief Contact event used for tracking contact tracing events
 *
 * Currently stores: contact and index case (both ints) and contact time (unsigned short int)
 * Thanks to igfoo
 */
struct ContactEvent
{
	int contact;
	int index;
	unsigned short int contact_time;
};

/**
 * @brief The global state of the model.
 *
 * TODO: Detailed explanation.
 */
struct PopVar
{
	int S, L, I, R, D, cumI, cumR, cumD, cumC, cumTC, cumFC, cumDC, trigDC;
	int cumH; //Added cumulative hospitalisation: ggilani 28/10/14
	int cumCT, cumCC, DCT, cumDCT; //Added total and cumulative contact tracing: ggilani 15/06/17, and equivalents for digital contact tracing: ggilani 11/03/20
	int cumC_country[MAX_COUNTRIES]; //added cumulative cases by country: ggilani 12/11/14
	int cumHQ, cumAC, cumAA, cumAH, cumACS, cumAPC, cumAPA, cumAPCS;
	//// age specific versions of above variables. e.g. cumI is cumulative infections. cumIa is cumulative infections by age group.
	int cumIa[NUM_AGE_GROUPS], cumCa[NUM_AGE_GROUPS], cumDa[NUM_AGE_GROUPS];
	int cumI_adunit[MAX_ADUNITS], cumC_adunit[MAX_ADUNITS], cumD_adunit[MAX_ADUNITS], cumT_adunit[MAX_ADUNITS], cumH_adunit[MAX_ADUNITS], cumDC_adunit[MAX_ADUNITS]; //added cumulative hospitalisation per admin unit: ggilani 28/10/14, cumulative detected cases per adunit: ggilani 03/02/15
	int cumCT_adunit[MAX_ADUNITS], cumCC_adunit[MAX_ADUNITS], trigDC_adunit[MAX_ADUNITS]; //added cumulative CT per admin unit: ggilani 15/06/17
	int cumDCT_adunit[MAX_ADUNITS], DCT_adunit[MAX_ADUNITS]; //added cumulative and overall digital contact tracing per adunit: ggilani 11/03/20
	int cumItype[INFECT_TYPE_MASK], cumI_keyworker[2], cumC_keyworker[2], cumT_keyworker[2];
	Infection *inf_queue[MAX_NUM_THREADS]; // the queue (i.e. list) of infections. 1st index is thread, 2nd is person.
	int n_queue[MAX_NUM_THREADS]; 	// number of infections in inf_queue
	int* p_queue[NUM_PLACE_TYPES], *pg_queue[NUM_PLACE_TYPES], np_queue[NUM_PLACE_TYPES];		// np_queue is number of places in place queue (by place type), p_queue, and pg_queue is the actual place and place-group queue (i.e. list) of places. 1st index is place type, 2nd is place.
	int NumPlacesClosed[NUM_PLACE_TYPES], n_mvacc, mvacc_cum;
	float* cell_inf;  //// List of spatial infectiousnesses by person within cell.
	double sumRad2, maxRad2, cumT, cumV, cumVG, cumUT, cumTP, cumV_daily, cumVG_daily; //added cumVG, cumVG_daily
	int* CellMemberArray, *CellSuscMemberArray;
	int** InvAgeDist;
	int* mvacc_queue;
	int nct_queue[MAX_ADUNITS]; // queue for contact tracing: ggilani 12/06/17
	ContactEvent* dct_queue[MAX_ADUNITS]; //queues for digital contact tracing: ggilani 14/04/20
	int ndct_queue[MAX_ADUNITS]; //queues for digital contact tracing: ggilani 10/03/20
	int contact_dist[MAX_CONTACTS+1]; //added this to store contact distribution: ggilani 13/04/20
	double* origin_dest[MAX_ADUNITS]; //added intermediate storage for calculation of origin-destination matrix: ggilani 02/02/15

	///// Prevalence quantities (+ by admin unit)
	int Mild, ILI, SARI, Critical, CritRecov, /*cumulative incidence*/ cumMild, cumILI, cumSARI, cumCritical, cumCritRecov;
	int Mild_adunit[MAX_ADUNITS], ILI_adunit[MAX_ADUNITS], SARI_adunit[MAX_ADUNITS], Critical_adunit[MAX_ADUNITS], CritRecov_adunit[MAX_ADUNITS];
	/// cum incidence quantities. (+ by admin unit)
	int cumMild_adunit[MAX_ADUNITS], cumILI_adunit[MAX_ADUNITS], cumSARI_adunit[MAX_ADUNITS], cumCritical_adunit[MAX_ADUNITS], cumCritRecov_adunit[MAX_ADUNITS];
	int Mild_age[NUM_AGE_GROUPS], ILI_age[NUM_AGE_GROUPS], SARI_age[NUM_AGE_GROUPS], Critical_age[NUM_AGE_GROUPS], CritRecov_age[NUM_AGE_GROUPS];
	/// cum incidence quantities. (+ by age group)
	int cumMild_age[NUM_AGE_GROUPS], cumILI_age[NUM_AGE_GROUPS], cumSARI_age[NUM_AGE_GROUPS], cumCritical_age[NUM_AGE_GROUPS], cumCritRecov_age[NUM_AGE_GROUPS];

	int cumDeath_ILI, cumDeath_SARI, cumDeath_Critical;		// tracks cumulative deaths from ILI, SARI & Critical severities
	int cumDeath_ILI_adunit[MAX_ADUNITS], cumDeath_SARI_adunit[MAX_ADUNITS], cumDeath_Critical_adunit[MAX_ADUNITS];		// tracks cumulative deaths from ILI, SARI & Critical severities
	int cumDeath_ILI_age[NUM_AGE_GROUPS], cumDeath_SARI_age[NUM_AGE_GROUPS], cumDeath_Critical_age[NUM_AGE_GROUPS]; 

	//// above quantities need to be amended in following parts of code:
	//// i) InitModel (set to zero); Done
	//// ii) RecordSample: (collate from threads);
	//// iii) RecordSample: add to incidence / Timeseries).
	//// iv) SaveResults
	//// v) SaveSummaryResults
	///// need to update these quantities in InitModel (DONE), Record Sample (DONE) (and of course places where you need to increment, decrement).

};

/**
 * @brief Recorded time-series variables (typically populated from the `POPVAR` state)
 *
 * In the function `RecordSample` we transform (copy parts, calculate summary statistics)
 * of the `POPVAR` state into a time-stamped `RESULTS` structure.
 *
 * NOTE: This struct must contain only doubles (and arrays of doubles) for the TSMean
 * 	     averaging code to work.
 */
struct Results
{
	double t, S, L, I, R, D, incC, incTC, incFC, incI, incR, incD, incDC ;
	double incH; //added incidence of hospitalisation: ggilani 28/10/14
	double CT, incCT, incCC, DCT, incDCT; //added total numbers being contact traced and incidence of contact tracing: ggilani 15/06/17, and for digital contact tracing: ggilani 11/03/20
	double incC_country[MAX_COUNTRIES]; //added incidence of cases
	double cumT, cumUT, cumTP, cumV, cumTmax, cumVmax, cumDC, extinct, cumVG; //added cumVG
	double incHQ, incAC, incAH, incAA, incACS, incAPC, incAPA, incAPCS;
	double incIa[NUM_AGE_GROUPS], incCa[NUM_AGE_GROUPS], incDa[NUM_AGE_GROUPS];
	double incItype[INFECT_TYPE_MASK], Rtype[INFECT_TYPE_MASK], Rage[NUM_AGE_GROUPS], Rdenom;
	double rmsRad, maxRad, PropPlacesClosed[NUM_PLACE_TYPES], PropSocDist;
	double incI_adunit[MAX_ADUNITS], incC_adunit[MAX_ADUNITS], cumT_adunit[MAX_ADUNITS], incD_adunit[MAX_ADUNITS], cumD_adunit[MAX_ADUNITS], incH_adunit[MAX_ADUNITS], incDC_adunit[MAX_ADUNITS]; //added incidence of hospitalisation per day: ggilani 28/10/14, incidence of detected cases per adunit,: ggilani 03/02/15
	double incCT_adunit[MAX_ADUNITS], incCC_adunit[MAX_ADUNITS], incDCT_adunit[MAX_ADUNITS], DCT_adunit[MAX_ADUNITS]; //added incidence of contact tracing and number of people being contact traced per admin unit: ggilani 15/06/17
	double incI_keyworker[2], incC_keyworker[2], cumT_keyworker[2];

	///@{
	/** Severity States track the COVID-19 states (e.g., mild, critical, etc.) */
	double Mild, ILI, SARI, Critical, CritRecov;				// Prevalence				//// Must be: i) initialised to zero in SetUpModel. ii) outputted in SaveResults iii) outputted in SaveSummaryResults
	double incMild, incILI, incSARI, incCritical, incCritRecov;	// Incidence				//// Must be: i) initialised to zero in SetUpModel. ii) calculated in RecordSample iii) outputted in SaveResults.
	double cumMild, cumILI, cumSARI, cumCritical, cumCritRecov;	// cumulative incidence		//// Must be: i) initialised to zero in SetUpModel. ii) outputted in SaveResults
	double incDeath_ILI, incDeath_SARI, incDeath_Critical;		// tracks incidence of death from ILI, SARI & Critical severities
	double cumDeath_ILI, cumDeath_SARI, cumDeath_Critical;		// tracks cumulative deaths from ILI, SARI & Critical severities
	///@}

	/////// Severity States by admin unit
	double Mild_adunit[MAX_ADUNITS], ILI_adunit[MAX_ADUNITS], SARI_adunit[MAX_ADUNITS], Critical_adunit[MAX_ADUNITS], CritRecov_adunit[MAX_ADUNITS];				// Prevalence by admin unit
	double incMild_adunit[MAX_ADUNITS], incILI_adunit[MAX_ADUNITS], incSARI_adunit[MAX_ADUNITS], incCritical_adunit[MAX_ADUNITS], incCritRecov_adunit[MAX_ADUNITS];	// incidence by admin unit
	double cumMild_adunit[MAX_ADUNITS], cumILI_adunit[MAX_ADUNITS], cumSARI_adunit[MAX_ADUNITS], cumCritical_adunit[MAX_ADUNITS], cumCritRecov_adunit[MAX_ADUNITS]; // cumulative incidence by admin unit
	double incDeath_ILI_adunit[MAX_ADUNITS], incDeath_SARI_adunit[MAX_ADUNITS], incDeath_Critical_adunit[MAX_ADUNITS];		// tracks incidence of death from ILI, SARI & Critical severities
	double cumDeath_ILI_adunit[MAX_ADUNITS], cumDeath_SARI_adunit[MAX_ADUNITS], cumDeath_Critical_adunit[MAX_ADUNITS];		// tracks cumulative deaths from ILI, SARI & Critical severities

	/////// Severity States by age group
	double Mild_age[NUM_AGE_GROUPS], ILI_age[NUM_AGE_GROUPS], SARI_age[NUM_AGE_GROUPS], Critical_age[NUM_AGE_GROUPS], CritRecov_age[NUM_AGE_GROUPS];				// Prevalence by admin unit
	double incMild_age[NUM_AGE_GROUPS], incILI_age[NUM_AGE_GROUPS], incSARI_age[NUM_AGE_GROUPS], incCritical_age[NUM_AGE_GROUPS], incCritRecov_age[NUM_AGE_GROUPS];	// incidence by admin unit
	double cumMild_age[NUM_AGE_GROUPS], cumILI_age[NUM_AGE_GROUPS], cumSARI_age[NUM_AGE_GROUPS], cumCritical_age[NUM_AGE_GROUPS], cumCritRecov_age[NUM_AGE_GROUPS]; // cumulative incidence by admin unit
	double incDeath_ILI_age[NUM_AGE_GROUPS], incDeath_SARI_age[NUM_AGE_GROUPS], incDeath_Critical_age[NUM_AGE_GROUPS];		// tracks incidence of death from ILI, SARI & Critical severities
	double cumDeath_ILI_age[NUM_AGE_GROUPS], cumDeath_SARI_age[NUM_AGE_GROUPS], cumDeath_Critical_age[NUM_AGE_GROUPS];		// tracks cumulative deaths from ILI, SARI & Critical severities

	/////// possibly need quantities by age (later)
	//// state variables (S, L, I, R) and therefore (Mild, ILI) etc. changed in i) SetUpModel (initialised to zero); ii)

	//// above quantities need to be amended in following parts of code:
	//// i) SetUpModel (set to zero);
	//// ii) RecordSample: add to incidence / Timeseries).
	//// iii) SaveResults and SaveSummary results.
	///// need to update these quantities in InitModel (DONE), Record Sample (DONE) (and of course places where you need to increment, decrement).

};

/**
 * Supports producing individual infection events from the simulation (and is not used that
 * much because it was developed for Ebola, and slows the simulation).
 *
 * Added Events struct to allow us to log and write out infection events: ggilani 10/10/14
 */
struct Events
{
	double infectee_x, infectee_y, t, t_infector;
	int run, infectee_ind, infector_ind, type, infectee_adunit, listpos, infectee_cell, infector_cell, thread;
};

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
struct IndexList
{
	int id;
	float prob;
};

/**
 * @brief Airport state.
 *
 * Not used for COVID-19 right now. Might be more relevant for USA and
 * other countries that have lots of internal flights. Slows the simulation.
 */
struct Airport
{
	int num_mcell, num_place, Inv_prop_traffic[129], Inv_DestMcells[1025], Inv_DestPlaces[1025];
	unsigned short int num_connected, *conn_airports;
	float total_traffic, loc_x, loc_y;
	float* prop_traffic;
	IndexList* DestMcells, *DestPlaces;
};

/**
 * @brief The basic unit of the simulation and is associated to a geographical location.
 *
 * Interventions (e.g., school closures) are tracked at this level. It contains a list of its
 * members (people), places (schools, universities, workplaces etc.), road networks, links to
 * airports etc.
 */
struct Microcell
{
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
	IndexList* AirportList;
};

/**
 * @brief Holds microcells.
 *
 * Keeps track of susceptible, latent and infected people (in addition to details like who
 * is vaccinated, treated etc.) Also contains data for the spatial gravity model for social
 * interactions (probability distributions).
*/
struct Cell
{
	int n, S, L, I, R, D, cumTC, S0, tot_treat, tot_vacc;
	int* members, *susceptible, *latent, *infected; //// pointers to people in cell. e.g. *susceptible identifies where the final susceptible member of cel is.
	int* InvCDF;
	float tot_prob, *cum_trans, *max_trans;
	short int CurInterv[MAX_INTERVENTION_TYPES];
};

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
struct Place
{
	int n, mcell;
	unsigned short int ng, treat, control_trig, country;
	unsigned short int close_start_time, close_end_time, treat_end_time;
	unsigned short int* AvailByAge;
	unsigned short int Absent[MAX_ABSENT_TIME], AbsentLastUpdateTime;
	float loc_x, loc_y;
	float ProbClose;
	int* group_start, *group_size, *members;
};

/**
 * @brief Deprecated intervention mechanism.
 *
 * Not currently being used, but may be reinstated.
 */
struct Intervention
{
	int InterventionType, DoAUThresh, NoStartAfterMin;
	double StartTime, StopTime, MinDuration, RepeatInterval, TimeOffset;
	double StartThresholdHigh, StartThresholdLow, StopThreshold, Level, LevelCellVar, LevelAUVar, LevelCountryVar, ControlParam, LevelClustering;
	unsigned int MaxRounds, MaxResource;
};

/**
 * @brief A political entity that administers a geographical area.
 */
struct AdminUnit
{
	int id, cnt_id, NI, n; //added n - number of people in admin unit: ggilani 05/01/15
	Intervention InterventionList[MAX_INTERVENTIONS_PER_ADUNIT];
	char cnt_name[100], ad_name[200];
	int NP, place_close_trig;
	double CaseIsolationTimeStart, HQuarantineTimeStart, DigitalContactTracingTimeStart;
	double SocialDistanceTimeStart, PlaceCloseTimeStart; //added these to admin unit in the hope of getting specific start times for Italy: ggilani 16/03/20
	//adding in admin level delays and durations for admin units: ggilani 17/03/20
	double SocialDistanceDelay, HQuarantineDelay, CaseIsolationDelay, PlaceCloseDelay, DCTDelay;
	double SocialDistanceDuration, HQuarantineDuration, CaseIsolationPolicyDuration, PlaceCloseDuration, DCTDuration;
	int* dct, ndct; //arrays for admin unit based digital contact tracing: ggilani 10/03/20
	double* origin_dest; //storage for origin-destination matrix between admin units: ggilani 28/01/15
};

#pragma pack(pop)

extern Person* Hosts;
extern Household* Households;
extern PopVar State, StateT[MAX_NUM_THREADS];
extern Cell* Cells, ** CellLookup;
extern Microcell* Mcells, ** McellLookup;
extern Place** Places;
extern AdminUnit AdUnits[MAX_ADUNITS];

//// Time Series defs:
//// TimeSeries is an array of type results, used to store (unsurprisingly) a time series of every quantity in results. Mostly used in RecordSample.
//// TSMeanNE and TSVarNE are the mean and variance of non-extinct time series. TSMeanE and TSVarE are the mean and variance of extinct time series. TSMean and TSVar are pointers that point to either extinct or non-extinct.
extern Results* TimeSeries, *TSMean, *TSVar, *TSMeanNE, *TSVarNE, *TSMeanE, *TSVarE; //// TimeSeries used in RecordSample, RecordInfTypes, SaveResults. TSMean and TSVar

extern Airport* Airports;
extern Events* InfEventLog;
extern int* nEvents;


extern double inftype[INFECT_TYPE_MASK], inftype_av[INFECT_TYPE_MASK], infcountry[MAX_COUNTRIES], infcountry_av[MAX_COUNTRIES], infcountry_num[MAX_COUNTRIES];
extern double indivR0[MAX_SEC_REC][MAX_GEN_REC], indivR0_av[MAX_SEC_REC][MAX_GEN_REC];
extern double inf_household[MAX_HOUSEHOLD_SIZE + 1][MAX_HOUSEHOLD_SIZE + 1], denom_household[MAX_HOUSEHOLD_SIZE + 1];
extern double inf_household_av[MAX_HOUSEHOLD_SIZE + 1][MAX_HOUSEHOLD_SIZE + 1], AgeDist[NUM_AGE_GROUPS], AgeDist2[NUM_AGE_GROUPS];
extern double case_household[MAX_HOUSEHOLD_SIZE + 1][MAX_HOUSEHOLD_SIZE + 1], case_household_av[MAX_HOUSEHOLD_SIZE + 1][MAX_HOUSEHOLD_SIZE + 1];
extern double PropPlaces[NUM_AGE_GROUPS * AGE_GROUP_WIDTH][NUM_PLACE_TYPES];
extern double PropPlacesC[NUM_AGE_GROUPS * AGE_GROUP_WIDTH][NUM_PLACE_TYPES], AirTravelDist[MAX_DIST];
extern double PeakHeightSum, PeakHeightSS, PeakTimeSum, PeakTimeSS;

extern int DoInitUpdateProbs;

#endif // COVIDSIM_MODEL_H_INCLUDED_
