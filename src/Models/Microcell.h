#ifndef COVIDSIM_MODELS_MICRO_CELL_H_INCLUDED_
#define COVIDSIM_MODELS_MICRO_CELL_H_INCLUDED_

#include "../IndexList.h"


enum struct TreatStat { // treatment status

	//// Untreated
	Untreated = 0,
	//// Untreated but flagged for treatment
	ToBeTreated = 1,
	//// Treated
	Treated = 2,
	//// Do not treat again (flag in TreatSweep there only to avoid code blocks being called again).
	DontTreatAgain = 3,
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
	/* Note use of short int here limits max run time to USHRT_MAX*ModelTimeStep - e.g. 65536*0.25=16384 days=44 yrs.
	   Global search and replace of 'unsigned short int' with 'int' would remove this limit, but use more memory.
	*/

	int n; // Number of people in microcell
	int adunit; // index of admin unit that this microcell belongs to. Will take value \in 1,..., P.NumAdunits (after P.NumAdunits set in SetupModel.cpp::SetupPopulation). NOTE: Distinct from Model.h::AdminUnit::id, which gives identifier of admin unit with reference to population file admin unit file.
	int* members; // array of members/hosts of microcell

	int* places[MAX_NUM_PLACE_TYPES]; // list of places (of various place types) within microcell
	unsigned short int NumPlacesByType[MAX_NUM_PLACE_TYPES]; // number of places (of various place types) within mircocell
	unsigned short int keyworkerproph, move_trig, place_trig, socdist_trig, keyworkerproph_trig;
	unsigned short int move_start_time, move_end_time;
	unsigned short int place_end_time, socdist_end_time, keyworkerproph_end_time;
	TreatStat moverest, treat, vacc, socdist, placeclose;
	unsigned short int treat_trig, vacc_trig;
	unsigned short int treat_start_time, treat_end_time;
	unsigned short int vacc_start_time;
	IndexList* AirportList;
};

#endif
