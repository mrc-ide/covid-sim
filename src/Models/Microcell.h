#pragma once

#include "../IndexList.h"

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
