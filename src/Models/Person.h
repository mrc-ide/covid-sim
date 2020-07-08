#pragma once

#include <cstdint>
#include <limits>

#include "../Country.h"
#include "../InfStat.h"

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
	unsigned int care_home_resident : 1;
	unsigned int quar_comply : 2;		// can be 0, 1, or 2
	unsigned int digitalContactTracingUser : 1;
	unsigned int digitalContactTraced : 1;
	unsigned int index_case_dct : 2;

	unsigned char Travelling;	// Range up to MAX_TRAVEL_TIME
	unsigned char age;
	unsigned char num_treats;		// set to 0 and tested < 2. but never modified?
	unsigned short int PlaceGroupLinks[NUM_PLACE_TYPES];	// These can definitely get > 255

		short int infect_type;		// INFECT_TYPE_MASK
	InfStat inf;
	Severity Severity_Current, Severity_Final; //// Note we allow Severity_Final to take values: Severity_Mild, Severity_ILI, Severity_SARI, Severity_Critical (not e.g. Severity_Dead or Severity_RecoveringFromCritical)

	unsigned short int detected_time; //added hospitalisation flag: ggilani 28/10/2014, added flag to determined whether this person's infection is detected or not
	unsigned short int absent_start_time, absent_stop_time;
	unsigned short int isolation_start_time;
	unsigned short int infection_time, latent_time;		// Set in DoInfect function. infection time is time of infection; latent_time is a misnomer - it is the time at which person become infectious (i.e. infection time + latent period for this person). latent_time will also refer to time of onset with ILI or Mild symptomatic disease.
	unsigned short int recovery_or_death_time;	// set in DoIncub function
	unsigned short int SARI_time, Critical_time, RecoveringFromCritical_time; //// /*mild_time, ILI_time,*/ Time of infectiousness onset same for asymptomatic, Mild, and ILI infection so don't need mild_time etc.
	unsigned short int treat_start_time, treat_stop_time, vacc_start_time;  //// set in TreatSweep function.
	unsigned short int dct_start_time, dct_end_time, dct_trigger_time, dct_test_time; //digital contact tracing start and end time: ggilani 10/03/20
	int ncontacts; //added this in to record total number of contacts each index case records: ggilani 13/04/20

};

struct PersonQuarantine
{
	uint8_t  comply;		// can be 0, 1, 2
	uint16_t start_time;	// timestep quarantine is started

	// don't remove the extra parentheses around std::numeric_limits<uint16_t>::max
	// because it conflicts with the max() preprocessor macro in MSVC builds
	PersonQuarantine() : comply(2), start_time((std::numeric_limits<uint16_t>::max)()-1) {} 
};
