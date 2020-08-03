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

	/** \brief  Query whether a host should be included in mass vaccination.
	*           The conditions for not being vaccinated are either: the host is dead,
	*           or the host is a current case, or the host has recovered from
	*           a symptomatic infection.
	*  \return  TRUE if this host should not be vaccinated for the above reasons.
	*/
	bool do_not_vaccinate() const;

	/** \brief  Query whether a host is alive. This includes all states except the
	*           two death states (death through symptomatic, or asymptomatic illness).
	*  \return  TRUE if this host is alive.
	*/
	bool is_alive() const;

	/** \brief  Query whether a host is a current case. (InfStat::Case)
	*   \return TRUE if the host is a current case.
    */
	bool is_case() const;

	/** \brief  Query whether a host is dead. This includes death through either
	*           symptomatic or asymptomatic illess.
	*   \return TRUE if the host is dead
	*/
	bool is_dead() const;

	/** \brief  Query whether a host died with asymptomatic illness. (InfStat::Dead_WasAsymp)
	*   \return TRUE if the host is dead and was asymptomatic.
	*/
	bool is_dead_was_asymp() const;

	/** \brief  Query whether a host died with symptomatic illness. (InfStat::Dead_WasSymp)
	*   \return TRUE if the host is dead and was asymptomatic.
	*/
	bool is_dead_was_symp() const;

	/** \brief  Query whether a host was immune at the start of the epidemic. (InfStat::ImmuneAtStart)
	*   \return TRUE if the host was already immune.
	*/
	bool is_immune_at_start() const;

	/** \brief  Query whether an infected host is asymptomatic, hence not a case. (InfStat::InfectiousAsymptomaticNotCase)
	*   \return TRUE if a host is infectious but asymptomatic.
	*/
	bool is_infectious_asymptomatic_not_case() const;

	/** \brief  Query whether an infected host is about to become symptomatic. (InfStat::InfectiousAlmostSymptomatic)
	*   \return TRUE if the host is about to become symptomatic.
	*/
	bool is_infectious_almost_symptomatic() const;

	/** \brief  Query whether an infected host is in the latent period. (InfStat::Latent)
	*   \return TRUE if the host is in latent period.
	*/
	bool is_latent() const;

	/** \brief  Query whether an infected host is (currently) symptomic, or ever has been. Acceptable states are:
	*           Latent, Infectious Asymptomatic, Infectious Almost Symptomatic, Recovered from Asymptomatic infection,
	*           or died from Asymptomatic infection.
	*   \return TRUE if the host is not, and never was symptomatic.
	*/
	bool is_never_symptomatic() const;

	/** \brief  Query whether an infected host is not yet symptomatic, but could become so. Acceptable states are latent, susceptible, or
	*           Infectious Almost Symptomatic. 
	*   \return TRUE if the host is infected but not yet symptomatic.
	*/
	bool is_not_yet_symptomatic() const;

	/** \brief  Query whether an infected host is a current or potential host for the epidemic. This excludes all who have recovered,
	*           died, or were immune at the start.
	*   \return TRUE if the host is susceptible, or currently infected.
	*/
	bool is_susceptible_or_infected() const;

	/** \brief  Query whether a host has recovered from infection. This includes recovery from both symptomatic and asymptomatic
	*           infections.
	*   \return TRUE if the host has recovered from symptomatic or asymptomatic infection.
	*/
	bool is_recovered() const;

	/** \brief  Query whether a host has recovered from a symptomatic infection.
	*   \return TRUE only if the host has recovered from symptomatic infection.
	*/
	bool is_recovered_symp() const;

	/** \brief  Query whether a host is susceptible.
	*   \return TRUE only if the host is in the susceptible state.
	*/
	bool is_susceptible() const;

	/** \brief  Set host to be a case. (See InfStat::Case)
	*/
	void set_case();

	/** \brief  Set host to have died.
	 *  For infections that were symptomatic, (ie, they passed through the InfStat::Case state),
	 *  the death state is InfStat::DeadWasSymp. Otherwise (ie, they passed through InfStat::InfectiousAsymptomaticNotCase),
	 *  the death state is InfStat::DeadWasAsymp.
	*/
	void set_dead();

	/** \brief  Set host to be initially immune. (See InfStat::ImmuneAtStart)
	*/
	void set_immune_at_start();

	/** \brief  Set host to be infectious, and about to become symptomatic. (See InfStat::InfectiousAlmostSymptomatic)
	*/
	void set_infectious_almost_symptomatic();

	/** \brief  Set host to be infectious, but asymptomatic, so not defined as a case. (See InfStat::InfectiousAsymptomaticNotCase)
	*/
	void set_infectious_asymptomatic_not_case();

	/** \brief  Set host to be in latent stage. (See InfStat::Latent)
	*/
	void set_latent();

	/** \brief  Set host to be recovered.
	*  For infections that were symptomatic, (ie, they passed through the InfStat::Case state),
	*  the recovery state is InfStat::RecoveredFromSymp. Otherwise (ie, they passed through InfStat::InfectiousAsymptomaticNotCase),
	*  the recovery state is InfStat::RecoveredFromAsymp.
	*/
	void set_recovered();

	/** \brief  Set host to be susceptible. (See InfStat::Susceptible)
	*/
	void set_susceptible();

	
	


private:
	InfStat inf;

};

struct PersonQuarantine
{
	uint8_t  comply;		// can be 0, 1, 2
	uint16_t start_time;	// timestep quarantine is started

	// don't remove the extra parentheses around std::numeric_limits<uint16_t>::max
	// because it conflicts with the max() preprocessor macro in MSVC builds
	PersonQuarantine() : comply(2), start_time((std::numeric_limits<uint16_t>::max)()-1) {} 
};
