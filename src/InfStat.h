#ifndef COVIDSIM_INFSTAT_H_INCLUDED_
#define COVIDSIM_INFSTAT_H_INCLUDED_

//// Infection Status definitions / labels (generally positive value indicates asymptomatic infection,
//// negative value indicates symptomatic infection).

enum struct InfStat {

//// Note - DO NOT CHANGE these definitions without accounting for "Quarantined not Infected" /
//// "Quarantined not symptomatic" calculation: relies on value below being negative for symptomatic people.

	//// Susceptible
	Susceptible = 0,
	//// Neither infectious nor symptomatic (E or L).
	Latent = 1,
	//// Infectious and about to become a case.
	InfectiousAlmostSymptomatic = -1,
	//// Not just asymptomatic, but also will not become symptomatic (i.e. a case)
	InfectiousAsymptomaticNotCase = 2,
	//// Case. Infectious and symptomatic.
	Case = -2,
	//// Recovered from asymptomatic infection
	RecoveredFromAsymp = 3,
	//// Recovered from symptomatic infection
	RecoveredFromSymp = -3,
	//// InfStat_Recovered (will use this for abs() values) so code reads correctly
	Recovered = 3,
	//// Immune at start of epidemic - used to model partially immune population. Distinct from recovered, who recovered during runtime. Doesn't take negative values.
	ImmuneAtStart = 4,
	//// Dead was asymptomatic
	Dead_WasAsymp = 5,
	//// Dead was symptomatic
	Dead_WasSymp = -5,
	//// Dead (will use this for abs() values) so code reads correctly
	Dead = 5
};

bool not_yet_symptomatic(InfStat x);
bool is_dead(InfStat x);
bool is_recovered(InfStat x);
bool pre_recovered(InfStat x);
bool never_symptomatic(InfStat x);

//// SeverityClass definitions / labels (numbers arbitrary but must be different to each other).
enum struct Severity {
	//// Flag value.
	Asymptomatic,
	Mild,
	ILI,
	SARI,
	Critical,
	//// Recovering from Critical (not recovered yet).
	RecoveringFromCritical,
	//// label to avoid double counting. Not sure you need.
	Dead,
	//// label to avoid double counting. Not sure you need.
	Recovered
};

#endif // COVIDSIM_INFSTAT_H_INCLUDED_
