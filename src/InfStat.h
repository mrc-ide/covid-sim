#ifndef COVIDSIM_INFSTAT_H_INCLUDED_
#define COVIDSIM_INFSTAT_H_INCLUDED_

//// Infection Status definitions / labels (generally positive value indicates asymptomatic infection,
//// negative value indicates symptomatic infection).

enum struct InfStat {

//// Note - DO NOT CHANGE these definitions without accounting for "Quarantined not Infected" /
//// "Quarantined not symptomatic" calculation: relies on value below being negative for symptomatic people.

	// Further note: August 2020 - See https://github.com/mrc-ide/covid-sim/pull/445  - which refactors the
	// numerical comparisons or abs() that were used on states. See the functions in Person.cpp. Hence,
	// the values are no longer critical; the old duplicates (Dead and Recovered) are no longer wanted,
	// so any unique values below will work equivalently. For now I will leave them as they are below to allow
	// some historic continuity, as they are "familiar" numbers.

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
	//// Immune at start of epidemic - used to model partially immune population. Distinct from recovered, who recovered during runtime. Doesn't take negative values.
	ImmuneAtStart = 4,
	//// Dead was asymptomatic
	Dead_WasAsymp = 5,
	//// Dead was symptomatic
	Dead_WasSymp = -5,
};

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
