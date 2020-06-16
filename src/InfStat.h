#ifndef COVIDSIM_INFSTAT_H_INCLUDED_
#define COVIDSIM_INFSTAT_H_INCLUDED_

//// Infection Status definitions / labels (generally positive value indicates asymptomatic infection, 
//// negative value indicates symptomatic infection).

enum InfStat {	
  
//// Note - DO NOT CHANGE these definitions without accounting for "Quarantined not Infected" / 
//// "Quarantined not symptomatic" calculation: relies on value below being negative for symptomatic people.

	//// Susceptible
	InfStat_Susceptible = 0,
	//// E or L (neither infectious nor symptomatic).
	InfStat_Latent = 1,
	//// Infectious and about to become a case.
	InfStat_InfectiousAlmostSymptomatic = -1,
	//// Not just asymptomatic, but also will not become symptomatic (i.e. a case.)
	InfStat_InfectiousAsymptomaticNotCase = 2,
	//// case. Infectious and symptomatic.
	InfStat_Case = -2,
	//// Recovered from asymptomatic infection
	InfStat_RecoveredFromAsymp = 3,
	//// Recovered from symptomatic infection
	InfStat_RecoveredFromSymp = -3,
	//// InfStat_Recovered (will use this for abs() values) so code reads correctly
	InfStat_Recovered = 3,
	//// Immune at start of epidemic - used to model partially immune population. Distinct therefore from recovered, who recovered during runtime. Doesn't take negative values.
	InfStat_ImmuneAtStart = 4,
	//// Dead was asymptomatic
	InfStat_Dead_WasAsymp = 5,
	//// Dead was symptomatic
	InfStat_Dead_WasSymp = -5,
	//// Dead (will use this for abs() values) so code reads correctly
	InfStat_Dead = 5
};

//// SeverityClass definitions / labels (numbers arbitrary but must be different to each other).
enum Severity {
	//// Flag value.
	Severity_Asymptomatic = -1,
	Severity_Mild = 0,
	Severity_ILI = 1,
	Severity_SARI = 2,
	Severity_Critical = 3,
	//// Recovering from Critical (not recovered yet).
	Severity_RecoveringFromCritical = 4,
	//// label to avoid double counting. Not sure you need.
	Severity_Dead = 5,
	//// label to avoid double counting. Not sure you need.
	Severity_Recovered = 6
};

#endif // COVIDSIM_INFSTAT_H_INCLUDED_
