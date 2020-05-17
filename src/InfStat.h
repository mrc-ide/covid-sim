#ifndef COVIDSIM_INFSTAT_H_INCLUDED_
#define COVIDSIM_INFSTAT_H_INCLUDED_

//// Infection Status definitons / labels (generally positive value indicates asymptomatic infection, negative value indicates symptomatic infection).
enum InfStat {
	InfStat_Susceptible = 0,										//// Susceptible
	InfStat_Latent = 1,												  //// E or L (neither infectious nor symptomatic).
	InfStat_InfectiousAlmostSymptomatic =-1,		//// Infectious and about to become a case.
	InfStat_InfectiousAsymptomaticNotCase = 2,	//// Not just asymptomatic, but also will not become symptomatic (i.e. a case.)
	InfStat_Case = -2,													//// case. Infectious and symptomatic.
	InfStat_RecoveredFromAsymp = 3,							//// Recovered from asymptomatic infection
	InfStat_RecoveredFromSymp = -3,							//// Recovered from symptomatic infection
	InfStat_Recovered = 3,											//// InfStat_Recovered (will use this for abs() values) so code reads correctly
	InfStat_ImmuneAtStart = 4,									//// Immune at start of epidemic - used to model partially immune population. Distinct therefore from recovered, who recovered during runtime. Doesn't take negative values.
	InfStat_Dead_WasAsymp = 5,									//// Dead was asymptomatic
	InfStat_Dead_WasSymp = -5,									//// Dead was symptomatic
	InfStat_Dead = 5														//// Dead (will use this for abs() values) so code reads correctly
};

// SeverityClass defintions / labels (numbers arbitrary but must be different to each other).
enum Severity {
	Severity_Asymptomatic = -1,									//// Flag value.
	Severity_Mild = 0,
	Severity_ILI = 1,
	Severity_SARI = 2,
	Severity_Critical = 3,
	Severity_RecoveringFromCritical = 4,				//// Recovering from Critical (not recovered yet).
	Severity_Dead = 5,													//// label to avoid double counting. Not sure you need.
	Severity_Recovered = 6											//// label to avoid double counting. Not sure you need.
};

#endif // COVIDSIM_INFSTAT_H_INCLUDED_
