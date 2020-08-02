#include "InfStat.h"

bool not_yet_symptomatic(InfStat x)
{
	return (x == InfStat::Susceptible || x == InfStat::Latent || x == InfStat::InfectiousAlmostSymptomatic);
}

// Equivalent to forms like abs(Hosts[i].inf) == Dead
bool is_dead(InfStat x)
{
	return (x == InfStat::Dead || x == InfStat::Dead_WasSymp || x == InfStat::Dead_WasAsymp);
}

// Equivalent to forms like abs(Hosts[i].inf) == Recovered
bool is_recovered(InfStat x)
{
	return (x == InfStat::Recovered) || (x == InfStat::RecoveredFromSymp) || (x == InfStat::RecoveredFromAsymp);
}

// Equivalent to forms like abs(Hosts[i].inf]) < Recovered
bool pre_recovered(InfStat x)
{
	return	(x == InfStat::Susceptible) ||
		(x == InfStat::Latent) || (x == InfStat::InfectiousAlmostSymptomatic) ||
		(x == InfStat::InfectiousAsymptomaticNotCase) || (x == InfStat::Case);
}

bool never_symptomatic(InfStat x)
{
	return (x == InfStat::Latent) || (x == InfStat::InfectiousAsymptomaticNotCase) ||
		(x == InfStat::RecoveredFromAsymp) || (x == InfStat::ImmuneAtStart) ||
		(x == InfStat::Dead_WasAsymp);

}
