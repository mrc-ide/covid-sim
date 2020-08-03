#include "Models\Person.h"
#include "InfStat.h"

bool Person::is_susceptible() const
{
	return (this->inf == InfStat::Susceptible);
}

bool Person::is_immune_at_start() const
{
  return (this->inf == InfStat::ImmuneAtStart);
}

bool Person::is_case() const
{
	return (this->inf == InfStat::Case);
}

void Person::set_susceptible() {
	this->inf = InfStat::Susceptible;
}

bool Person::not_yet_symptomatic() const
{
	return (this->inf == InfStat::Susceptible || this->inf == InfStat::Latent || this->inf == InfStat::InfectiousAlmostSymptomatic);
}

// Equivalent to forms like abs(Hosts[i].inf) == Dead
bool Person::is_dead() const
{
	return (this->inf == InfStat::Dead_WasSymp || this->inf == InfStat::Dead_WasAsymp);
}

bool Person::is_alive() const {
	return !this->is_dead();
}

bool Person::is_recovered_symp() const {
	return this->inf == InfStat::RecoveredFromSymp;
}
bool Person::is_dead_symp() const {
	return this->inf == InfStat::Dead_WasSymp;
}

bool Person::is_dead_was_asymp() const {
	return this->inf == InfStat::Dead_WasAsymp;
}

bool Person::is_infectious_almost_symptomatic() const {
	return this->inf == InfStat::InfectiousAlmostSymptomatic;
}

bool Person::is_infectious_asymptomatic_not_case() const {
	return this->inf == InfStat::InfectiousAsymptomaticNotCase;
}

bool Person::do_not_vaccinate() const {

	// Originally: inf < InfStat::InfectiousAlmostSymptomatic) || (Hosts[ai].inf >= InfStat::Dead_WasAsymp)
	// ie, inf is either < -1, or >= 5. (ie, -2, -3, -5, or 5)
	// 5 or -5 is "dead". -2 is case. -3 is RecoveredFromSymp.

	return this->is_dead() || this->is_case() || this->is_recovered_symp();
}

// Equivalent to forms like abs(Hosts[i].inf) == Recovered
bool Person::is_recovered() const
{
	return (this->inf == InfStat::RecoveredFromSymp) || (this->inf == InfStat::RecoveredFromAsymp);
}

bool Person::is_latent() const
{
	return (this->inf == InfStat::Latent);
}

void Person::set_immune_at_start() {
	this->inf = InfStat::ImmuneAtStart;
}

void Person::set_latent() {
	this->inf = InfStat::Latent;
}

void Person::set_case() {
	this->inf = InfStat::Case;
}

void Person::set_recovered() {
	this->inf = (this->inf == InfStat::Case) ? InfStat::RecoveredFromSymp : InfStat::RecoveredFromAsymp;
}

void Person::set_dead() {
	this->inf = (this->inf == InfStat::Case) ? InfStat::Dead_WasSymp : InfStat::Dead_WasAsymp;
}

void Person::set_infectious_almost_symptomatic() {
	this->inf = InfStat::InfectiousAlmostSymptomatic;
}

void Person::set_infectious_asymptomatic_not_case() {
	this->inf = InfStat::InfectiousAsymptomaticNotCase;
}

// Equivalent to forms like abs(Hosts[i].inf]) < Recovered
bool Person::pre_recovered() const
{
	return	(this->inf == InfStat::Susceptible) ||
		(this->inf == InfStat::Latent) || (this->inf == InfStat::InfectiousAlmostSymptomatic) ||
		(this->inf == InfStat::InfectiousAsymptomaticNotCase) || (this->inf == InfStat::Case);
}

bool Person::never_symptomatic() const
{
	return (this->inf == InfStat::Latent) || (this->inf == InfStat::InfectiousAsymptomaticNotCase) ||
		(this->inf == InfStat::RecoveredFromAsymp) || (this->inf == InfStat::ImmuneAtStart) ||
		(this->inf == InfStat::Dead_WasAsymp);

}
