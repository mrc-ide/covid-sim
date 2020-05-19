#include "Model.h"
#include <cstdlib>

bool Person::is_alive() const {
	return !this->is_dead();
}

bool Person::is_dead() const {
	return std::abs(this->infectionState) == InfStat_Dead;
}

void Person::make_susceptible() {
	this->infectionState = InfStat_Susceptible;
}

void Person::make_immune() {
	this->infectionState = InfStat_ImmuneAtStart;
}

void Person::make_infected() {
	this->infectionState = InfStat_Latent;
}

void Person::make_infectious_almost_symptomatic() {
	this->infectionState = InfStat_InfectiousAlmostSymptomatic;
}

void Person::make_infectious_asymptomatic() {
	this->infectionState = InfStat_InfectiousAsymptomaticNotCase;
}

void Person::make_case() {
	this->infectionState = InfStat_Case;
}

void Person::make_dead() {
	this->infectionState = (InfStat)(InfStat_Dead * this->infectionState / std::abs(this->infectionState));
}

void Person::make_recovered() {
	this->infectionState = (InfStat)(InfStat_Recovered * this->infectionState / std::abs(this->infectionState));
}

bool Person::is_susceptible() const {
	return this->infectionState == InfStat_Susceptible;
}

bool Person::is_latent() const {
	return this->infectionState == InfStat_Latent;
}

bool Person::is_symptomatic() const {
	return this->infectionState < InfStat_Susceptible;
}

bool Person::is_asymptomatic() const {
	return this->infectionState >= InfStat_Susceptible;
}

bool Person::is_case() const {
	return this->infectionState == InfStat_Case;
}

bool Person::is_exposed() const {
	return std::abs(this->infectionState) == InfStat_Latent;
}

bool Person::is_infectious() const {
	return this->infectionState == InfStat_InfectiousAsymptomaticNotCase
		|| this->infectionState == InfStat_Case;
}

bool Person::is_infected() const {
	return this->infectionState != InfStat_Susceptible
		&& this->infectionState != InfStat_ImmuneAtStart;
}

bool Person::has_recovered() const {
	return std::abs(this->infectionState) == InfStat_Recovered;
}

bool Person::has_not_recovered() const {
	return std::abs(this->infectionState) < InfStat_Recovered;
}

bool Person::is_positive() const {
	return this->is_infectious()
		|| this->infectionState == InfStat_InfectiousAlmostSymptomatic;
}